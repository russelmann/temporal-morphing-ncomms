bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled

% Read dry and wet models and reconstruct final output
load(fullfile('experiments','dry_fit.mat'));
load(fullfile('experiments','wet_fit.mat'));

vars = {'thick' 'length' 'load' 'time'};

min_thick = 0.3;
step_thick = 0.05;
num_thick = 8;
min_length = 4 - bracket_length_correction;
step_length = 0.25;
num_length = 21;
min_time = 0;
step_time = 0.1;
num_time = 21;

brackets = struct('thick',min_thick + step_thick * (0:num_thick-1), ...
                  'length',min_length + step_length * (0:num_length-1));
time = min_time + step_time * (0:num_time-1);

% Resampling

grid = make_grid(brackets.thick, brackets.length, 0:0.01:10);
grid = sortrows(grid);
Xdry = polyco(grid, dryfit.pcf);
res_dry = Xdry * dryfit.beta;

res_wet = cell(size(grid,1),1);

for i = 1:size(grid,1)
    data = [ones(numel(time),3)  time'];
    data(:,1:3) = data(:,1:3) .* grid(i,:);
    Xwet = polyco(data, wetfit.pcf);
    ds = grid(i,3) .* exp(Xwet * wetfit.beta);
    dt = [time(2:end) - time(1:end-1)  -1]';
    res = res_dry(i) + cumsum([0; ds .* dt]);
    res(res > 10) = NaN;
    res_wet{i} = [data  res(1:end-1)];
end
wet = array2table(cell2mat(res_wet), 'VariableNames',[vars 'pos']);
wet = wet(~isnan(wet.pos),:); % Filter large strains
wet = wet(wet.pos < wet.length,:); % Filter large strains
wet = wet(wet.load < 10,:); % Filter large stress

% Fitting

% Define polynomial terms
ppow = 4;
pcf = fliplr(tril(ones(ppow)));

brgrid = sortrows(make_grid(brackets.thick, brackets.length, time));

Beta = zeros(ppow, size(brgrid,1));

for i = 1:size(brgrid,1)
    fltr = wet.thick == brgrid(i,1) & wet.length == brgrid(i,2) & wet.time == brgrid(i,3);
    fwet = wet(fltr,:);
    data = fwet.pos;
    
    if numel(data) == 1
        beta = zeros(ppow,1);
        wet.xload(fltr) = 0;
    else
        Y = fwet.load;

%         figure; hold on;
%         plot(fwet.pos, Y, 'x');

        % Construct polynomial
        X = data .^ (1:ppow);

    %     % Linear regression
    %     beta0 = (X' * X) \ (X' * Y); % LS
    %     plot3(fwet.pos, X * beta0, fwet.time, 'x')

        % Quadratic programming
%         options.Display = 'iter';
        options.Display = 'final';

        % Linear constraints for derivatives
        dataA = (0:0.01:brgrid(i,2))';
        A1 = polyco(dataA, pcf, 1); % displacement
        B1 = polyco(dataA, pcf, 1, 2); % displacement second
        A = [A1; -B1];
        beta = quadprog(sparse(X'*X), -Y'*X, -A, A(:,1) * 0, [], [], [], [], [], options);

%         plot(data(:,1), X * beta, 'o');
        
        wet.xload(fltr) = X * beta;
    end

    Beta(:,i) = beta;
end

writetable(wet, fullfile('experiments','stress_fit.csv'));

save(fullfile('experiments','stress_refit.mat'), 'brackets', 'time', 'Beta');

% Conversion to energy (analytic integration)
wBeta = Beta ./ (2:ppow+1)';

% output final energy polynominal
jout = {};
jout.bracket_model = {};

jout.bracket_model.grid_min = {};
jout.bracket_model.grid_min.rows = 3;
jout.bracket_model.grid_min.cols = 1;
jout.bracket_model.grid_min.data = [min_thick min_length min_time];

jout.bracket_model.grid_step = {};
jout.bracket_model.grid_step.rows = 3;
jout.bracket_model.grid_step.cols = 1;
jout.bracket_model.grid_step.data = [step_thick step_length step_time];

jout.bracket_model.grid_num = {};
jout.bracket_model.grid_num.rows = 3;
jout.bracket_model.grid_num.cols = 1;
jout.bracket_model.grid_num.data = [num_thick num_length num_time];

jout.bracket_model.beta = {}; %N.B. starts from power 2!
jout.bracket_model.beta.rows = size(wBeta,1);
jout.bracket_model.beta.cols = size(wBeta,2);
jout.bracket_model.beta.data = wBeta(:) / 4; % since there are 4 brackets per specimen
fout = fopen('bracket_model.json','w');
fprintf(fout, jsonencode(jout));
fclose(fout);
