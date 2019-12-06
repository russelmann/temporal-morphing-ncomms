bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled

dry_name = 'dry_20181102';

dry = readtable(fullfile('experiments', [dry_name '.csv']));
vars = {'thick' 'length' 'load'};

ithk = find(contains(vars,'thick'));
ilng = find(contains(vars,'length'));
ilod = find(contains(vars,'load'));

% Compute strain from displacement
dry.strain = dry.pos ./ dry.length;

% Find distinct values
thk = unique(dry.thick);
lng = unique(dry.length);

% Filter data
fltr = ismember(dry.thick, thk) & ismember(dry.length, lng);
dry = dry(fltr,:);

% Plot data-observations
figure; hold on;
xlim([0 11]);
scatter3(dry.load, dry.pos, dry.length, 3, dry.thick);

% Calculations
n = numel(vars); % number of variables
m = height(dry); % number of observations

% Define weights
W = eye(m);

% Construct polynomial
powers = cell2struct(cell(n,1),vars);
powers.thick = 1;
powers.length = 1;
powers.load = 2;
pcf = genpoly(cell2mat(struct2cell(powers))'); % thickness length load
pcf = [pcf  ones(size(pcf,1),1) * ilod]; % multiply by load

X = polyco(table2array(dry(:,vars)), pcf);
[~,jb] = rref(X);
pcf = pcf(jb,:); % remove linearly dependant polynomial terms (why are they?)
pcf = sort(pcf, 2);
X = polyco(table2array(dry(:,vars)), pcf);

% Extract observations
Y = dry.pos;

% Linear regression
% beta0 = (X' * W * X) \ (X' * W * Y); % LS weighted
% plot3(dry.load, X * beta0, dry.length, '.')

% Constrained regression
grid = cell2struct(cell(n,1),vars);
grid.thick = 0.3:0.05:0.65;
grid.length = (4:1:9) - bracket_length_correction;
grid.load = 0:1:10;
grid = struct2cell(grid);
[RTHK,RLNG,RLOD] = ndgrid(grid{:});
dryA = [RTHK(:) RLNG(:) RLOD(:)];
A = [
   -polyco(dryA, pcf, ithk);
    polyco(dryA, pcf, ilng);
    polyco(dryA, pcf, ilod);
%     -polyco(dryA, pcf, ilod, 2); % load second
];
options.Display = 'iter';
beta = quadprog(sparse(X' * W * X), -Y' * W * X, -A, A(:,1) * 0, [], [], [], [], [], options);
% plot3(dry.load, X * beta, dry.length, '.')

% Plotting
for i = 1:length(thk)
    for j = 1:length(lng)
        pfltr = (dry.thick == thk(i)) & (dry.length == lng(j));
        cth = (thk(i) - min(thk)) / (max(thk) - min(thk));
        col = [0  1-cth/2  1-cth/2];
        plot3(dry(pfltr,:).load, X(pfltr,:) * beta, dry(pfltr,:).length, '-', 'LineWidth',2 , 'Color',col)
    end
end

%Plot grid
grid = cell2struct(cell(n,1),vars);
grid.thick = 0.3:0.05:0.65;
grid.length = (4:0.05:9) - bracket_length_correction;
grid.load = 0:0.1:10;
grid = struct2cell(grid);
[mX, mY, mZ] = meshgrid(grid{:});
drygrid = array2table([mX(:) mY(:) mZ(:)], 'VariableNames',vars);
Xgrid = polyco(table2array(drygrid), pcf);
scatter3(drygrid.load, Xgrid * beta, drygrid.length, 1, drygrid.thick);

xlabel('load')
ylabel('displacement')
zlabel('length')

dryfit.pcf = pcf;
dryfit.beta = beta;
save(fullfile('experiments','dry_fit.mat'), 'dryfit');
