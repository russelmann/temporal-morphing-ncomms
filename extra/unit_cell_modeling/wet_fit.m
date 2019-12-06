bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled

wet = readtable(fullfile('experiments','wet_20181119.csv'));
vars = {'thick' 'length' 'load' 'time'};

ithk = find(contains(vars,'thick'));
ilng = find(contains(vars,'length'));
ilod = find(contains(vars,'load'));
itme = find(contains(vars,'time'));

% Compute strain from displacement
wet.strain = wet.pos ./ wet.length;

% Find distinct values
thk = unique(wet.thick);
lng = unique(wet.length);
lod = unique(wet.load);

% Filter data
fltr = ismember(wet.thick, thk) & ismember(wet.length, lng);
wet = wet(fltr,:);

% Compute pos rates
wet = sortrows(wet, {'thick' 'length' 'load' 'time'});
wet.dt = [wet.time(2:end) - wet.time(1:end-1); -1];
ds = [wet.pos(2:end) - wet.pos(1:end-1); 0];
wet.strate = log(max(1e-2,ds ./ wet.dt) ./ wet.load);
wet = wet(wet.dt > 0,:); % remove indefinite pos rates

% Plot data-observations
% fig_pos = figure; hold on;
% xlabel('time,min')
% ylabel('pos')
% zlabel('stress,N')
% xlim([0 2]);
% ylim([0 5]);
% plot3(wet.time, wet.pos, wet.load, '.')

% Calculations
n = length(vars); % number of variables
m = height(wet); % number of observations

% Define weights
[~,~,cid] = unique(wet(:,{'thick' 'length' 'load'}),'rows');
nc = sum(cid == 1:max(cid))';
wc = min(nc) ./ nc;
W = diag(sparse(wc(cid)));
% W = eye(m);

% Plot different data slices
% for i = 1:length(thk)
%     for j = 1:length(lod)
%         subplot(3,3,j+i*3-3);
%         fltr = (wet(:,2) == thk(i)) & (wet(:,4) == lod(j));
%         plot(wet(fltr,1), Y(fltr), '.');
%         xlim([0 2]);
%         ylim([0 1]);
%     end
% end

% Plot data-observations
figure; hold on;
xlabel('time,min')
ylabel('log-pos rate')
zlabel('stress,N')
xlim([0 2]);
plot3(wet.time, wet.strate, wet.load, '.')

% Construct polynomial
pcf = genpoly([4 4 4 4]);
X = polyco(table2array(wet(:,vars)), pcf);
[~,jb] = rref(X);
pcf = pcf(jb,:); % remove linearly dependent polynomial terms (why are they?)
pcf = sort(pcf, 2);
X = polyco(table2array(wet(:,vars)), pcf);

% Extract observations
Y = wet.strate;

% Linear regression
% beta0 = (X' * W * X) \ (X' * W * Y); % LS weighted
% plot3(wet.time, X * beta0, wet.load, '.')

% Quadratic programming
options.Display = 'iter';

% linear constraints for second derivatives
grid = cell2struct(cell(n,1),vars);
grid.thick = 0.3:0.05:0.65;
grid.length = (4:1:9) - bracket_length_correction;
grid.load = 0:1:10;
grid.time = 0:0.1:2;
grid = struct2cell(grid);
[RTHK,RLNG,RLOD,RTME] = ndgrid(grid{:});
wetA = [RTHK(:) RLNG(:) RLOD(:) RTME(:)];

% wetA = table2array(wet(:,vars));
% wetA(:,4) = wetA(:,4) + 5;
% wetA = [table2array(wet(:,vars)); wetA];
A1 = polyco(wetA, pcf, ithk); % thick
A2 = polyco(wetA, pcf, ilng); % length
A3 = polyco(wetA, pcf, ilod); % load
A4 = polyco(wetA, pcf, itme); % time
A = [-A1; A2; A3; A4];
% A = [A1; -A2; A4];

beta = quadprog(sparse(X'*W*X), -Y'*W*X, -A, A(:,1) * 0, [], [], [], [], [], options);
plot3(wet.time, X * beta, wet.load, '.')

% Output fitted polynomial
wetfit.pcf = pcf;
wetfit.beta = beta;
save(fullfile('experiments','wet_fit.mat'), 'wetfit');
