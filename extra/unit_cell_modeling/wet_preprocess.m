bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled

wet_name = 'wet_20181119';
vars = {'thick' 'length' 'load' 'pos' 'time'};
var_types = cell(numel(vars),1);
var_types(:,1) = {'double'};

load(fullfile('experiments', [wet_name '.mat']));
Wet = eval(wet_name);

for i = 1:numel(Wet)
    data = Wet{i}.data;
    
    % Convert time to minutes
    data.time = data.time / 60;
    
    % Cut tails
    d = (data.pos(2:end) - data.pos(1:end-1)) ./ (data.time(2:end) - data.time(1:end-1));
    tail_start = find(movmean(d, 10) < -0.1,1);
    if tail_start, data = data(1:tail_start,:); end
    
    % Cut extreme strains
    thlenlod = str2double(strsplit(Wet{i}.name,'_'));
    leng = thlenlod(2) - bracket_length_correction;
    data = data(data.pos / leng < 0.7,:);
    
    Wet{i}.data = data;
end

figure; hold on;
% xlim([0 11]);

% Pack into one large table
wet = table('Size',[sum(cellfun(@(x) size(x.data,1), Wet)) numel(vars)], 'VariableTypes',var_types, 'VariableNames',vars);

ind = 0;
for i = 1:numel(Wet)
    data = Wet{i}.data;
    
    tdata = table('Size',[size(data,1) numel(vars)], 'VariableTypes',var_types, 'VariableNames',vars);
    thlenlod = str2double(strsplit(Wet{i}.name,'_'));
    tdata.thick = ones(size(tdata,1),1) * thlenlod(1);
    tdata.length = ones(size(tdata,1),1) * thlenlod(2) - bracket_length_correction;
    tdata.load = ones(size(tdata,1),1) * thlenlod(3); %data.load; %N.B. we ignore error in load magnitude!
    tdata.pos = data.pos;
    tdata.time = data.time;
    wet(ind+1:ind+size(data,1),:) = tdata;
    ind = ind + size(data,1);
    
    plot3(tdata.time, tdata.pos, tdata.load);
end

writetable(wet, fullfile('experiments', [wet_name '.csv']));
