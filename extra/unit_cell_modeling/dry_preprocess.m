bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled
end_time = 5.25;

dry_name = 'dry_20181102';
vars = {'thick' 'length' 'load' 'pos'};
var_types = cell(numel(vars),1);
var_types(:,1) = {'double'};

load(fullfile('experiments', [dry_name '.mat']));
Dry = eval(dry_name);

for i = 1:numel(Dry)
    data = Dry{i}.data;
    
    % Cut tails
    data = data(data.time < end_time, :);
    
    Dry{i}.data = data;
end

figure; hold on;
xlim([0 11]);

% Pack into one large table
dry = table('Size',[sum(cellfun(@(x) size(x.data,1), Dry)) numel(vars)], 'VariableTypes',var_types, 'VariableNames',vars);

ind = 0;
for i = 1:numel(Dry)
    data = Dry{i}.data;
    
    tdata = table('Size',[size(data,1) numel(vars)], 'VariableTypes',var_types, 'VariableNames',vars);
    thlen = str2double(strsplit(Dry{i}.name,'_'));
    tdata.thick = ones(size(tdata,1),1) * thlen(1);
    tdata.length = ones(size(tdata,1),1) * thlen(2) - bracket_length_correction;
    tdata.load = data.load;
    tdata.pos = data.pos;
    dry(ind+1:ind+size(data,1),:) = tdata;
    ind = ind + size(data,1);
    
    plot3(tdata.load, tdata.pos, tdata.length);
end

writetable(dry, fullfile('experiments', [dry_name '.csv']));
