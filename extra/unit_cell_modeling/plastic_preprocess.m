bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled
% end_time = 5.25;

plastic_name = 'plastic_20181119';
vars = {'thick' 'length' 'load' 'pos' 'time'};
var_types = cell(numel(vars),1);
var_types(:,1) = {'double'};

load(fullfile('experiments', [plastic_name '.mat']));
Plastic = eval(plastic_name);

for i = 1:numel(Plastic)
    data = Plastic{i}.data;
    
    Plastic{i}.data = data;
end

figure; hold on;
% xlim([0 11]);

% Pack into one large table
plastic = table('Size',[sum(cellfun(@(x) size(x.data,1), Plastic)) numel(vars)], 'VariableTypes',var_types, 'VariableNames',vars);

pf = array2table(zeros(numel(Plastic),4));
pf.Properties.VariableNames = {'thick' 'length' 'maxpos' 'factor'};

ind = 0;
for i = 1:numel(Plastic)
    data = Plastic{i}.data;
    
    tdata = table('Size',[size(data,1) numel(vars)], 'VariableTypes',var_types, 'VariableNames',vars);
    thlenlod = str2double(strsplit(Plastic{i}.name,'_'));
    tdata.thick = ones(size(tdata,1),1) * thlenlod(1);
    tdata.length = ones(size(tdata,1),1) * thlenlod(2) - bracket_length_correction;
    tdata.load = data.load;
    tdata.pos = data.pos;
    tdata.time = data.time;
    plastic(ind+1:ind+size(data,1),:) = tdata;
    ind = ind + size(data,1);

    pf.thick(i) = tdata.thick(1);
    pf.length(i) = tdata.length(1);
    pf.maxpos(i) = max(tdata.pos);
    pf.factor(i) = tdata.pos(end) / max(tdata.pos);
    
    plot3(tdata.time, tdata.pos, tdata.load);
end

writetable(plastic, fullfile('experiments', [plastic_name '.csv']));

figure;
scatter3(pf.thick, pf.length, pf.factor, 10, pf.maxpos)
