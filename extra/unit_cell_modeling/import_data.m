dry_names = {'dry_20181102'};
dry_export_name = 'dry_20181102';
wet_names = {'wet_1N_20181119' 'wet_2N_20181119' 'wet_3N_20181119' 'wet_4and5N_20181119'};
wet_export_name = 'wet_20181119';
plastic_names = {'plastic_20181119'};
plastic_export_name = 'plastic_20181119';
wet_speed_names = {'wet_speed_20181102'};
wet_speed_export_name = 'wet_speed_20181102';

vars = {'time' 'pos' 'load'};

import_xls(dry_names, dry_export_name, vars);
import_xls(wet_names, wet_export_name, vars);
import_xls(plastic_names, plastic_export_name, vars);
import_xls(wet_speed_names, wet_speed_export_name, vars);

function import_xls(test_names, export_name, vars)
    numex = zeros(numel(test_names),1);
    for i = 1:numel(test_names)
        fname = fullfile('experiments', [test_names{i} '.xls']);
        [~,sheet_names] = xlsfinfo(fname);
        numex(i) = numel(sheet_names);
    end
    imp_cell = cell(sum(numex),1);
    ind = 1;
    for test_name = test_names
        fname = fullfile('experiments', [test_name{1} '.xls']);
        [~,sheet_names] = xlsfinfo(fname);
        for k = 1:numex(i)
            data_mat = xlsread(fname, sheet_names{k});
            imp_cell{ind}.name = sheet_names{k};
            imp_cell{ind}.data = array2table(data_mat, 'VariableNames',vars);
            ind = ind + 1;
        end
    end
    save_cell.(export_name) = imp_cell;
    save(fullfile('experiments', [export_name '.mat']), '-struct', 'save_cell');
end