set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontSize',8)
colmap = color_setup(8);

% test = struct('name','reg_bend_20191029.xls' , 'title','bending', 'xlim',25);
% test = struct('name','reg_compress_20191029.xls' , 'title','compression', 'xlim',25);
test = struct('name','reg_shear_20191029.xls' , 'title','shearing', 'xlim',10);
% test = struct('name','reg_stretch_20191029.xls' , 'title','stretching', 'xlim',50);

data = readtable(fullfile('experiments',test.name), 'Sheet','Test 1');
data = table2array(data(3:end,:));
data = cellfun(@str2num, data);
data = downsample(data, 100);

figure; hold on; grid on; colormap(colmap);
ax = gca; ax.FontSize = 8;
set(gca,'TickLabelInterpreter', 'latex');
title(test.title);
xlabel('strain, \%');
ylabel('stress, MPa');
xlim([0 test.xlim]);
%ylim([0 inf]);

plot(data(:,1), data(:,2), '-', 'Color',colmap(3,:), 'LineWidth',2);

plotf_size(6, 5);
