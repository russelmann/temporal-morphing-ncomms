set(0,'defaulttextinterpreter','latex')
colmap = color_setup(5);

data = table2cell(readtable('rubber_stress_relaxation.csv'));
data = cell2mat(data(:,1:end-1));
figure; hold on; grid on;
% title(sprintf('$l = %.f$, $h = %.2f$ mm', 8, 0.4));
xlabel('time, h');
ylabel('force, N');
%set(gca,'xtick',0:3:24);
plot(data(:,1)/60, data(:,8), 'Color',colmap(3,:), 'LineWidth',2);
xlim([0 60]);
ylim([0 5]);
set(gca,'TickLabelInterpreter', 'latex');

plotf_size(6, 5);