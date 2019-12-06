set(0,'defaulttextinterpreter','latex')
colmap = color_setup(5);

plastic_name = 'wet_speed_20181102';
vars = {'thick' 'length' 'load' 'pos' 'time'};
var_types = cell(numel(vars),1);
var_types(:,1) = {'double'};

load(fullfile('experiments', [plastic_name '.mat']));
Plastic = eval(plastic_name);

figure; hold on; grid on;
title(sprintf('$l = %.f$, $h = %.2f$ mm', 8, 0.4));
xlabel('time, s');
ylabel('displacement, mm');
set(gca,'xtick',0:30:120);
for i = 1:numel(Plastic)
    data = Plastic{i}.data;
    plot(data.time, data.pos, 'Color',colmap(i,:), 'LineWidth',2);
end
plot([0 120], 0.4 * [1 1], 'k--');
xlim([0 120]);
ylim([0 2.3]);
set(gca,'TickLabelInterpreter', 'latex');

plotf_size(6, 5);