set(0,'defaulttextinterpreter','latex')
colmap = color_setup(5);
bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled

wet = readtable(fullfile('experiments','wet_20181119.csv'));
wet.length = bracket_length_correction + wet.length;
wet.time = wet.time * 60;

lng = unique(wet.length);
figure; hold on; colormap(colmap);
xlabel('thickness, mm');
ylabel('$\log(dx/dt)$');
xlim([0.3 0.65]);
set(gca,'xtick',0.3:0.1:0.65);
ax = gca; ax.FontSize = 8;
set(gca,'TickLabelInterpreter', 'latex');
grid on;
for i = 1:numel(lng)
    wet_lng = wet(wet.length == lng(i),:);
    lod = unique(wet_lng.load);
    for j = 1:numel(lod)
        rates = zeros(numel(thk),1);
        wet_lod = wet_lng(wet_lng.load == lod(j),:);
        thk = unique(wet_lod.thick);
        for k = 1:numel(thk)
            wet_thk = wet_lod(wet_lod.thick == thk(k),:);
            wet_thk = downsample(wet_thk, 10);
            rates(k) = wet_thk.pos(end) / wet_thk.time(end);
        end
    end
    ci = 1 + round( (lng(i) - min(lng)) / (max(lng) - min(lng)) * (size(colmap,1)-1) );
    plot(thk, log(rates), 'Color',colmap(ci,:), 'LineWidth',2);
end

caxis([5 10]);
hcb = colorbar('Ticks',5:1:9);
colorTitleHandle = get(hcb,'Title');
colorTitleHandle.Interpreter = 'latex';
set(colorTitleHandle ,'String','$l$, mm');

hcb.TickLabelInterpreter = 'latex';

plotf_size(6, 5);
