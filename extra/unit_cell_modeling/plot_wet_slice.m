set(0,'defaulttextinterpreter','latex')
colmap = color_setup(4);
bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled

wet = readtable(fullfile('experiments','wet_20181119.csv'));

grid_ranges = struct('length',unique(wet.length), 'thick',unique(wet.thick), 'load',unique(wet.load), 'time',0:0.1:2);
fwt = strain_sample(grid_ranges);

lng = unique(wet.length);
figure; hold on; colormap(colmap);
xlabel('time, s');
ylabel('displacement, mm');
xlim([0 120]);
set(gca,'xtick',0:30:120);
set(gca,'ytick',0:1:5);
ax = gca; ax.FontSize = 8;
set(gca,'TickLabelInterpreter', 'latex');
grid on;
for i = 2%1:numel(lng)
    ylim([0 lng(i)]);
    wet_lng = wet(wet.length == lng(i),:);
    fwt_lng = fwt(fwt.length == lng(i),:);
    thk = unique(wet_lng.thick);
    for j = 1:numel(thk)
        wet_thk = wet_lng(wet_lng.thick == thk(j),:);
        fwt_thk = fwt_lng(fwt_lng.thick == thk(j),:);
        lod = unique(wet_thk.load);
        for k = 3%:numel(lod)
            title(sprintf('$l = %.f$ mm, $F = %.f$ N', lng(i) + bracket_length_correction, lod(k)));
            wet_lod = wet_thk(wet_thk.load == lod(k),:);
            fwt_lod = fwt_thk(fwt_thk.load == lod(k),:);
            wet_lod = downsample(wet_lod, 10);
%             ci = 1 + floor(min(1, (thk(j) - min(wet.thick))/(max(wet.thick) - min(wet.thick))) * (size(colmap,1) - 1));
            ci = 1 + round( (thk(j) - min(thk)) / (max(thk) - min(thk)) * (size(colmap,1)-1) );
            plot(wet_lod.time * 60, wet_lod.pos, 'Color',colmap(ci,:), 'LineWidth',2);
            plot(fwt_lod.time * 60, fwt_lod.pos, '--', 'Color',colmap(ci,:), 'LineWidth',0.5);
        end
    end
end
caxis([0.3 0.7]);
hcb = colorbar('Ticks',0.3:0.1:0.6);
colorTitleHandle = get(hcb,'Title');
colorTitleHandle.Interpreter = 'latex';
set(colorTitleHandle ,'String','$h$, mm');

hcb.TickLabelInterpreter = 'latex';

plotf_size(6, 5);
