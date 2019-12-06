clear singlelh
%singlelh = [2 2]; % plot single with length index singlelh(1) and thickness index singlelh(2); comment to plot all

set(0,'defaulttextinterpreter','latex')
colmap = color_setup(5);
bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled

wet = readtable(fullfile('experiments','wet_20181119.csv'));
min_load = min(wet.load);
rng_load = max(wet.load) - min_load + 1;

grid_ranges = struct('length',unique(wet.length), 'thick',unique(wet.thick), 'load',unique(wet.load), 'time',0:0.1:2);
fwt = strain_sample(grid_ranges);

lng = unique(wet.length);
if exist('singlelh', 'var'), lng = lng(singlelh(1)); end %%% single plot
figure; hold on; colormap(colmap);
for i = 1:numel(lng)
    wet_lng = wet(wet.length == lng(i),:);
    fwt_lng = fwt(fwt.length == lng(i),:);
    thk = unique(wet_lng.thick);
    if exist('singlelh', 'var'), thk = thk(singlelh(2)); end %%% single plot
    for j = 1:numel(thk)
        wet_thk = wet_lng(wet_lng.thick == thk(j),:);
        fwt_thk = fwt_lng(fwt_lng.thick == thk(j),:);
        subplot(numel(lng),numel(thk),j+(i-1)*numel(thk)); hold on; grid on;
        set(gca,'xtick',0:30:120);
        set(gca,'ytick',0:1:5);
        ax = gca; ax.FontSize = 8;
        set(gca,'TickLabelInterpreter', 'latex');
        title(sprintf('$l = %.f$, $h = %.2f$ mm', lng(i) + bracket_length_correction, thk(j)));
        if (i == numel(lng)), xlabel('time, s'); end
        if j == 1, ylabel('displacement, mm'); end
        xlim([0 120]);
        ylim([0 lng(i)]);
        lod = unique(wet_thk.load);
        for k = 1:numel(lod)
            wet_lod = wet_thk(wet_thk.load == lod(k),:);
            fwt_lod = fwt_thk(fwt_thk.load == lod(k),:);
            wet_lod = downsample(wet_lod, 10);
            ci = lod(k);
            plot(wet_lod.time * 60, wet_lod.pos, 'Color',colmap(ci,:), 'LineWidth',2);
            plot(fwt_lod.time * 60, fwt_lod.pos, '--', 'Color',colmap(ci,:), 'LineWidth',0.5);
        end
    end
end
if exist('singlelh', 'var')
    hcb = colorbar('Ticks',1:5);
else
    hp4 = get(subplot(numel(lng),numel(thk),numel(lng)*numel(thk)),'Position');
    hcb = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)  0.02  hp4(2)], 'Ticks',1:5);
end
caxis([1 6]);
colorTitleHandle = get(hcb,'Title');
colorTitleHandle.Interpreter = 'latex';
set(colorTitleHandle ,'String','load, N');
hcb.TickLabelInterpreter = 'latex';

plotf_size(18, 26);
if exist('singlelh', 'var'), plotf_size(6, 5); end %%% single plot
