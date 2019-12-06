set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontSize',8)
colmap = color_setup(8);
bracket_length_correction = 2 * cos(37/180*pi) * 2; % base radius times bracket angle cosine, doubled

dry = readtable(fullfile('experiments','dry_20181102.csv'));
min_thk = min(dry.thick);
rng_thk = max(dry.thick) - min_thk;

vars = {'thick' 'length' 'load'};
grid_ranges = struct('thick',unique(dry.thick), 'length',unique(dry.length), 'load',0:0.1:10);
load(fullfile('experiments','dry_fit.mat'), 'dryfit');
grd = make_grid(grid_ranges.thick, grid_ranges.length, grid_ranges.load);
grd = sortrows(grd);
Xdry = polyco(grd, dryfit.pcf);
fdr = [grd Xdry * dryfit.beta];
fdr = array2table(fdr, 'VariableNames',[vars 'pos']);

lng = unique(dry.length);
figure; hold on; colormap(colmap);
for i = 1:numel(lng)
    dry_lng = dry(dry.length == lng(i),:);
    fdr_lng = fdr(fdr.length == lng(i),:);
    subplot(1,5,i); hold on; grid on;
    set(gca,'xtick',0:2:10);
    set(gca,'ytick',0:0.05:0.2);
    ax = gca; ax.FontSize = 8;
    set(gca,'TickLabelInterpreter', 'latex');
    title(sprintf('$l = %.f$ mm', lng(i) + bracket_length_correction));
    xlabel('load, N');
    if i == 1, ylabel('strain'); end
    xlim([0 10]);
    ylim([0 0.2]);
    thk = unique(dry_lng.thick);
    for k = 1:numel(thk)
        dry_thk = dry_lng(dry_lng.thick == thk(k),:);
        fdr_thk = fdr_lng(fdr_lng.thick == thk(k),:);
        dry_thk = downsample(dry_thk, 10);
        ci = 1 + round( (thk(k) - min_thk) / rng_thk * (size(colmap,1)-1) );
        plot(dry_thk.load, dry_thk.pos / lng(i), 'Color',colmap(ci,:), 'LineWidth',2);
        plot(fdr_thk.load, fdr_thk.pos / lng(i), '--', 'Color',colmap(ci,:), 'LineWidth',0.5);
    end
end
hp4 = get(subplot(1,5,5),'Position');
caxis([min_thk  min_thk + rng_thk + 0.05]);
hcb = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)+0.1  0.03  hp4(2)*5], 'Ticks', 0.3:0.05:0.65);
colorTitleHandle = get(hcb,'Title');
colorTitleHandle.Interpreter = 'latex';
set(colorTitleHandle ,'String','$h$, mm');
hcb.TickLabelInterpreter = 'latex';

plotf_size(18, 5);
