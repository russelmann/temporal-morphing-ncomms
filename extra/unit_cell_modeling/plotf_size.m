% Resize figure to width x height in centimeters
function plotf_size(width, height)
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height], 'PaperUnits', 'centimeters', 'PaperSize', [width, height]);
end