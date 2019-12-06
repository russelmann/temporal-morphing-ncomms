function colmap = color_setup(numcol)
    set(0,'defaulttextinterpreter','latex')
    colorMap = cbrewer('seq','GnBu',(numcol * 2) - 1);
    colmap = colorMap(numcol-2 + (1:numcol),:);
end

