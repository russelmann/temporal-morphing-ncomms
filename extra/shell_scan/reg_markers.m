function [tform, pairs, V3x] = reg_markers(V3c, V3s)
    pcs = pointCloud(V3s);
    pss = pointCloud(V3c);
    [tform,movingReg,~] = pcregistericp(pcs, pss);

    V3x = movingReg.Location;
    D = pdist2(V3c,V3x);
    ASSIGN = munkres(D');
    pairs = [(1:size(V3x,1)); ASSIGN]';
end
