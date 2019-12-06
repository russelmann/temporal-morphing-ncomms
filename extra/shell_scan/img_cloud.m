% Lift markers from flat to 3D
% V3 joins close clusters of vertices from W3 (W3 returned for verification)
function [V3, W2, W3] = img_cloud(img, V, F, UV, TF)
    UV(:,1) = UV(:,1) * size(img,2);
    UV(:,2) = (1 - UV(:,2)) * size(img,1);

    CC = bwconncomp(1 - img);
    W2 = zeros(CC.NumObjects,2);
    for i = 1:CC.NumObjects
        [X,Y] = ind2sub(CC.ImageSize,CC.PixelIdxList{i});
        W2(i,:) = mean([Y X]);
    end
    
%     regions = detectMSERFeatures(img);
%     W22 = double(regions.Location);

    TR3d = triangulation(F, V);
    TR2d = triangulation(TF, UV);
    ID = pointLocation(TR2d, W2);
    W2 = W2(~isnan(ID),:);
    ID = ID(~isnan(ID));
    [ID, ia, ~] = unique(ID);
    W2 = W2(ia,:);

    W3 = barycentricToCartesian(TR3d,ID,cartesianToBarycentric(TR2d,ID,W2));
    
    % Join close vertices
    dst_join = 3;
    D = squareform(pdist(W3));
    D = D + eye(size(D)) * 100;
    [I, J] = ind2sub(size(D),find(D < dst_join));
    pairs = unique(sort([I J],2),'rows');
    clust = zeros(size(W3,1),1);
    ind = 1;
    for i = 1:size(pairs,1)
        if clust(pairs(i,1)) == 0
            if clust(pairs(i,2)) == 0
                clust(pairs(i,:)) = ind;
                ind = ind + 1;
            else
                clust(pairs(i,1)) = clust(pairs(i,2));
            end
        else
            if clust(pairs(i,2)) == 0
                clust(pairs(i,2)) = clust(pairs(i,1));
            else
                clust(clust == clust(pairs(i,2))) = clust(pairs(i,1));
            end
        end
    end
    for i = 1:numel(clust)
        if clust(i) == 0
            clust(i) = ind;
            ind = ind + 1;
        end
    end
    V3 = zeros(max(clust),3);
    cnt = zeros(max(clust),1);
    for i = 1:numel(clust)
        V3(clust(i),:) = V3(clust(i),:) + W3(i,:);
        cnt(clust(i)) = cnt(clust(i)) + 1;
    end
    V3 = V3 ./ cnt;
end