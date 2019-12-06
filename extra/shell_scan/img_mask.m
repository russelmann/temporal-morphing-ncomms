% black out all pixels not used by texture map (UV, TF)
function img = img_mask(imgname, UV, TF)
    imgin = imread(imgname);

    UV(:,1) = 1 + UV(:,1) * size(imgin,2);
    UV(:,2) = 1 + (1 - UV(:,2)) * size(imgin,1);

    TR2d = triangulation(TF,UV);

    [X,Y] = meshgrid(1:size(imgin,2),1:size(imgin,1));
    points = [X(:) Y(:)];
    pp = pointLocation(TR2d,points);
    pp = reshape(pp,size(imgin,1),size(imgin,2));
    pp = isnan(pp);

    img = imgin;
    pp = repmat(pp,1,1,size(imgin,3));
    img(pp) = 0;
end