function img = img_proc(imgin, whue, wsat, wval)
    % Transform to HSV space
    imhsv = rgb2hsv(imgin);

    hue = imhsv(:,:,1);
    sat = imhsv(:,:,2);
    val = imhsv(:,:,3);

    % Filter marker color
    hue(whue(1) <= hue & hue <= whue(2)) = 1;
    sat(wsat(1) <= sat & sat <= wsat(2)) = 1; 
    val(wval(1) <= val & val <= wval(2)) = 1;

    hue(hue < 1) = 0;
    sat(sat < 1) = 0;
    val(val < 1) = 0;

    img = imhsv .* 0;
    img(:,:,1) = hue;
    img(:,:,2) = sat;
    img(:,:,3) = val;
    img = prod(img,3);
    img = imfill(img,'holes');

    img = 1 - bwareaopen(img,10);
end
