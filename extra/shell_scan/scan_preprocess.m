folder = 'data';

smp_shell = {};

% Specify shell name by commenting the rest %

smp_shell.sname = 'saddle';
% smp_shell.sname = 'scale';
% smp_shell.sname = 'spiral';
% smp_shell.sname = 'weave';

% Specific options for shells %

if strcmp(smp_shell.sname, 'saddle')
    smp_shell.whue = [0.5 0.6];
    smp_shell.wsat = [0.15 1];
    smp_shell.wval = [0.2 1];
end

if strcmp(smp_shell.sname, 'scale')
    smp_shell.whue = [0.1 1];
    smp_shell.wsat = [0.15 1];
    smp_shell.wval = [0.2 1];
end

if strcmp(smp_shell.sname, 'spiral')
    smp_shell.whue = [0.1 1];
    smp_shell.wsat = [0.15 1];
    smp_shell.wval = [0.2 1];
end

if strcmp(smp_shell.sname, 'weave')
    smp_shell.whue = [0.11 1];
    smp_shell.wsat = [0.15 1];
    smp_shell.wval = [0.1 1];
end

% Scan processing %

sname = smp_shell.sname;
whue = smp_shell.whue;
wsat = smp_shell.wsat;
wval = smp_shell.wval;

addpath('..');
colmap = color_setup(100);

[V,F,UV,TF,~,~] = readOBJ([folder '/' sname '_scan.obj']);

% Preprocess image (convert to binary for marker detection)
img_masked = img_mask([folder '/' sname '_scan.jpg'], UV, TF);
imwrite(img_masked,[folder '/' sname '_scan.png']);
img_marker = img_proc(img_masked, whue, wsat, wval);
imwrite(img_marker,[folder '/' sname '_scan.png']);
imshow(img_masked);
hold on
h = imshow(img_marker);
set(h,'AlphaData',0.5)

[V3s, W2s, W3s] = img_cloud(img_marker, V, F, UV, TF);
fid = fopen([folder '/' sname '_scan_markers.obj'],'w');
fprintf(fid, 'v %f\t%f\t%f\n', V3s');

figure
imshow(img_marker)
hold on
scatter(W2s(:,1), W2s(:,2))

figure; hold on;
scatter3(W3s(:,1), W3s(:,2), W3s(:,3))
scatter3(V3s(:,1), V3s(:,2), V3s(:,3))
axis equal
