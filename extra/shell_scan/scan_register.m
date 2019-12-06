folder = 'data';

smp_shell = {};

% Specify shell name by commenting the rest %

% smp_shell.sname = 'saddle';
smp_shell.sname = 'scale';
% smp_shell.sname = 'spiral';
% smp_shell.sname = 'weave'; % DO NOT USE, see scan_register_weave.m

sname = smp_shell.sname;

[V3s] = readOBJ([folder '/' sname '_scan_markers.obj']);

V3c = dlmread([folder '/' sname '_markers.txt']);

[tform, pairs, V3x] = reg_markers(V3c, V3s);

% Minimize distance for the matched pairs
V3z = V3c(pairs(:,2),:);
x0 = zeros(6,1);
x = fminunc(@(x) fdist(x, V3z, V3s), x0);
[~, V3s] = fdist(x, V3z, V3s);
dst = sqrt(sum((V3z - V3s).^2,2));

PD = pdist(V3c);
diameter = max(PD(:));
fprintf('Absolute mean error: %.2f\nRelative mean error: %.2f%%\nStandart deviation: %.2f\n', mean(dst), mean(dst) / diameter * 100, std(dst));

% Plot markers and distances
figure; hold on; axis equal
scatter3(V3c(:,1),V3c(:,2),V3c(:,3));
scatter3(V3x(:,1),V3x(:,2),V3x(:,3));
for i = 1:size(pairs,1)
    seg = [V3x(pairs(i,1),:); V3c(pairs(i,2),:)];
    plot3(seg(:,1), seg(:,2), seg(:,3),'b');
end

% Plot error histogram
% set(0,'defaulttextinterpreter','latex');
% figure
% set(gca,'TickLabelInterpreter', 'latex');
% hist(dst);
% xlabel('error, mm');
% width = 6;
% height = 5;
% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height], 'PaperUnits', 'centimeters', 'PaperSize', [width, height]);
% 
% h = findobj(gca,'Type','patch');
% h.FaceColor = [0 0.5 0.5];
% h.EdgeColor = 'w';

% Plot stencil with color-coded errors
% [V,F] = readOBJ([folder '/' 'saddle.obj']);
% V(:,3) = V(:,3) - min(V(:,3));
% TR = triangulation(F,V);
% C = -ones(size(F,1),1);
% C(pairs(:,2)) = dst;
% F0 = F(C < 0,:);
% F = F(C >= 0,:);
% C = C(C >= 0);
% figure; hold on;
% colormap(colmap);
% trisurf(F,V(:,1),V(:,2),V(:,3),C);
% h = trisurf(F0,V(:,1),V(:,2),V(:,3), F0(:,1)*0);
% h.FaceColor = [1 1 1] * 0.75;
% set(gca, 'xtick', [0 30]);
% set(gca, 'ytick', 0);
% set(gca, 'ztick', 0);
% ax = gca; ax.FontSize = 8; set(gca,'TickLabelInterpreter', 'latex');
% axis equal; grid on;

% hcb = colorbar;
% colorTitleHandle = get(hcb,'Title');
% colorTitleHandle.Interpreter = 'latex';
% set(colorTitleHandle ,'String','$\epsilon$, mm');
% hcb.TickLabelInterpreter = 'latex';
% hcb.Ticks = 0:0.5:1.5;
% 
% width = 10;
% height = 6;
% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, width, height], 'PaperUnits', 'centimeters', 'PaperSize', [width, height]);
