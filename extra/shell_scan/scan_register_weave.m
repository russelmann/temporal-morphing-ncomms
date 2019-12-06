folder = 'data/weave';

sname = 'weave';

[V3c, Fc] = readOBJ([folder '/' sname '_markers_edg.obj']);
[V3s, Fs] = readOBJ([folder '/' sname '_scan_markers_edg.obj']);

[V3r, F] = readOBJ([folder '/' sname '_scan_markers_edg_r.obj']);
D = pdist2(V3c,V3r);
p = [1:size(V3s,1); munkres(D')]'; 
V3z = V3c(p(:,2),:);

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
V = V3s;
scatter3(V(:,1),V(:,2),V(:,3));
for i = 1:size(V3s,1)
    seg = [V(i,:); V3z(i,:)];
    plot3(seg(:,1), seg(:,2), seg(:,3),'b');
end
