function [E, V] = fdist(x, V3z, V3s)
    R = eye(3);
    if 0 < norm(x(4:6))
        R = vrrotvec2mat([x(4:6) / norm(x(4:6)); norm(x(4:6))]);
    end
    V = V3s * R + ones(size(V3s,1),1) * x(1:3)';
    E = sum(sqrt(sum((V3z - V).^2,2)));
end