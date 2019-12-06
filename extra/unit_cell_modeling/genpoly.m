function pcf = genpoly(pows)
% Generate polynomial coefficients up to given powers
% pows   row vector of maximal powers for each variable
    n = size(pows,2);
    maxpow = max(pows);
    R = cell(maxpow,1);
    [R{:}] = ndgrid(0:n);
    R = reshape(cell2mat(R),(1+n)*maxpow,[]);
    R = reshape(R',[],maxpow);
    pcf = unique(sort(R,2),'rows');

    fltr = true(size(pcf,1),1);
    for i = 1:n
        fltr = fltr & (sum(pcf==i,2) <= pows(i));
    end
    pcf = pcf(fltr,:);
end

