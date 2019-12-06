function vars = polyco(x, pcf, varargin)
% generate polynomial
% x     data columns
% pcf   multiplicator indexes from
% varargin{1} > 0   derivative of the polynomial wrt the variable
% varargin{1} < 0   integral of the polynomial wrt the variable
    d = 0;
    if nargin > 2, d = 1 + abs(varargin{1}); end
    
    if nargin == 2
        q = 0;
    elseif nargin == 3
        q = sign(varargin{1});
    else
        q = varargin{2};
    end
    if q < -1, error('Double or more integration is not implemented'); end
    
    pol = [ones(size(x,1),1) x];
    vars = zeros(size(pol,1),size(pcf,1));
    nd = 1;
    for i = 1:size(pcf,1)
        pc = 1 + pcf(i,:);
        if d > 0 && q > 0
            nd = sum(pc==d);
            nd = prod(nd - (0:q-1));
            pc(find(pc==d,q)) = 1;
        end
        if d > 0 && q == -1
            nd = 1 / (sum(pc==d) + 1);
            pc = [pc d];
        end
        vars(:,i) = nd * prod(pol(:,pc),2);
    end
end