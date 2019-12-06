function grid = make_grid(varargin)
    grid = cell(1,nargin);
    [grid{:}] = ndgrid(varargin{:});
    grid = cell2mat(cellfun(@(x) x(:),grid,'uniformoutput',false));
end
