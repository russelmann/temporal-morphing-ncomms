function wet = strain_sample(grid_ranges)
    % Read dry and wet models and reconstruct final output
    load(fullfile('experiments','dry_fit.mat'));
    load(fullfile('experiments','wet_fit.mat'));

    vars = {'thick' 'length' 'load' 'time'};

    time = grid_ranges.time;
    
    % Resampling

    grid = make_grid(grid_ranges.thick, grid_ranges.length, grid_ranges.load);
    grid = sortrows(grid);
    Xdry = polyco(grid, dryfit.pcf);
    res_dry = Xdry * dryfit.beta;

    res_wet = cell(size(grid,1),1);

    for i = 1:size(grid,1)
        data = [ones(numel(time),3)  time'];
        data(:,1:3) = data(:,1:3) .* grid(i,:);
        Xwet = polyco(data, wetfit.pcf);
        ds = grid(i,3) .* exp(Xwet * wetfit.beta);
        dt = [time(2:end) - time(1:end-1)  -1]';
        res = res_dry(i) + cumsum([0; ds .* dt]);
        res(res > 10) = NaN;
        res_wet{i} = [data  res(1:end-1)];
    end
    wet = array2table(cell2mat(res_wet), 'VariableNames',[vars 'pos']);
    wet = wet(~isnan(wet.pos),:); % Filter large strains
    wet = wet(wet.pos < wet.length,:); % Filter large strains
    wet = wet(wet.load < 10,:); % Filter large stress
end

