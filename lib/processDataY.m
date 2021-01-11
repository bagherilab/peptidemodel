function [cY, norm] = processDataY(data_file, receptor, responses)

d = load(data_file); 
Yfull = d.(['response_' receptor]);
nResps = length(responses);
cY = cell(nResps,1);
mu = zeros(1, nResps);
sigma = zeros(1, nResps);

for iY = 1:nResps
    % Get response data.
    response = strsplit(responses{iY}, '_');

    % Load data combination.
    Y = Yfull.(response{1});

    % Filter data by binary for EC50 and EC10.
    if strcmp(response{1}(1), 'E')
        Y = Y(Yfull.BINARY == 1);
    end

    % Try different variations on response data.
    if length(response) > 1
        if strcmp(response{2}, 'log')
            fun = @(val) log(val);
        elseif strcmp(response{2}, 'inv')
            fun = @(val) 1./val;
        end

        Y = fun(Y);
    end

    % Normalize columns.
    [cY{iY} mu(iY) sigma(iY)] = zscore(Y);
end

norm.mu = mu;
norm.sigma = sigma;

end