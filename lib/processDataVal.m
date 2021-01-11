function [cX, norm] = processDataVal(data_file, receptor, top_inds, responses)

d = load(data_file); 
Xfull = d.predictor;
Yfull = d.(['response_' receptor]);

nResps = length(responses);
cX = cell(1, nResps);
mu = cell(1, nResps);
sigma = cell(1, nResps);

for iY = 1:nResps
    % Get response data.
    response = strsplit(responses{iY}, '_');

    % Load data combination.
    X = Xfull(:, top_inds);

    % Filter data by binary for EC50 and EC10.
    if strcmp(response{1}(1), 'E')
        X = X(Yfull.BINARY == 1, :);
    end

    % Normalize columns.
    [cX{iY}, mu{iY}, sigma{iY}] = zscore(X);
end

norm.mu = mu;
norm.sigma = sigma;

end