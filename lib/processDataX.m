function cX = processDataX(data_file, receptor, grouping, responses, vectorizations, map)

d = load(data_file); 
Xfull = d.predictor;
Yfull = d.(['response_' receptor]);

nVects = length(vectorizations);
nResps = length(responses);
cX = cell(nVects, nResps);

for iY = 1:nResps
    fprintf('.');
    
    for iX = 1:nVects
        % Get predictor data.
        v = vectorizations(iX,:);
        full_inds = getFullInds(v, grouping, map);

        % Get response data.
        response = strsplit(responses{iY}, '_');
        
        % Load data combination.
        X = Xfull(:, full_inds);

        % Filter data by binary for EC50 and EC10.
        if strcmp(response{1}(1), 'E')
            X = X(Yfull.BINARY == 1, :);
        end
        
        % Normalize columns.
        cX{iX, iY} = zscore(X);
    end
end

end