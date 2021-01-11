function s = parseResults(results_file, yids, nVec, cutoff, scaled_cutoff) 
 
n = (2^nVec - 1);
Q = zeros(length(yids), n);
VIP = cell(length(yids), n);
BETA = cell(length(yids), n);

% Step through each response type.
for iY = 1:length(yids)
    fprintf('.');
    yid = yids{iY};

    % Load data.
    d = load(results_file, yid);
    d = d.(yid);
    
    % Step through each vectorization by binary string index
    for iX = 1:n
        xid = sprintf(['x%0' num2str(nVec) 's'], dec2bin(iX)); % vectorization code
        [Q(iY,iX), ind] = max(d.(xid).q2); % get q2
        VIP{iY,iX} = d.(xid).vip(:, ind); % vip scores
        BETA{iY,iX} = d.(xid).beta(:, ind); % beta coefficients
    end
end

% Get best response type.
counts = sum(Q >= cutoff, 2); % count models with q2 above cutoff
[~, top_resp_ind] = max(counts); % top response type

% Get top values.
[top_q2_sorted, top_inds] = sort(Q(top_resp_ind, :), 'descend');
top_q2_scaled = top_q2_sorted/top_q2_sorted(1);
thresh = find(top_q2_scaled <= scaled_cutoff, 1) - 1;

% Store values in output structure.
s.q2 = Q;
s.q2_cutoff_counts = counts;
s.q2_cutoff_perc = counts/n; 
s.top_resp = top_resp_ind;
s.top_q2_sorted = top_q2_sorted;
s.top_vects_sorted = top_inds;
s.top_q2_scaled = top_q2_scaled;
s.q2_thresh = thresh;
s.VIPs = VIP(top_resp_ind, top_inds(1:thresh));
s.BETAs = BETA(top_resp_ind, top_inds(1:thresh));

end