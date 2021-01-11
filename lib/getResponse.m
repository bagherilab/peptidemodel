function val = getResponse(seq, iRecep, iResp, receptor, response, thresh, dict_file, val_file)
%% Calculate the selected response for the selected receptor given a new
%% peptide sequence.

% Check sequence length and take only the first 13 if longer.
if length(seq) > 13
    seq = seq(1:13);
end

% Get selected top features.
d = load(val_file);
top_inds = d.top_feat_inds{d.thresholds == thresh, iRecep};

% Load selected regression
plsr = d.(['Results_' receptor]).(response);
[~, ind] = max(plsr.q2);
betas = plsr.beta(:,ind);

% Get normalization parameters.
data = d.(['Data_' receptor]);
mux = data.Xnorm.mu{iResp};
muy = data.Ynorm.mu(iResp);
sigmax = data.Xnorm.sigma{iResp};
sigmay = data.Ynorm.sigma(iResp);

% Get selected indices for predictors.
fullVect = vectorize(seq, ones(1, 8), dict_file);
X = fullVect(top_inds);
x = (X - mux)./sigmax;

% Calculate value.
val = x*betas(2:end)*sigmay + muy;

end