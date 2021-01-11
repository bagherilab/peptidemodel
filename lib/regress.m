function s = regress(x, y, nComp)

% Calculate fits and r2 metric.
betas = zeros(size(x,2) + 1, nComp);
RSS = zeros(1, nComp);
nP = size(x,1);

for iComp = 1:min(nComp,size(x,2))
    [~,~,~,~,beta,~,MSE] = plsregress(x, y, iComp, 'cv', nP);
    [~,~,~,~,~,W,~,yVar] = plsnipals(x, y, iComp);
    betas(:,iComp) = beta;
    y_fit = [ones(size(x,1),1) x]*beta;
    RSS(iComp) = sum((y_fit - y).^2);
    
    % Catches cases where model has extremely large MSE.
    if MSE(2,iComp + 1) > 1e25 || isnan(MSE(2, iComp + 1)) || ...
            isnan(yVar(iComp)) || yVar(iComp) <= 1e-3
        % Delete unused columns from betas and RSS vectors.
        betas(:,(iComp + 1):end) = [];
        RSS((iComp + 1):end) = [];
        break;
    end
end

% Calculate cross validation metric.
PRESS = MSE(2,2:end)*nP;
TSS = sum((y - mean(y)).^2);

s.q2 = 1 - PRESS/TSS;
s.r2 = 1 - RSS/TSS;

% Keep betas and vip for all.
s.beta = betas;
s.vip = calcVip(W, yVar);
s.weights = W;

end