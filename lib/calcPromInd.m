function I = calcPromInd(ec50s)

N = length(ec50s);

% Cap values at EC50 of 10,000.
inds = ec50s > 10000;
ec50s(inds) = 10000;
e = 1./ec50s;

% Calculate sums.
m = 0;
for i = 1:N
    denom = sum(e);
    m = m + (e(i)/denom)*log(e(i)/denom);
end

I = -1/log(N)*m;

end 