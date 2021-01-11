function vip = calcVip(W, pctvar)

[n, p] = size(W);
vin = zeros(n,p);
vip = zeros(n,p);

for i = 1:p
    w = W(:,i);
    pv = pctvar(i);
    vin(:,i) = (w.^2)*pv;
    vip(:,i) = sqrt(((n*sum(vin(:,1:i), 2))/sum(pctvar(1:i))));
end  

end