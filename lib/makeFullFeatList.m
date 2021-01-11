function [C, S] = makeFullFeatList() 

nProps = 8;
nRes = 13;
C = cell(1, nProps*nRes);
S.R = zeros(1, nProps*nRes);
S.P = zeros(1, nProps*nRes);
i = 1;

for iP = 1:nProps
    for iR = 1:nRes
        C{i} = ['residue ' num2str(iR) ' : property ' num2str(iP)];
        S.R(i) = iR;
        S.P(i) = iP;
        i = i+1;
    end
end

end