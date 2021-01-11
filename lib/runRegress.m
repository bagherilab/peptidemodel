function s = runRegress(data, nComp, vectorizations, responses)

nVects = size(vectorizations, 1);
nResps = size(responses, 2);

for iX = 1:nVects
    fprintf('     Running v %d\n', iX);
    v = vectorizations(iX,:);
    v_str = num2str(v);
    xid = ['x' v_str(~isspace(v_str))];

    for iY = 1:nResps
        yid = responses{iY};
            
        % Run regression.
        x = data.predictor{iX, iY};
        y = data.response{iY};
        s.(yid).(xid) = regress(x, y, nComp);
    end
end
    
end