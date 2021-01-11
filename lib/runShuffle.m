function s = runShuffle(data, nComp, nReps, vectorizations)

nVects = size(vectorizations, 1);

for iX = 1:nVects
    fprintf('     Running v %d\n', iX);
    v = vectorizations(iX,:);
    v_str = num2str(v);
    xid = ['x' v_str(~isspace(v_str))];
    x = data.predictor{iX};

    for iY = 1:nReps
        yid = ['y' num2str(iY)];
        y = data.response{iY};
        s.(yid).(xid) = regress(x, y, nComp);
    end
end

end