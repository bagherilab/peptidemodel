function s = processPredictors(data_file, nVect, dict_file)

% Parse data from excel file.
[~, txt] = xlsread(data_file);

% Get sequences and vectorize.
s.keys = txt(2:end,1);
s.seqs = txt(2:end,2);

v = [];
for i = 1:length(s.seqs)
    v = [v; vectorize(s.seqs{i}, ones(1, nVect), dict_file)];
end

s.predictor = v;

end