function s = processDicts(dict_file)

%% Parse out single residue vectors.
[num, txt] = xlsread(dict_file, 'VHSE'); 
for iFeat = 1:8
   % Add feature name and min/max of values.
    feat_name = txt{1, iFeat + 1};
    s.single_labels{iFeat} = feat_name;
    range = [min(num(:,iFeat)) max(num(:,iFeat))];

    % Iterate through each residue and add to dictionary.
    dict = containers.Map;
    for iResidue = 1:20
        dict(txt{iResidue + 1,1}) = normalize(num(iResidue, iFeat), range);
    end

    s.(feat_name) = dict;
end

end