function single_vect = calcSingle(seq, single_dict)

% Iterate through entire sequence and find single residue values.
single_vect = zeros(1,length(seq));
for iAA = 1:length(seq)
    key = seq(iAA);
    
    if isKey(single_dict, key)
        single_vect(iAA) = single_dict(key);
    else
        single_vect(iAA) = 0;
    end
end

end