function V = vectorize(seq, vect_id, dict_file)

d = load(dict_file);
V = [];

for iFeat = 1:length(vect_id)
    dict = d.(d.single_labels{iFeat});
    V = [V calcSingle(seq, dict)];
end

end