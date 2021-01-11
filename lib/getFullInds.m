function full_inds = getFullInds(vect, group, map)

switch group
    case 1
        full_inds = reshape(map(vect ~= 0, :)', 1, []);
    case 2
        full_inds = sort(reshape(map(:, vect ~= 0), 1, []));
end

end