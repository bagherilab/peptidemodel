function V = normalize(v, range)
a = 0;
b = 1;
min_val = range(1);
max_val = range(2);
V = (b - a)*(v - min_val)/(max_val - min_val) + a;
end