function v = createVects(n)

v = zeros(2^n - 1,n);

for i = (2^n - 1):-1:1
    x = dec2bin(i);
    str = sprintf(['%0' num2str(n) 's'], x);
    v(i,:) = str - '0';
end

end