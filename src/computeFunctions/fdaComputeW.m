function W = computeW(data, bList)

bNames = unique(bList);
k = length(bNames);
Nt = length(data);

Ni = arrayfun(@(x)sum(bList == x, bNames));



















