function y = mediansplit(y)

wh = find(y > median(y));

y = -1 * ones(size(y));
y(wh) = 1;

% sometimes this can return all -1s or all 1s


return
