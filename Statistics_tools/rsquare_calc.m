function rsquare = rsquare_calc(X, Y)
%
% rsquare(X, y)
%
% Returns variance in each col. of Y explained by X

% make sure X has intercept
X = intercept(X, 'end');

for i = 1:size(Y, 2)

    y = Y(:,i);

    b = X\y;

    fits = X * b;

    rsquare(i) = var(fits) / var(y);

end

end