function X = hrf_saturation(X)
% X = hrf_saturation(X)
% nonlinear squashing function
% zero at x = 0, shrinks to -1 with higher values, elbow around 5

% u / (u+1)  = scaling factor from 0 to 1 for shrinkage
% u - u * [u / (u+1)] = shrink by scaling factor based on current level
% u - u^2 / (u + 1)

% satfun = @(u) (-u) ./ (u + 1);
% adjval = satfun(X);
% X = X + X .* adjval;
% does not work for blocks

% satfun = @(u) u - (u ./ (u + 1));
% %satfun = @(u) u - u.^2 ./ (u + 1);
% X = satfun(X);
%
% wh = intercept(X, 'which');
%
% X(:, wh) = 1;

% Simple bilinear saturation function for X
% hand-coded

for i = 1:size(X, 2)
    
    wh = X(:, i) > 3.5;
    X(wh, i) = X(wh, i) - .9 * (X(wh, i) - 3.5);
    
        wh = X(:, i) < -3.5;
    X(wh, i) = X(wh, i) - .9 * (X(wh, i) + 3.5);
    
    wh = X(:, i) > 2;
    X(wh, i) = X(wh, i) - .6 * (X(wh, i) - 2);
    
        wh = X(:, i) < 2;
    X(wh, i) = X(wh, i) - .6 * (X(wh, i) + 2);
    
    %X(X > 2, 1) = X(X > 2, 1) - .6 * (X(X > 2, 1) - 2);
    
end

end % function

