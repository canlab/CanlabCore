function [t2, S, p, t2crit, dfb, dfe, sc] = tsquaretest(X,pthresh, u0)
% Hotelling's t-square test that the multivariate sample mean of X is
% different from u0.
%
% X is an observations x variables data matrix
% If pthresh is entered, the critical t2 value is returned
% If u0 is not entered, the default null hypothesis is mean zero (the
% origin)
%
% ..
%    tor wager, july 7, 2006
% ..

[n,k] = size(X);

if nargin < 3
    % null hypothesis test value
    u0 = zeros(1,k);
end

% mean
xbar = mean(X);

% deviations
Xdev = X - repmat(xbar,n,1);

% sample covariance (follows a Wishart dist)
S = Xdev' * Xdev ./ (n-1);

% Hotelling's t-square
t2 = n * (xbar - u0) * inv(S) * (xbar - u0)';


% scaling of F dist; t2 follows sc*F(k,n-k)
sc = k * (n-1) ./ (n - k);

% degrees of freedom
dfb = k; dfe = n - k;

if nargin > 1
    % critical t2 value
    t2crit = sc * finv(1 - pthresh,dfb,dfe);
else
    t2crit = [];
end

% p-value
p = 1 - fcdf(t2 ./ sc,dfb,dfe);

return

