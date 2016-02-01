function [Fobs, p, dfb, dfe] = F_test_no_intercept(X,y,s)
% Test of the null hypothesis that the ENTIRE regression model X explains no
% variance in the data (y)
%
% :Usage:
% ::
%
%     [Fobs, p, dfb, dfe] = F_test_no_intercept(X,y,s)
%
% :Examples:
% Assess false positive rate with robust regression:
% ::
%
%    iter = 5000; 
%    warning off
%    fvals = zeros(1,iter); pvals = fvals;
%    for i = 1:iter
%      X = randn(10,1); y = randn(10,1);
%      [bb,stats]=robustfit(X,y);
%      % [bb,dev,stats] =glmfit(X,y);
%      [Fobs, p, dfb, dfe] = F_test_no_intercept(X,y,stats.s);
%      fvals(i) = Fobs; pvals(i) = p;
%    end
%    fpr = sum(pvals<.05) ./ iter
%    warning on
%
% :Based on:
% Johnson & Wichern, 5th ed., p. 371, Result 7.6
%
% for a regression model with k predictors, and q regression parameters 
% in reduced model, test the addition of q+1...k by extra sums of squares
% [n(s2red - s2full) / (k-q)]  / [n*s2full / (n-k-1)] ~ F(k-q),(n-k-1) 
%
% ..
%    WOrking function; in progress
%
%    i'm trying to generalize this to q = 0, when we care about the intercept
% ..

q = 0;                % number of params in reduced model; 0 for test of ENTIRE model, including intercept
interceptp = 0;       % number of intercept params; 0 for test of ENTIRE model, 1 for standard intercept

[n,k] = size(X);

if no_intercept(X)
    k = k + 1;
end

s2 = s ^ 2;
s2red = y' * y ./ n;   % full variance, zero parameters fit

varexp = s2red - s2;

dfb = k - q;            % + 1 because X doesn't contain intercept
dfe = n - k - interceptp;

Fnum = n * varexp / dfb;

Fden = n * s2 / dfe;

Fobs = Fnum / Fden;

p = 1 - fcdf(Fobs,dfb,dfe);

return





function val = no_intercept(X)

val = all(any(diff(X)));

return

