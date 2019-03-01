function [d_CI, P, P_CI] = effect_size_CI(d, n, varargin)
% Confidence interval for effect size (Cohen's d or Hedges' g) using noncentral t distribution
%
% d_CI = effect_size_CI(d, n, [dfe])
%
% - This function assumes dfe is n-1 unless specified otherwise
% - Assumes normally distributed data
%
% Inputs:
% d = Effect size estimate
% n = Sample size
%
% Outputs:
% d_CI = 95% confidence interval for d, 2-tailed
% P = P-value for the input d
% P_CI = 95% confidence interval for P, 2-tailed
%
% There are multiple methods for this and formulas are suggested by
% Hedges (1981), Theory for Glass's Estimator of Effect Size and Related Estimators
%
% This function uses the method suggested by David C. Howell, University of Vermont
% Using the CDF of the non-central t distribution for the observed t score
%
% See also:
% Hedges LV, Olkin I. Statistical methods for meta-analysis. Orlando: Academic Press Inc; 2014. p. 86.
% https://www.statisticshowto.datasciencecentral.com/hedges-g/
% https://books.google.com/books?id=5obZnfK5pbsC&pg=PA10&dq=hedges+g&hl=en&sa=X&ved=0ahUKEwjLl6eL4NPRAhVowYMKHQzGDpsQ6AEIGjAA#v=onepage&q=hedges%20g&f=false
%
% Formulas
% p2t = @(p, dfe, tails) abs(tinv(p ./ tails, dfe));
% p2d = @(p, dfmodel, dfe, tails) abs(tinv(p ./ tails, dfe)) ./ sqrt(dfmodel + dfe);
% d2t = @(d, dfmodel, dfe) d .* sqrt(dfmodel+dfe);
% d2p = @(d, dfmodel, dfe, tails) tails .* (1 - tcdf(d .* sqrt(dfmodel+dfe), dfe));
%
% Examples:
% ---------------------------------------
% e.g.,
% n = 100; d = 0.65;
% d_CI = effect_size_CI(d, n);
%
% By Tor Wager, 10/2018, alpha version. Comparisons with other methods and
% published examples should be performed, independent checks should be
% performed.

% t  = d * sqrt(N)
% we want critical value x where noncentral tcdf is 0.025 and 0.975 for 2-tailed CI
% or .05 .95 for 90% CI

tails = 2;

if length(varargin) > 0
    dfe = varargin{1};
else
    dfe = n - 1;
end

t = d * sqrt(n);

% Confidence interval on t

x = [0:.001:t];           % we want critical value x where noncentral tcdf is 0.025 and 0.975 for 2-tailed CI
p = nctcdf(t, dfe, x);  % Non-central t distribution CDF

wh = find(p < 0.975);   % 2.5% in each tail
t_CI = x(wh(1));    % bounds on t (noncentrality parameter)

x = [t:.001:4*t];  % we want critical value x where noncentral tcdf is 0.025 and 0.975 for 2-tailed CI
p = nctcdf(t, dfe, x);
wh = find(p > 0.025);

t_CI(2) = x(wh(end));

d_CI = t_CI ./ sqrt(n); % convert back to units of d

% Now get P-value and Confidence interval on P

% Convert d back to p
d2p = @(d, n, dfe, tails) tails .* (1 - tcdf(d .* sqrt(n), dfe));

P = d2p(d, n, dfe, tails);
P_CI = d2p(d_CI, n, dfe, tails);

end % main function


