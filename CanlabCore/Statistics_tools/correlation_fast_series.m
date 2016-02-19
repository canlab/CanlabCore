function [r, p, Tstat] = correlation_fast_series(Xi, Yi)
% :Usage:
% ::
%
%     [r, p, t] = correlation_fast_series(Xi, Yi)  
%
% Fast, memory-efficient way to get correlation between each column of Xi and Yi
%
%
% :Inputs:
%
%   **Xi:**
%        Matrix of observations (n instances by p variables)
%
%   **Yi:**
%        Matrix of observations (n instances by p variables)
%
%
% :Outputs:
%
%   **r:**
%        Pearson correlation coefficients
%
%   **p:**
%        Corresponding p-value for correlation coeffieints
%
%   **Tstat:**
%        T statistic for comparing correlation coefficients against 0
%
%
% :Examples:
% ::
%
%    r = ((Xi' * Yi) ./ N) ./ sqrt(diag(Xi' * Xi) ./ N * (Yi' * Yi) ./ N);
%    [r2, p, Tstat] = correlation_fast_series(Xi, Yi);
%    [r3, p2] = corrcoef([Yi, Xi]); r3 = r3(1, 2:end); p2 = p2(1, 2:end);
%    create_figure('test', 1, 3); plot(r3, r, 'kx'); axis equal; grid on; subplot(1, 3, 2); plot(r3, r2, 'kx');  axis equal; grid on; subplot(1, 3, 3); plot(p2, p, 'kx'); axis equal; grid on;
%
% ..
%    Tor Wager, June 2010
%
%    Notes: tested against corrcoef.m on June 28, 2010, tor
% ..

[N, m] = size(Xi); 

vx = var(Xi, 0)';
vy = var(Yi, 0);
Xi = scale(Xi, 1); Yi = scale(Yi, 1);
r = ((Xi' * Yi) ./ N) ./ sqrt(vx * vy);

% This works, but is too memory intensive because of Xi*Xi
% also, needs mean-centering??
%r = ((Xi' * Yi) ./ N) ./ sqrt(diag(Xi' * Xi) ./ N * (Yi' * Yi) ./ N);

% adjust for perfect correlations; avoid warnings
% done in corrcoef.m
r(r == 1) = 1 - 100*eps;
r(r == -1) = -1 + 100*eps;

Tstat = r .* sqrt((N-2) ./ (1 - r.^2));
   p = zeros(m, 1, class(Xi));
   p = 2 * tpvalue(-abs(Tstat), N - 2); % two-tailed p-value
   
   
   
% Borrowed from corrcoef.m, Matlab 7.5   
function p = tpvalue(x,v)
%TPVALUE Compute p-value for t statistic

normcutoff = 1e7;
if length(x)~=1 && length(v)==1
   v = repmat(v,size(x));
end

% Initialize P.
p = NaN(size(x));
nans = (isnan(x) | ~(0<v)); %  v==NaN ==> (0<v)==false

% First compute F(-|x|).
%
% Cauchy distribution.  See Devroye pages 29 and 450.
cauchy = (v == 1);
p(cauchy) = .5 + atan(x(cauchy))/pi;

% Normal Approximation.
normal = (v > normcutoff);
p(normal) = 0.5 * erfc(-x(normal) ./ sqrt(2));

% See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1
gen = ~(cauchy | normal | nans);
p(gen) = betainc(v(gen) ./ (v(gen) + x(gen).^2), v(gen)/2, 0.5)/2;

% Adjust for x>0.  Right now p<0.5, so this is numerically safe.
reflect = gen & (x > 0);
p(reflect) = 1 - p(reflect);

% Make the result exact for the median.
p(x == 0 & ~nans) = 0.5;
