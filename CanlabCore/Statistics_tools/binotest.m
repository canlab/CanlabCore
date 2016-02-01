function RES = binotest(X, p)
% Test the number of "hits" in each column of X against a null-hypothesis
% proportion p, using a binomial test.
%
% :Usage:
% ::
%
%     RES = binotest(X, p)
% 
% Assumes elements of each column of X are independent Bernoulli trials.
%
% :Inputs:
%
%   **X:**
%        is a matrix of "hits" and "misses", coded as 1s and 0s.
%
%   **p:**
%        is the null hypothesis proportion of "hits", e.g., often p = 0.5
%
% ..
%    Tor Wager, April 2011
% ..

X = double(X); % just in case

% Binomial test
% "By chance, prob of getting this many or more hits are..."
n = size(X, 1) .* ones(1, size(X, 2));
hits = sum(X);
prop = hits ./ n;

%binop = 2 * (1 - binocdf(hits, n, .5)); % two-tailed P-values; wrong
%because does not count prob of THIS MANY or more; counts prob of more.

binop = 2 * min(binocdf(hits, n, p), (1 - binocdf(hits - 1, n, p)));

if binop > 1, binop = 1.0; end % see https://github.com/canlab/CanlabCore/issues/9

binose = sqrt(prop .* (1-prop) ./ n); % based on estimated proportions

RES = struct('n', n, 'hits', hits, 'prop', prop, 'p_val', binop, 'SE', binose);

end
