function bf10 = binombf(k,n,p)
%
% bf10 = binombf(r,n,p)
%
% Calculates Bayes factor for binomial test with k successes, n trials and base probability p.
% This is based on the Bayes Factor example on Wikipedia (http://en.wikipedia.org/wiki/Bayes_factor).
% It is the same as the default prior assumed by Rouder's web applet (http://pcl.missouri.edu/bf-binomial).
%
warning off     'MATLAB:nchoosek:LargeCoefficient'
% Function to be integrated - updated after personal communication from EJW (old one didn't work with large N)
F = @(q,k,n,p) nchoosek(n,k) .* q.^k .* (1-q).^(n-k);

% Bayes factor calculation
bf01 = (nchoosek(n,k) .* p.^k .* (1-p).^(n-k)) / integral(@(q) F(q,k,n,p),0,1);

% Invert Bayes Factor
bf10 = 1 / bf01;