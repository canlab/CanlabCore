function bf10 = corrbf(r,n)
%
% bf10 = corrbf(r,n)
%
% Calculates JZS Bayes factor for correlation r and sample size n.
% This quantifies the evidence in favour of the alternative hypothesis.
% See Wetzels & Wagemakers, 2012, Psychon Bull Rev for details.
%

% Function to be integrated - updated after personal communication from EJW (old one didn't work with large N)
F = @(g,r,n) exp(((n-2)./2).*log(1+g)+(-(n-1)./2).*log(1+(1-r.^2).*g)+(-3./2).*log(g)+-n./(2.*g));

% Bayes factor calculation
bf10 = sqrt((n/2)) / gamma(1/2) * integral(@(g) F(g,r,n),0,Inf);