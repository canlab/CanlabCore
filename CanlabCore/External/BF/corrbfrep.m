function bf10 = corrbfrep(r,n)
%
% bf10 = corrbfrep(r,n)
%
% Calculates replication correlation Bayes Factor using a uniform prior
% between 0-1 (one sided test) as used by Boekel et al, 2014, Cortex.
%

% Function to be integrated
F = @(rho,r,n) ((1-rho.^2).^(1./2*(n-1))) ./ (1-rho*r).^(n-3./2);

% Bayes factor calculation
bf10 = integral(@(rho) F(rho,r,n), 0, 1);