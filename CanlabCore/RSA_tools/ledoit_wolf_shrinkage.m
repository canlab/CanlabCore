function [sigma, shrinkage] = ledoit_wolf_shrinkage(X)
% ledoit_wolf_shrinkage  Regularized covariance via Ledoit-Wolf shrinkage.
%
% Fallback for rsa.stat.covdiag when the Kriegeskorte rsatoolbox is not on
% the path. Shrinks the sample covariance toward a diagonal target with the
% optimal shrinkage intensity from Ledoit & Wolf (2004), J. Multivariate
% Analysis 88: 365-411.
%
% Inputs
%   X  [n x p] data (n observations x p variables)
%
% Outputs
%   sigma      [p x p] regularized covariance estimate
%   shrinkage  scalar in [0, 1] — applied shrinkage intensity
%
% Algorithm follows the rsa.stat.covdiag source: shrinks toward
% target = mean(diag(sample)) * I (constant diagonal). For our use the
% common "diagonal of sample, zero off-diagonal" target is also acceptable;
% we use the constant-diagonal form to match covdiag's conventions.

[n, p] = size(X);
if n < 2
    error('ledoit_wolf_shrinkage:tooFewObs', 'Need at least 2 observations.');
end

% Center
Xc = X - mean(X, 1);

% Sample covariance with ML normalization (matches Ledoit-Wolf derivation)
sample = (Xc' * Xc) / n;

% Target: mean variance on diagonal, 0 off-diagonal
mean_var = mean(diag(sample));
target   = mean_var * eye(p);

% Estimate phi = sum_ij Var(sample_ij)
Y2 = Xc.^2;
phi_mat = (Y2' * Y2) / n - sample.^2;
phi = sum(phi_mat(:));

% Estimate gamma = ||sample - target||_F^2
gamma = norm(sample - target, 'fro')^2;

% Optimal shrinkage intensity (clip to [0, 1])
if gamma > 0
    kappa     = phi / gamma;
    shrinkage = max(0, min(1, kappa / n));
else
    shrinkage = 1;
end

sigma = shrinkage * target + (1 - shrinkage) * sample;

end
