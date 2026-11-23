function [Xw, info] = whiten_within_subject(X, method)
% whiten_within_subject  Apply within-subject whitening to a replicate matrix.
%
% Generalizes the recipe from `10252024 RSA Whitening Walkthrough.mlx`:
%   1) Estimate regularized covariance via covdiag (or stock Ledoit-Wolf fallback).
%   2) Compute the inverse square root of the covariance via SVD.
%   3) Apply: X_whitened = X * sigma^{-1/2}.
%
% Inputs
%   X       [n x p] data (rows = replicates / observations, cols = features)
%   method  'covdiag' (use rsa.stat.covdiag if available) | 'diag' (variance only) | 'none'
%
% Outputs
%   Xw    [n x p] whitened data
%   info  struct with .method, .shrinkage, .sigma

if nargin < 2 || isempty(method), method = 'covdiag'; end
info = struct('method', method, 'shrinkage', [], 'sigma', []);

switch lower(method)
    case 'none'
        Xw = X;
        return
    case 'diag'
        v = var(X, 0, 1);
        v(~isfinite(v) | v <= eps) = 1;
        Xw = X ./ sqrt(v);
        info.sigma = diag(v);
        info.shrinkage = 1;
        return
    case 'covdiag'
        % Try rsa.stat.covdiag first, then fall back
        sigma = [];
        shrinkage = [];
        try
            if ~isempty(which('rsa.stat.covdiag'))
                [sigma, shrinkage] = rsa.stat.covdiag(X);
            elseif ~isempty(which('covdiag'))
                [sigma, shrinkage] = covdiag(X);
            end
        catch
            sigma = [];
        end
        if isempty(sigma)
            [sigma, shrinkage] = ledoit_wolf_shrinkage(X);
        end
    otherwise
        error('whiten_within_subject:badMethod', ...
            'method must be ''covdiag'', ''diag'', or ''none''; got %s.', method);
end

% Inverse square root via SVD (numerically stable)
[U, S, V] = svd(sigma);
s = diag(S);
% Guard against zero eigenvalues
s(s <= eps) = eps;
sigma_inv_sqrt = U * diag(1 ./ sqrt(s)) * V';

Xw = (sigma_inv_sqrt * X')';

info.sigma = sigma;
info.shrinkage = shrinkage;

end
