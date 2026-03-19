function [lambda_single, lambda_mean, stats] = reliability_estimate_lambda(W, varargin)
% reliability_estimate_lambda
%
% Estimate reliability (lambda) under a classical errors-in-variables (EIV) model
% using repeated measurements.
%
% -------------------------------------------------------------------------
% THEORY
% -------------------------------------------------------------------------
% Assume repeated measurements:
%
%   w_ij = x_i + delta_ij
%
% where:
%   x_i        = true latent value for subject i
%   delta_ij   = measurement error (mean 0, variance sigma_delta^2)
%
% Then reliability (lambda) is:
%
%   lambda = Var(x) / (Var(x) + Var(delta))
%
% This is estimated using an intraclass correlation coefficient (ICC):
%
%   ICC_single = (MS_between - MS_within) / (MS_between + (k - 1)*MS_within)
%
% Reliability of the mean across k measurements:
%
%   lambda_mean = Var(x) / (Var(x) + Var(delta)/k)
%               = (k * lambda) / (1 + (k - 1)*lambda)
%
% Relationship with ICC (intraclass correlation coefficient):
% lambda_single is equal to ICC(1, 1) 
% Use this if you are using a single measure in association with other 
% variables going forward
%
% lambda_mean is equal to ICC(1, k) 
% Use this if you are using the average over measures in association 
% with other variables going forward
%
% -------------------------------------------------------------------------
% RELATION TO ATTENUATION BIAS
% -------------------------------------------------------------------------
% If y is regressed on noisy predictor w:
%
%   beta_hat = beta_true * lambda
%
% So:
%
%   beta_true ≈ beta_hat / lambda
%
% For correlations:
%
%   r_obs = rho_true * sqrt(lambda_x * lambda_y)
%
% So:
%
%   rho_true ≈ r_obs / sqrt(lambda_x * lambda_y)
%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% W : n x k matrix
%     rows = observations (subjects)
%     cols = repeated measurements
%
% -------------------------------------------------------------------------
% OPTIONAL INPUTS (varargin)
% -------------------------------------------------------------------------
% 'type' : {'oneway', 'consistency'}
%
%   'oneway' (default):
%       ICC(1,1)-type estimate (exchangeable measurements)
%
%   'consistency':
%       ICC(3,1)-type estimate (removes column mean differences)
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% lambda_single : reliability of a single measurement
% lambda_mean   : reliability of the mean across k measurements
% stats         : struct with MS_between, MS_within, k
%
% -------------------------------------------------------------------------
% EXAMPLE 1: TWO REPEATS (k = 2), ATTENUATION CORRECTION
% -------------------------------------------------------------------------
% n = 200;
% x = randn(n,1);
% y = 0.6 * x + randn(n,1)*0.5;
%
% % noisy measurements of x
% w1 = x + randn(n,1)*0.8;
% w2 = x + randn(n,1)*0.8;
%
% W = [w1 w2];
%
% [lambda_single, lambda_mean] = reliability_estimate_lambda(W);
%
% % observed correlation
% r_obs = corr(mean(W,2), y);
%
% % "Ground truth" in this sim, to compare against:
% r_true = corr(x, y); 
%
% % corrected correlation (using mean reliability)
% rho_est = r_obs / sqrt(lambda_mean);
%
% lam_icc_single = ICC(1, 'single', W);
% lam_icc_mean = ICC(1, 'k', W);
% rho_est_icc = r_obs / sqrt(lam_icc_mean);
%
% fprintf('Observed r = %.3f\n', r_obs);
% fprintf('Corrected rho = %.3f, using icc = %.3f, true = %.3f\n', rho_est, rho_est_icc, r_true);
%
% -------------------------------------------------------------------------
% EXAMPLE 2: k = 5 REPEATS
% -------------------------------------------------------------------------
% n = 200;
% k = 5;
%
% x = randn(n,1);
% y = 0.5 * x + randn(n,1)*0.6;
%
% W = repmat(x,1,k) + randn(n,k)*0.7;
%
% [lambda_single, lambda_mean] = reliability_estimate_lambda(W);
%
% % predictor = mean across repeats
% wbar = mean(W,2);
%
% r_obs = corr(wbar, y);
%
% % "Ground truth" in this sim, to compare against:
% r_true = corr(x, y); 
%
% % correct for attenuation
% rho_est = r_obs / sqrt(lambda_mean);
%
% lam_icc_single = ICC(1, 'single', W);
% lam_icc_mean = ICC(1, 'k', W);
% rho_est_icc = r_obs / sqrt(lam_icc_mean);
%
% fprintf('Observed r = %.3f\n', r_obs);
% fprintf('Corrected rho = %.3f, using icc = %.3f, true = %.3f\n', rho_est, rho_est_icc, r_true);
%
% -------------------------------------------------------------------------
%
% See also: icc.m, plot_correlation_matrix.m, line_plot_multisubject.m

% Defaults
type = 'oneway';

% Parse varargin
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'type'
            type = varargin{i+1};
    end
end

[n, k] = size(W);

% Remove column means for consistency ICC (ICC(3))
if strcmpi(type, 'consistency')
    W = W - mean(W,1);
end

% Means
subject_means = mean(W,2);
grand_mean = mean(W(:));

% Sum of squares
SS_between = k * sum((subject_means - grand_mean).^2, 'omitnan');
SS_within  = sum(sum((W - subject_means).^2, 'omitnan'), 'omitnan');

% Degrees of freedom
df_between = n - 1;
df_within  = n * (k - 1);

% Mean squares
MS_between = SS_between / df_between;
MS_within  = SS_within / df_within;

% ICC (single measure)
lambda_single = (MS_between - MS_within) / ...
                (MS_between + (k - 1)*MS_within);

% ICC (average measure)
lambda_mean = (MS_between - MS_within) / MS_between;

% Clamp to valid range
lambda_single = max(min(lambda_single, 1), 0);
lambda_mean   = max(min(lambda_mean, 1), 0);

% Output stats
stats = struct('MS_between', MS_between, ...
               'MS_within', MS_within, ...
               'k', k);

end