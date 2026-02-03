function [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F] = fit_gls2(y, X, c, p)
% Fit a linear model using FGLS with an AR(p) innovations model.
%
% Inputs:
%   y : (T x 1) dependent variable
%   X : (T x K) design matrix
%   c : (K x Q) contrast matrix (each column is a contrast). Can be [].
%   p : AR order
%
% Outputs:
%   beta      : (K x 1) coefficients
%   t, pvals  : (K x 1) t-stats and p-values for beta
%   convals   : (Q x 1) contrast estimates
%   con_t     : (Q x 1) contrast t-stats
%   con_pvals : (Q x 1) contrast p-values
%   sigma     : residual std estimate (SSE/df)^0.5
%   Phi       : (p x 1) AR coefficients (standard sign convention)
%   df        : residual df = T - rank(X)
%   stebeta   : (K x 1) SE of beta
%   conste    : (Q x 1) SE of contrasts
%   F         : optional F-stat for joint test c' * beta = 0 (if c provided)

% Ensure inputs are double
y = double(y);
X = double(X);

T = size(X, 1);
K = size(X, 2);

% --- Step 1: FGLS estimation (AR(p) innovations) ---
[beta, ~, EstCoeffCov] = fgls(X, y, InnovMdl="AR", ARLags=p, Intercept=false);

% If fgls returned a reduced parameter vector (e.g., due to rank issues),
% try to align dimensions defensively.
beta = beta(:);

if size(beta, 1) ~= K
    % Attempt: keep nonzero entries (matches your original logic)
    beta = beta(beta ~= 0);

    nonZeroRows = any(EstCoeffCov ~= 0, 2);
    nonZeroCols = any(EstCoeffCov ~= 0, 1);
    EstCoeffCov = EstCoeffCov(nonZeroRows, nonZeroCols);
end

% --- Step 2: residuals and AR coefficients (for debugging/return only) ---
resid = y - X * beta;

% aryule returns polynomial A = [1 a1 ... ap] for A(z)*x = e.
% Usual AR form is x_t = sum(phi_k x_{t-k}) + e_t, with phi = -a.
A = aryule(resid, p);
Phi = -A(2:end).';
Phi = Phi(:);

% --- Step 3: SEs, t-stats ---
stebeta = sqrt(diag(EstCoeffCov));
t = beta ./ stebeta;

% --- Step 4: residual df (THIS is the df for first-level GLM t-tests) ---
df = T - rank(X);
df = max(df, 1);

% Residual variance estimate (consistent with df)
SSE = sum(resid.^2);
sigma = sqrt(SSE / df);

% --- Step 5: coefficient p-values ---
pvals = 2 * (1 - tcdf(abs(t), df));

% --- Step 6: contrasts (if provided) ---
convals   = [];
conste    = [];
con_t     = [];
con_pvals = [];
F         = [];

if ~isempty(c)
    c = double(c);
    % Contrast estimate and SE
    convals = (c' * beta);
    conste  = sqrt(diag(c' * EstCoeffCov * c));
    con_t   = convals ./ conste;

    % Use the SAME df as above (residual df)
    con_pvals = 2 * (1 - tcdf(abs(con_t), df));

    % Optional joint F-test for H0: c' * beta = 0
    % Works best when c has multiple columns (Q > 1).
    Q = size(c, 2);
    if Q >= 1
        CVC = c' * EstCoeffCov * c;
        % Guard against singularity
        if rcond(CVC) > 1e-12
            F = (convals' / CVC * convals) / Q;
        else
            F = NaN;
        end
    end
end

end
