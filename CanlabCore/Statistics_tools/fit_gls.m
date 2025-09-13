function [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F] = fit_gls(y, X, c, p)
% Fit a linear model using FGLS with an AR(p) model for innovations
%
% :Usage:
% ::
%     [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F] = fit_gls(y, X, c, p)
%
% :Inputs:
%   **y:**     Dependent variable (T x 1 vector)
%   **X:**     Design matrix (T x param matrix)
%   **c:**     Contrast vector(s) (param x # contrasts matrix)
%   **p:**     Order of AR model for innovations.
%
% :Outputs:
%   **beta:**  Coefficient estimates
%   **t:**     t-values for coefficients
%   **pvals:** p-values for coefficients
%   **convals, conste, con_t, con_pvals:** Contrast-related outputs
%   **sigma:** Residual standard deviation
%   **Phi:**   AR coefficients of residuals
%   **df:**    Degrees of freedom using Satterthwaite approximation
%   **stebeta:** Standard errors of coefficients
%   **F:**     F-statistic for the model

% Ensure inputs are double
y = double(y);
X = double(X);

% Step 1: Estimate coefficients using FGLS
[beta,~, EstCoeffCov] = fgls(X, y, InnovMdl="AR", ARLags=p);

if size(beta, 1) ~= size(X, 2)
    beta = beta(beta ~= 0);
    nonZeroRows = any(EstCoeffCov ~= 0, 2); 
    nonZeroCols = any(EstCoeffCov ~= 0, 1); 
    EstCoeffCov = EstCoeffCov(nonZeroRows, nonZeroCols);
end


% Step 2: Calculate residuals and AR(p) coefficients
resid = y - X * beta;

% Convert residuals to AR coefficients (optional for debugging)
Phi = arma2ar(resid, zeros(1, p)); % AR coefficients using arma2ar

% Step 3: Compute standard errors, t-values, and residual variance
stebeta = sqrt(diag(EstCoeffCov));
t = beta ./ stebeta;
sigma = sqrt(var(resid)); % Residual standard deviation

% Step 4: Contrast analysis
convals = [];
conste = [];
con_t = [];
con_pvals = [];
F = [];

if ~isempty(c)
    % Contrast values
    convals = (c' * beta);
    conste = sqrt(diag(c' * EstCoeffCov * c));
    con_t = convals ./ conste;
    con_pvals = 2 * (1 - tcdf(abs(con_t), size(X, 1) - size(X, 2)));
end

% Step 5: Degrees of freedom (Satterthwaite approximation)
df = (trace(EstCoeffCov)^2) / trace(EstCoeffCov * EstCoeffCov);

% Step 6: Model F-statistic (optional)
if nargout > 6
    SSE = sum(resid.^2);
    SST = sum((y - mean(y)).^2);
    SSM = SST - SSE;
    dfSSM = size(c, 2) - 1; % Degrees of freedom for contrasts
    F = (SSM / dfSSM) / (SSE / df);
end

% Step 7: p-values for coefficients
pvals = 2 * (1 - tcdf(abs(t), df));

end

