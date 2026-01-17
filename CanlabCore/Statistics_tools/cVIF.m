function cVIF = compute_cVIF(X, c)
% compute_cVIF  Calculate contrast-based variance inflation factor (cVIF)
%
%   cVIF = compute_cVIF(X, c)
%
%   This function computes the contrast-based variance inflation factor
%   (cVIF) for one or more contrasts of parameter estimates in a general
%   linear model, as defined in Mumford et al. (2025, Imaging Neuroscience).
%   cVIF quantifies how much the variance of a contrast is inflated by
%   collinearity among regressors in X.
%
%   The implementation follows Appendix C of Mumford et al.:
%
%       cVIF_{cβ̂} = c (X'X)^(-1) c'  /  c (X'X ∘ I)^(-1) c'
%
%   where X is the (centered and scaled) design matrix, c is a row vector
%   (or matrix) of contrast weights, I is the identity matrix, and ∘
%   denotes element-wise multiplication. The denominator corresponds to a
%   “best-case” design in which all off-diagonal elements (correlations)
%   between regressors are set to zero.
%
%   Usage:
%       cVIF = compute_cVIF(X, c);
%
%   Inputs:
%   -------
%   X : [N × P] matrix
%       Design matrix with N observations and P regressors. This function
%       will (a) remove any constant column (intercept) and (b) mean-center
%       and standardize each remaining column before computing cVIF, as
%       recommended in Mumford et al. (Appendix C).
%
%   c : [K × P] or [1 × P] vector/matrix
%       Contrast(s) over the P regressors. Each row is a contrast. If c is
%       provided as a column vector, it will be transposed internally.
%
%   Outputs:
%   --------
%   cVIF : [K × 1] vector
%       Contrast-based variance inflation factor for each contrast (one
%       value per row of c). Values > 1 indicate variance inflation due to
%       collinearity among regressors in X.
%
%   Example:
%   --------
%       % Suppose X is an N×P design matrix and we want cVIF for a
%       % condition-difference contrast [0 1 -1 0 ...]
%       c  = [0 1 -1 0];
%       v  = compute_cVIF(X, c);
%
%       % Create a sample fMRI design matrix and test VIF and cVIF
%
%        ons = create_random_onsets(100, 3, [.2 .2 .2 .2], 2);
%        [X,d,out,handles] = plotDesign(ons,[], 1.3);
%        getvif(X)
%        cVIF(X, [1 1 -1 -1 0])
%        cVIF(X, [.5 .5 -.5 -.5 0])
%        C = create_orthogonal_contrast_set(4)
%        C = [C [0 0 0]'];
%        cVIF(X, C)
%       
%       % Create a perfectly colinear design matrix and test
%        X = [X(:, 1) -X(:, 1)]
%        getvif(X)
%        cVIF(X, [1 -1])
%        cVIF(X, [-.5 .5])
%        cVIF(X, [1 1])
% 
%   Notes:
%   ------
%   • X should contain one column per condition/regressor; parametrically
%     modulated regressors with only two levels should be reparameterized
%     as separate columns per level.
%   • The function uses pinv (Moore-Penrose pseudoinverse) for numerical
%     stability, so it can be applied to mildly rank-deficient designs.
%
%   References:
%   -----------
%   Mumford, J. A., Demidenko, M. I., Bjork, J. M., et al. (2025).
%   Unintended bias in the pursuit of collinearity solutions in fMRI
%   analysis. Imaging Neuroscience, 3, 1–24.
%
%   See also: pinv, cov

% -------------------------------------------------------------------------
% Validate inputs
% -------------------------------------------------------------------------

if nargin < 2
    error('compute_cVIF requires two inputs: X (design) and c (contrast).');
end

if ~ismatrix(X)
    error('X must be a 2-D matrix (N × P).');
end

[N, P] = size(X);

% Ensure c is 2-D, with one contrast per row
c = squeeze(c);
if isvector(c)
    c = c(:)';  % make row
end

if size(c, 2) ~= P
    error('Number of columns in c (%d) must match number of regressors in X (%d).', ...
          size(c,2), P);
end

% -------------------------------------------------------------------------
% Remove intercept (constant) column if present
% -------------------------------------------------------------------------

col_std = std(X, 0, 1);
is_const = col_std < (sqrt(eps) * max(col_std));

if any(is_const)
    X = X(:, ~is_const);
    c = c(:, ~is_const);
    P = size(X, 2);
end

% -------------------------------------------------------------------------
% Center and scale regressors (as in Appendix C)
% -------------------------------------------------------------------------

% Mean-center
X = bsxfun(@minus, X, mean(X, 1));

% Standardize by column standard deviation
col_std = std(X, 0, 1);
% Avoid division by zero
col_std(col_std == 0) = 1;
X = bsxfun(@rdivide, X, col_std);

% -------------------------------------------------------------------------
% Compute X'X and its “decorrelated” version
% -------------------------------------------------------------------------

XX = X' * X;                  % P × P
invXX = pinv(XX);             % worst-case (true collinearity)

% Best-case: zero out off-diagonal elements of X'X
XX_diag = diag(diag(XX));     % P × P diagonal matrix
invXX_diag = pinv(XX_diag);   % best-case (no collinearity)

% -------------------------------------------------------------------------
% Compute cVIF for each contrast
% -------------------------------------------------------------------------

% Numerator: c (X'X)^(-1) c'
V_num = c * invXX * c';       % K × K (symmetric)
% Denominator: c (X'X ∘ I)^(-1) c'
V_den = c * invXX_diag * c';  % K × K (diagonal if c rows independent)

% Extract per-contrast variances (diagonal elements)
v_num = diag(V_num);
v_den = diag(V_den);

% Protect against numerical issues
if any(v_den <= 0)
    warning('compute_cVIF:NonPositiveDenominator', ...
        ['One or more best-case contrast variances are non-positive. ', ...
         'cVIF may be unreliable for those contrasts.']);
end

cVIF = v_num ./ v_den;

end % function
