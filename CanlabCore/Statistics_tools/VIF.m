function [V, contrast_V] = VIF(X, varargin)
% Compute variance inflation factors (VIFs) for regressors in a design matrix
%
% :Usage:
% ::
%
%     [V, contrast_V] = VIF(X, [contrast matrix C])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2026 Tor Wager, 3/11/2026
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
% ..
%
% :Inputs:
%
%   **X:**
%        [n x p] design matrix containing regressors. May or may not include
%        an intercept column, though an intercept is typical for the vast
%        majority of GLM design matrices
%
% :Outputs:
%
%   **V:**
%        Row vector of variance inflation factors for regressors in X.
%        VIFs are returned in the same order as the input columns.
%
%
% -------------------------------------------------------------------------
% Background
%
% The variance inflation factor (VIF) quantifies how much the variance of a
% regression coefficient is inflated because of collinearity among predictors.
%
% For predictor j:
%
%       VIF_j = 1 / (1 - R_j^2)
%
% where R_j^2 is the R² obtained by regressing predictor j on all other
% predictors.
%
% Interpretation: 
% VIFs are a multiplier for error variance of regression parameter
% estimates (beta-hats). If VIF = 3 for a regressor, var(beta_hat) = 3x
% larger due to multicollinearity.
% 
% The values below are rough rules of thumb.
%
%       VIF = 1        orthogonal predictors; VIFs are always >= 1
%       VIF ~ 1–3      minimal to moderate collinearity
%       VIF ~ 3–5      moderate to strong collinearity
%       VIF > 10       extreme multicollinearity (rule-of-thumb)
%       VIF = Inf      exact linear dependence
%
% Note: Higher VIFs may be more acceptable in large-sample observational
% designs, depending on study goals.  In experimental designs, VIFs should
% be close to 1. When adding observed correlates (e.g., head motion
% regressors in fMRI), VIFs will increase. We raise our eyebrows when they
% get higher than 3 or so.  
% 
% Note: classical VIFs are formulated for individual predictors, 
% but Jeanette Mumford has formulated an extension for contrasts across predictors, 
% which are typically important in fMRI studies (and other applications).
% See cVIF.m
%
% -------------------------------------------------------------------------
% The Intercept VIF
%
% The intercept also has a VIF that measures how strongly the constant
% vector can be explained by the regressors.  This is not a typical
% formulation in traditional statistics, but is important for many
% applications, particularly in fMRI studies. Researchers often want to
% interpret task activation estimates (beta-hats) compared to 0, i.e., 
% compared to the implicit resting baseline. The stability of those beta-hats 
% depends on the stability of the intercept. The Intercept VIF assesses
% this.
%
% Conceptually:
%
%       VIF_intercept = 1 / (1 - R_0^2)
%
% where R_0^2 is obtained by regressing the intercept column (all ones)
% on the predictors.
%
% If the intercept is unstable (large VIF), the regression coefficients
% become numerically unstable and even their **signs may change** under
% small perturbations of the model or data.
%
% Mean-centering predictors stabilizes the intercept by forcing
%
%       sum(X_centered) = 0
%
% which makes regressors orthogonal to the intercept. This pushes the
% collinearity entirely into the predictor space.
%
% The individual predictor coefficients may still be unstable, but
% contrasts among them are often stable and interpretable.
%
%
% -------------------------------------------------------------------------
% Computational methods for VIFs
%
% There are three essentially equivalent approaches.
%
% 1. Auxiliary regression definition
%
%       VIF_j = 1 / (1 - R_j^2)
%
% Example:
%
%       y = X(:,j);
%       Xother = X(:,setdiff(1:p,j));
%       b = glmfit(Xother, y);
%       yhat = [ones(n,1) Xother]*b;
%       R2 = 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
%       VIF(j) = 1/(1-R2);
%
%
% 2. Inverse correlation matrix
%
%       VIF = diag(inv(corr(intercept(X,'remove'))))'
%
%       Here intercept.m is a CANlab utility for identifying, adding, and
%       removing intercepts.  The intercept must be removed here for
%       inv(corr(X)) to be defined.
%
% 3. Cross-product solve
%
%       XtX = X'*X
%       VIF = ( diag(XtX \ eye(size(XtX))) .* diag(XtX) )'
%
% This function uses the **cross-product solve method** after centering
% predictors, which matches the regression definition exactly. Note that
% non-intercept predictors must be centered for this to match the
% regression and inverse correlation methods.
%
% -------------------------------------------------------------------------
% IMPLEMENTATION NOTES
%
% • Predictors are mean-centered to ensure equivalence with the regression
%   definition of VIF.
%
% • An intercept is added if not present to ensure correct variance
%   partitioning.
%
% • VIFs are returned in the original regressor order.
%
% • If the input design did not include an intercept, it is temporarily added,
% % and its VIF is removed from the output.
%
% -------------------------------------------------------------------------
% Contrast variance inflation factors (cVIF)
%
% In many analyses the scientific question concerns contrasts among
% regression coefficients rather than individual parameters. Even when the
% design matrix contains substantial collinearity and individual
% coefficients are unstable, certain contrasts may still be estimable and
% have well-behaved variance.
%
% The contrast variance inflation factor (cVIF) quantifies how much
% collinearity inflates the variance of a linear contrast of parameters.
%
% For a contrast vector c:
%
%       cVIF_{cβ̂} = c (X'X)^(-1) c'  /  c (X'X ∘ I)^(-1) c'
%
% where:
%
%   X     = centered and standardized design matrix
%   I     = identity matrix
%   ∘     = element-wise (Hadamard) multiplication
%
% The denominator corresponds to a hypothetical "best-case" design in which
% all regressors are orthogonal (i.e., all correlations among regressors are
% set to zero).
%
% Interpretation:
%
%       cVIF = 1     : no variance inflation
%       cVIF > 1     : variance inflated due to collinearity
%       large cVIF   : unstable contrast estimate
%
% Importantly, contrasts may remain estimable even when the design matrix
% itself is rank deficient. For example, in a dummy-coded one-hot model
% with four condition regressors and an intercept, the individual
% coefficients are not uniquely identifiable, but pairwise differences
% between conditions are still estimable contrasts.
%
% If contrasts are supplied to VIF(), they are passed to cVIF() and the
% resulting values are returned in contrast_V.
%
% Reference:
%   Mumford et al. (2025) Imaging Neuroscience.
%
% -------------------------------------------------------------------------
%
% See also: getvif (deprecated except use in GA opt), cVIF, scn_spm_design_check
%
% -------------------------------------------------------------------------
% :Examples:
% ::
%
%    %% Example 1: Dummy-coded 4-condition design
%
%    n = 40;
%    condvec = repelem(1:4,10)';
%
%    D = zeros(n,4);
%    for k=1:4
%        D(:,k) = condvec==k;
%    end
%
%    X = [ones(n,1) D];
%    VIF(X)
%
%   %% Example 1b: Pairwise contrasts in dummy-coded design
% 
%   % Create contrast matrix testing all pairwise differences between 4 conditions
%   C = [
%     0  1 -1  0  0
%     0  1  0 -1  0
%     0  1  0  0 -1
%     0  0  1 -1  0
%     0  0  1  0 -1
%     0  0  0  1 -1
%     ];
% 
%   [V, contrast_V] = VIF(X, C)
% 
% % These contrasts remain estimable even though the design is collinear
% % because they compare conditions relative to one another.
%
%    %% Example 2: Effects-coded version
%
%    C = zeros(n,3);
%    for k=1:3
%        C(:,k) = (condvec==k) - (condvec==4);
%    end
%
%    Xeff = [ones(n,1) C];
%    VIF(Xeff)
%
%    %% Example 3: Random regressors
%
%    X = rand(100,4); % mean is positive here, overlapping with intercept
%    X(:, end+1) = 1; % intercept
%    VIF(X)
%    disp('The intercept VIF is high here because the regressors together can approximate it; its estimate is unstable')
%
%    %% Centered version
%
%    Xc = X(:, 1:end-1) - mean(X(:, 1:end-1));
%    Xc(:, end+1) = 1; % intercept
%    VIF(Xc)
%    disp('The intercept VIF is 1 here because the collinearity has been forced into the regressors by centering')
%
%    %% Example 4: Random regressors, but one regressor perfectly
%    % correlated with a combo of two others
%
%    X = rand(100,4); % mean is positive here, overlapping with intercept
%    X(:, 3) = .4 * X(:, 1) + .6 * X(:, 2);
%    X(:, end+1) = 1; % intercept
%    VIFs = VIF(X)
%    VIFs(4:5) % The independent random regressor and intercept are OK, no perfect collinearity
% -------------------------------------------------------------------------

[n,p] = size(X);

% Optional contrasts
contrast_V = [];
C = [];
if ~isempty(varargin), C = varargin{1}; end

% Detect intercept column (if any)
ic = intercept(X,'which');

% Track whether intercept was originally present
had_intercept = ~isempty(ic);

% Remove intercept to isolate predictor block
Xp = intercept(X,'remove');

% ---------------------------------------------------------------------
% Compute intercept VIF using algebraic projection BEFORE centering
% ---------------------------------------------------------------------

if had_intercept

    % Response vector for intercept regression (column of ones)
    y = X(:, ic);

    % Cross-product inverse for predictors
    XtX_inv = inv(Xp' * Xp);

    % Algebraic R^2 for regression of intercept on predictors
    R2_intercept = (y' * Xp * XtX_inv * Xp' * y) / (y' * y);

    % Convert to VIF
    if abs(1 - R2_intercept) < 1e-12
        intercept_VIF = Inf;
    else
        intercept_VIF = 1 / (1 - R2_intercept);
    end

end

% ---------------------------------------------------------------------
% Compute VIFs for regressors using centered design matrix X with intercept enforced
% ---------------------------------------------------------------------

% If no intercept exists, add one temporarily
if ~had_intercept
    X = [X ones(n,1)];
    ic = size(X,2);
end

% Mean-center predictors to ensure correct VIF calculation
Xp = Xp - mean(Xp);

% Reassemble design matrix with intercept as final column
Xc = [Xp X(:,ic)];

% Compute cross-product matrix
XtX = Xc' * Xc;

% Solve system instead of computing matrix inverse
% This is numerically more stable
Vall = ( diag(XtX \ eye(size(XtX))) .* diag(XtX) )';

% If original model did not contain an intercept,
% remove the intercept VIF before returning
if ~had_intercept
    V = Vall(1:end-1);
else
    V = Vall;
    V(end) = intercept_VIF;  % this will be moved later if needed
end

% If original model contained an intercept and it was not the final column,
% rearrange the VIFs so the intercept VIF is in the same location in the
% design matrix

if had_intercept && ic ~= p

    % Vall currently assumes intercept is last column
    % Reconstruct original order

    % Move intercept VIF to its original column location
    V_reordered = zeros(size(V));

    idx = setdiff(1:p, ic);
    V_reordered(idx) = V(1:end-1);
    V_reordered(ic) = V(end);

    V = V_reordered;

end

% ---------------------------------------------------------------------
% Compute contrast VIFs if contrasts were provided
% ---------------------------------------------------------------------

if nargin > 1 && ~isempty(C)

    % Ensure contrasts are row-wise
    if size(C, 2) ~= p
        error('Contrast matrix must have same number of columns as X; Each row is a contrast vector');
    end

    % Compute contrast-based VIFs using cVIF
    contrast_V = cVIF(X, C);

end

end % VIF 