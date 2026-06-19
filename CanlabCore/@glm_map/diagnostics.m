function obj = diagnostics(obj, varargin)
% Compute and report design diagnostics for a glm_map object.
%
% Evaluates the conditioning of the design matrix X (and contrasts C):
% variance inflation factors (VIF) per regressor, contrast VIFs (cVIF),
% per-observation leverage, condition number, rank deficiency, and a
% redundant/near-collinear row report. Results are stored back into the
% object and (optionally) printed.
%
% :Usage:
% ::
%
%     obj = diagnostics(obj, varargin)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object with a design matrix available (obj.X non-empty).
%
% :Optional Inputs:
%
%   **'doverbose' / 'noverbose':**
%        Print a summary table (default true).
%
%   **'vif_thresh', [value]:**
%        VIF warning threshold (default 4).
%
% :Outputs:
%
%   **obj:**
%        glm_map with obj.diagnostics populated (fields:
%        Variance_inflation_factors, Leverages -- same names as
%        fmri_data.regress out.diagnostics -- plus
%        Contrast_variance_inflation_factors, Cooks_distance (per-observation,
%        only when residuals are available), condition_number, rank_deficient,
%        collinearity_report, vif_threshold) and warnings appended.
%
% :Examples:
% ::
%
%     g = glm_map('X', [ones(30,1) zscore((1:30)')], 'level', 2);
%     g = diagnostics(g);
%
% :See also:
%   - VIF, cVIF, fmri_data.regress
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation. Operates on the design only (no fit
%    required), so it can be used to screen a design before fitting.
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
doverbose = true;
if any(strcmpi(varargin, 'noverbose')), doverbose = false; end

vif_thresh = 4;
wh = find(strcmpi(varargin, 'vif_thresh'));
if ~isempty(wh), vif_thresh = varargin{wh(1) + 1}; end

X = obj.X;
if isempty(X)
    error('glm_map:NoDesign', 'No design matrix available (obj.X is empty). Build or supply a design first.');
end

mywarnings = {};

% Diagnostics are collected into a local struct, then stored in
% obj.diagnostics (a nested property) at the end.
d = obj.diagnostics;
if ~isstruct(d), d = struct(); end
d.vif_threshold = vif_thresh;

% -------------------------------------------------------------------------
% Variance inflation factors (per regressor). Field name matches
% fmri_data.regress out.diagnostics.Variance_inflation_factors.
% -------------------------------------------------------------------------
d.Variance_inflation_factors = VIF(X);

if any(d.Variance_inflation_factors > vif_thresh)
    mywarnings{end + 1} = sprintf(['Design multicollinearity: %d regressor(s) have VIF > %g. ' ...
        'Check obj.diagnostics.Variance_inflation_factors and obj.regressor_names.'], ...
        sum(d.Variance_inflation_factors > vif_thresh), vif_thresh);
end

% -------------------------------------------------------------------------
% Contrast variance inflation factors (per contrast), if contrasts defined
% -------------------------------------------------------------------------
d.Contrast_variance_inflation_factors = [];
if ~isempty(obj.contrasts)
    if size(obj.contrasts, 1) ~= size(X, 2)
        mywarnings{end + 1} = sprintf(['Contrast matrix has %d rows but design has %d regressors; ' ...
            'skipping contrast VIFs.'], size(obj.contrasts, 1), size(X, 2));
    else
        % cVIF expects one contrast per row -> transpose [P x K] to [K x P]
        d.Contrast_variance_inflation_factors = cVIF(X, obj.contrasts');
    end
end

% -------------------------------------------------------------------------
% Leverage (per observation)
% -------------------------------------------------------------------------
H = X * pinv(X);
d.Leverages = diag(H)';

if any(abs(zscore(d.Leverages)) >= 3)
    mywarnings{end + 1} = ['Some observations have high leverage (abs(z(leverage)) >= 3); ' ...
        'the fit may be unstable. Check obj.diagnostics.Leverages.'];
end

% -------------------------------------------------------------------------
% Cook's distance (per observation). Measures the influence of each
% observation (image) on the fitted values. Requires residuals, so it is
% computed only for a fitted object whose residuals and sigma are available
% (fit with the 'residuals' option, or via glm_map's automatic Cook's-D
% pass). In the voxelwise GLM it is summarized per observation as the mean
% across voxels of the per-voxel Cook's distance:
%   D_i = mean_v [ r_iv^2 / (k * s_v^2) ] * h_i / (1 - h_i)^2
% with r = residuals, s_v^2 = residual variance (sigma^2), h_i = leverage,
% k = number of parameters.
% -------------------------------------------------------------------------
d.Cooks_distance = [];
have_resid = ~isempty(obj.residuals) && isa(obj.residuals, 'fmri_data') && ~isempty(obj.residuals.dat);
have_sigma = ~isempty(obj.sigma)     && isa(obj.sigma, 'fmri_data')     && ~isempty(obj.sigma.dat);

if have_resid && have_sigma
    h = d.Leverages(:);                     % [n x 1] leverage per observation
    k = size(X, 2);                         % number of parameters
    r = obj.residuals.dat;                  % [voxels x n] residuals
    s2 = obj.sigma.dat(:) .^ 2;             % [voxels x 1] residual variance per voxel

    if size(r, 2) == numel(h)
        good = isfinite(s2) & s2 > 0;       % use voxels with valid residual variance
        if any(good)
            msr  = mean( (r(good, :) .^ 2) ./ s2(good), 1 );   % [1 x n] mean across voxels
            hfac = (h ./ (1 - h) .^ 2)';                       % [1 x n] leverage factor
            d.Cooks_distance = (msr / k) .* hfac;

            if any(d.Cooks_distance > 1)
                mywarnings{end + 1} = sprintf(['%d observation(s) have Cook''s distance > 1 ' ...
                    '(highly influential). Check obj.diagnostics.Cooks_distance.'], sum(d.Cooks_distance > 1));
            end
        end
    end
end

% -------------------------------------------------------------------------
% Conditioning / rank
% -------------------------------------------------------------------------
d.condition_number = cond(X);
d.rank_deficient   = rank(X) < size(X, 2);

if d.rank_deficient
    mywarnings{end + 1} = 'Design matrix X is rank deficient (rank(X) < number of regressors).';
end

% -------------------------------------------------------------------------
% Redundant / near-collinear column report
% -------------------------------------------------------------------------
report = struct();
report.vif_threshold      = vif_thresh;
report.high_vif_columns   = find(d.Variance_inflation_factors > vif_thresh);

% Duplicate (identical) columns
ncol = size(X, 2);
dup_pairs = [];
for a = 1:ncol - 1
    for b = a + 1:ncol
        if isequal(X(:, a), X(:, b))
            dup_pairs(end + 1, :) = [a b]; %#ok<AGROW>
        end
    end
end
report.duplicate_column_pairs = dup_pairs;
if ~isempty(dup_pairs)
    mywarnings{end + 1} = sprintf('%d pair(s) of identical design columns detected (see obj.collinearity_report).', size(dup_pairs, 1));
end

% Near-collinear pairs by pairwise correlation magnitude
R = corrcoef(X);
R(logical(eye(ncol))) = 0;
[ia, ib] = find(triu(abs(R) > 0.95, 1));
report.high_correlation_pairs = [ia ib];

d.collinearity_report = report;

% -------------------------------------------------------------------------
% Store diagnostics struct and warnings
% -------------------------------------------------------------------------
obj.diagnostics = d;
obj.warnings = [obj.warnings(:); mywarnings(:)]';
obj.history{end + 1} = 'diagnostics: computed VIF, cVIF, leverage, condition number, collinearity report';

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------
if doverbose
    rn = obj.regressor_names;
    fprintf('\n  glm_map diagnostics\n  %s\n', repmat('-', 1, 50));
    fprintf('  %-28s %s\n', 'Regressor', 'VIF');
    for i = 1:numel(d.Variance_inflation_factors)
        if i <= numel(rn) && ~isempty(rn{i}), name = rn{i}; else, name = sprintf('R%d', i); end
        flag = ''; if d.Variance_inflation_factors(i) > vif_thresh, flag = '  <-- high'; end
        fprintf('  %-28s %6.2f%s\n', name, d.Variance_inflation_factors(i), flag);
    end
    if ~isempty(d.Contrast_variance_inflation_factors)
        fprintf('  %s\n  Contrast VIFs:\n', repmat('-', 1, 50));
        for i = 1:numel(d.Contrast_variance_inflation_factors)
            if i <= numel(obj.contrast_names) && ~isempty(obj.contrast_names{i})
                name = obj.contrast_names{i};
            else
                name = sprintf('Con%d', i);
            end
            fprintf('  %-28s %6.2f\n', name, d.Contrast_variance_inflation_factors(i));
        end
    end
    fprintf('  %s\n', repmat('-', 1, 50));
    fprintf('  condition number : %.2f\n', d.condition_number);
    fprintf('  rank deficient   : %d\n', d.rank_deficient);
    if ~isempty(d.Cooks_distance)
        [maxcd, wmax] = max(d.Cooks_distance);
        fprintf('  max Cook''s dist  : %.3f (observation %d)\n', maxcd, wmax);
    end
    if ~isempty(mywarnings)
        fprintf('  %d warning(s):\n', numel(mywarnings));
        for i = 1:numel(mywarnings), fprintf('    - %s\n', mywarnings{i}); end
    end
    fprintf('\n');
end

end % diagnostics
