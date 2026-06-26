function obj = fit(obj, data, varargin)
% Fit a mass-univariate GLM and populate result maps in a glm_map object.
%
% Scikit-learn-style fit: builds the design (if event/1st-level and not yet
% built), runs fmri_data.regress on the data, and unpacks the outputs into
% the object's statistic_image result maps (betas, t, contrast_estimates,
% contrast_t), plus sigma, dfe, and design diagnostics.
%
% :Usage:
% ::
%
%     obj = fit(obj, data, varargin)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object with a design specified (either a .design event
%        model or a direct .X matrix), and optionally contrasts.
%
%   **data:**
%        An fmri_data object whose images (columns of .dat) correspond to the
%        rows of X (a within-run timeseries for level 1, contrast/subject
%        images for level 2).
%
% :Optional Inputs:
%
%   **'robust':**
%        Use robust (bisquare) regression. Passed through to fmri_data.regress.
%
%   **'AR', [order]:**
%        Use an autoregressive AR(p) error model of the given order. Only
%        valid when obj.is_timeseries is true. AR(4) is recommended for
%        typical BOLD time series (it captures the low-order serial
%        autocorrelation better than AR(1) and approximates an ARMA(1,1)
%        process). Requires two additional toolboxes, which are checked for
%        up front and reported clearly if missing:
%          - Econometrics Toolbox (fgls)        — the GLS fit
%          - Signal Processing Toolbox (aryule) — AR-coefficient estimation
%        Estimation loops over voxels, so whole-brain AR fits are slow.
%
%   **'pthresh', [p-value]:**
%        Initial threshold p-value applied to the t/contrast_t maps
%        (default 0.001).
%
%   **'thresh_type', 'unc' | 'fdr':**
%        Initial threshold type (default 'unc').
%
%   **'residuals':**
%        Also store residuals in obj.residuals.
%
%   **'display':**
%        Show fmri_data.regress orthviews (default off; glm_map suppresses it).
%
%   **'doverbose' / 'noverbose':**
%        Toggle verbose output (default true).
%
% :Outputs:
%
%   **obj:**
%        The glm_map object with result maps and diagnostics populated, and
%        is_fitted == true.
%
% :Examples:
% ::
%
%     dat = load_image_set('emotionreg');
%     g   = glm_map('X', [ones(30,1) zscore((1:30)')], 'level', 2, ...
%                   'regressor_names', {'intercept' 'cov'});
%     g   = add_contrasts(g, [0 1], {'cov_effect'});
%     g   = fit(g, dat);
%     table(g, 'contrast'); montage(g, 'contrast_t');
%
% :See also:
%   - fmri_data.regress, build_design, diagnostics, threshold
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation. Thin orchestration layer over
%    fmri_data.regress: glm_map owns the design/diagnostics/result container,
%    regress remains the compute engine.
% ..

% -------------------------------------------------------------------------
% Parse options
% -------------------------------------------------------------------------
doverbose   = ~any(strcmpi(varargin, 'noverbose'));
do_robust   = any(strcmpi(varargin, 'robust'));
do_resid    = any(strcmpi(varargin, 'residuals')) || any(strcmpi(varargin, 'residual'));
do_display  = any(strcmpi(varargin, 'display'));

pthresh = 0.001;
wh = find(strcmpi(varargin, 'pthresh'));
if ~isempty(wh), pthresh = varargin{wh(1) + 1}; end

thresh_type = 'unc';
wh = find(strcmpi(varargin, 'thresh_type'));
if ~isempty(wh), thresh_type = varargin{wh(1) + 1}; end
if ~ismember(lower(thresh_type), {'unc', 'fdr'})
    error('glm_map:BadThreshType', 'thresh_type must be ''unc'' or ''fdr''.');
end

ar_order = 0;
wh = find(strcmpi(varargin, 'AR'));
if ~isempty(wh), ar_order = varargin{wh(1) + 1}; end

% -------------------------------------------------------------------------
% Validate inputs
% -------------------------------------------------------------------------
if ~isa(data, 'fmri_data')
    error('glm_map:BadData', 'data must be an fmri_data object.');
end

if ar_order > 0 && ~obj.is_timeseries
    error('glm_map:ARnotTimeseries', ...
        'AR error models require obj.is_timeseries == true (within-run timeseries data).');
end

% Check AR-model dependencies up front and fail with a clear, actionable
% message rather than partway through the per-voxel fit. AR(p) estimation
% (fmri_data.regress -> fit_gls2) needs fgls (Econometrics Toolbox) for the
% GLS fit and aryule (Signal Processing Toolbox) for the AR coefficients.
if ar_order > 0
    missing = {};
    if isempty(which('fgls')),   missing{end + 1} = 'fgls (Econometrics Toolbox)'; end
    if isempty(which('aryule')), missing{end + 1} = 'aryule (Signal Processing Toolbox)'; end
    if ~isempty(missing)
        error('glm_map:ARDependencyMissing', ...
            ['AR(%d) error models require functions that are not on your path:\n  - %s\n' ...
             'Install/license the listed toolbox(es), or fit without ''AR'' ' ...
             '(OLS, or add ''robust'' to down-weight outlier time points).'], ...
            ar_order, strjoin(missing, '\n  - '));
    end
end

% Ensure a design matrix is available; build it for event/1st-level models
if isempty(obj.X)
    if obj.level == 1 && ~isempty(obj.design)
        obj = build_design(obj);
    else
        error('glm_map:NoDesign', ...
            'No design available. Supply obj.X (direct mode) or an obj.design to build (event mode).');
    end
end

X = obj.X;

if size(data.dat, 2) ~= size(X, 1)
    error('glm_map:SizeMismatch', ...
        'data has %d images but the design has %d rows. They must match.', size(data.dat, 2), size(X, 1));
end

% -------------------------------------------------------------------------
% Intercept and contrast preparation
% -------------------------------------------------------------------------
% Ensure the design has an intercept and that the regressor names and the
% contrast matrix line up with the design columns BEFORE fitting:
%   - A single overall intercept is enforced as the LAST column (intercept.m
%     'end' semantics: an existing intercept is moved, never duplicated; one is
%     added if none is present). Designs that already carry per-run baseline
%     columns (multi-session first-level) are left as they are.
%   - A contrast vector that is one (or more) elements short gets a 0 weight
%     for the intercept column(s), so contrasts can be entered over the
%     regressors of interest without an explicit intercept entry.
% The model is then fit first (all betas estimated) and contrasts are applied
% afterwards by regress to form the contrast maps.
rn = obj.regressor_names;
C  = obj.contrasts;                          % [regressors x contrasts]
[X, rn, C] = local_prepare_design(obj, X, rn, C);

% Persist the prepared design on direct-mode objects so g.X / g.regressor_names
% match the fit. For event designs, X comes from .design and is unchanged here.
if isempty(obj.design)
    obj.Xdirect = X;
    obj.regressor_names_direct = rn;
end
obj.contrasts = C;

% -------------------------------------------------------------------------
% Assemble regress() argument list
% -------------------------------------------------------------------------
data.X = X;

% Threshold p-value must precede the 'unc'/'fdr' keyword (regress convention).
% The intercept is handled above, so tell regress not to add another.
regress_args = {pthresh, thresh_type, 'nointercept'};

if do_robust,  regress_args{end + 1} = 'robust'; end
if ar_order > 0, regress_args = [regress_args, {'AR', ar_order}]; end
% Always request residuals so diagnostics can compute Cook's distance; they
% are dropped again below unless the caller asked to keep them (do_resid).
regress_args{end + 1} = 'residual';

% Suppress regress's own orthviews unless explicitly requested
if ~do_display, regress_args{end + 1} = 'nodisplay'; end
if ~doverbose,  regress_args{end + 1} = 'noverbose'; end

% Regressor names (length must match number of design columns)
if ~isempty(rn) && numel(rn) == size(X, 2)
    regress_args = [regress_args, {'names', rn}];
end

% Contrasts: regress expects C as [regressors x contrasts] (our storage)
if ~isempty(C)
    regress_args = [regress_args, {'C', C}];
    if ~isempty(obj.contrast_names)
        regress_args = [regress_args, {'contrast_names', obj.contrast_names}];
    end
end

if ~isempty(obj.analysis_name)
    regress_args = [regress_args, {'analysis_name', obj.analysis_name}];
end

% -------------------------------------------------------------------------
% Run the regression. fmri_data.regress now returns a glm_map object whose
% properties mirror the historical out-struct fields.
% -------------------------------------------------------------------------
out = regress(data, regress_args{:});

% -------------------------------------------------------------------------
% Unpack result maps and the nested input structs
% -------------------------------------------------------------------------
obj.betas = out.betas;
obj.t     = out.t;
obj.sigma = out.sigma;
obj.df    = out.df;
obj.input_parameters     = out.input_parameters;
obj.input_image_metadata = out.input_image_metadata;
obj.contrast_summary_table = out.contrast_summary_table;

% Error degrees of freedom (scalar summary; per-voxel df lives in obj.df)
if ~isempty(out.df) && ~isempty(out.df.dat)
    obj.dfe = double(median(out.df.dat(:)));
else
    obj.dfe = size(X, 1) - size(X, 2);
end
obj.t.dfe = obj.dfe;

if ~isempty(out.contrast_estimates)
    obj.contrast_estimates = out.contrast_estimates;
    obj.contrast_t         = out.contrast_t;
    obj.contrast_t.dfe     = obj.dfe;
end

% Stash residuals so run_diagnostics() can compute Cook's distance. If the caller
% did not ask to keep them, they are cleared after diagnostics (below).
if ~isempty(out.residuals)
    obj.residuals = out.residuals;
end

% Sync names back from regress (it may have generated/added them)
if isempty(obj.design)
    obj.regressor_names_direct = out.regressor_names;
end
if ~isempty(out.contrast_names)
    obj.contrast_names = out.contrast_names;
end

% -------------------------------------------------------------------------
% Diagnostics and provenance
% -------------------------------------------------------------------------
obj.warnings = [obj.warnings(:); out.warnings(:)]';

obj.fit_parameters = struct( ...
    'robust', do_robust, ...
    'ar_order', ar_order, ...
    'is_timeseries', obj.is_timeseries, ...
    'pthresh', pthresh, ...
    'thresh_type', thresh_type, ...
    'do_resid', do_resid);

% Compute the full diagnostic set (adds cVIF, Cook's distance, condition
% number, collinearity report; uses canonical VIF/cVIF rather than getvif)
obj = run_diagnostics(obj, 'noverbose');

% Drop the (potentially large) residual maps unless the caller kept them
if ~do_resid
    obj.residuals = [];
end

obj.history{end + 1} = sprintf('fit: regress (%s%s), p<%g %s, dfe=%g', ...
    local_iif(do_robust, 'robust', 'OLS'), ...
    local_iif(ar_order > 0, sprintf(', AR(%d)', ar_order), ''), ...
    pthresh, thresh_type, obj.dfe);

if doverbose
    fprintf('\n  glm_map fit complete: %d regressors, %d contrasts, dfe = %g.\n', ...
        obj.num_regressors, obj.num_contrasts, obj.dfe);
    if ~isempty(obj.warnings)
        fprintf('  %d warning(s); see obj.warnings.\n', numel(obj.warnings));
    end
end

end % fit


% =====================================================================
function s = local_iif(tf, a, b)
% Inline if: return a if tf else b.
if tf, s = a; else, s = b; end
end % local_iif


% =====================================================================
function [X, names, C] = local_prepare_design(obj, X, names, C)
% Ensure the design has a single overall intercept as the LAST column (unless
% it already carries per-run baseline columns), keep regressor names aligned,
% and pad/reorder the contrast matrix C ([regressors x contrasts]) to match.
% Uses intercept.m semantics ('which'/'end'): an existing all-ones intercept is
% moved to the end (never duplicated); one is added if none is present.
k = size(X, 2);

% Pad names to the number of design columns
if ~iscell(names), names = {}; end
names = names(:)';
while numel(names) < k, names{end + 1} = sprintf('R%d', numel(names) + 1); end %#ok<AGROW>
names = names(1:k);

wh_allones = intercept(X, 'which');          % columns that are exactly all ones
n_block    = sum(obj.wh_intercept);          % per-run / baseline intercept columns

if numel(wh_allones) == 1
    % A single overall intercept: enforce it as the last column.
    if wh_allones ~= k
        order = [setdiff(1:k, wh_allones), wh_allones];
        X = X(:, order);
        names = names(order);
        if ~isempty(C) && size(C, 1) == k, C = C(order, :); end
    end

elseif isempty(wh_allones) && n_block == 0
    % No intercept of any kind: add one at the end.
    X = [X, ones(size(X, 1), 1)];
    names{end + 1} = 'Intercept';
    if ~isempty(C), C = [C; zeros(1, size(C, 2))]; end
end
% (multiple all-ones columns, or per-run block intercepts present: leave as-is)

% Pad contrast rows with 0 weight for any trailing (e.g. intercept) columns
if ~isempty(C)
    if size(C, 1) < size(X, 2)
        C = [C; zeros(size(X, 2) - size(C, 1), size(C, 2))];
    elseif size(C, 1) > size(X, 2)
        error('glm_map:ContrastSize', ...
            'Contrast matrix has %d rows but the design has %d columns (incl. intercept).', ...
            size(C, 1), size(X, 2));
    end
end

end % local_prepare_design
