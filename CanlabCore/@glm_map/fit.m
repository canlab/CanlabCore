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
%        Use an autoregressive error model of the given order. Only valid
%        when obj.is_timeseries is true.
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
% Assemble regress() argument list
% -------------------------------------------------------------------------
data.X = X;

% Threshold p-value must precede the 'unc'/'fdr' keyword (regress convention)
regress_args = {pthresh, thresh_type};

if do_robust,  regress_args{end + 1} = 'robust'; end
if ar_order > 0, regress_args = [regress_args, {'AR', ar_order}]; end
if do_resid,   regress_args{end + 1} = 'residual'; end

% Suppress regress's own orthviews unless explicitly requested
if ~do_display, regress_args{end + 1} = 'nodisplay'; end
if ~doverbose,  regress_args{end + 1} = 'noverbose'; end

% Regressor names (length must match number of design columns)
rn = obj.regressor_names;
if ~isempty(rn) && numel(rn) == size(X, 2)
    regress_args = [regress_args, {'names', rn}];
end

% Contrasts: regress expects C as [regressors x contrasts] (our storage)
if ~isempty(obj.contrasts)
    regress_args = [regress_args, {'C', obj.contrasts}];
    if ~isempty(obj.contrast_names)
        regress_args = [regress_args, {'contrast_names', obj.contrast_names}];
    end
end

if ~isempty(obj.analysis_name)
    regress_args = [regress_args, {'analysis_name', obj.analysis_name}];
end

% -------------------------------------------------------------------------
% Run the regression
% -------------------------------------------------------------------------
out = regress(data, regress_args{:});

% -------------------------------------------------------------------------
% Unpack result maps
% -------------------------------------------------------------------------
obj.betas = out.b;
obj.t     = out.t;
obj.sigma = out.sigma;

% Error degrees of freedom (scalar summary; per-voxel df lives in out.df)
if isfield(out, 'df') && ~isempty(out.df) && ~isempty(out.df.dat)
    obj.dfe = double(median(out.df.dat(:)));
else
    obj.dfe = size(X, 1) - size(X, 2);
end
obj.t.dfe = obj.dfe;

if isfield(out, 'contrast_images') && ~isempty(out.contrast_images)
    obj.contrast_estimates = out.contrast_images;
    obj.contrast_t         = out.con_t;
    obj.contrast_t.dfe     = obj.dfe;
end

if do_resid && isfield(out, 'resid')
    obj.residuals = out.resid;
end

% Sync names back from regress (it may have generated/added them)
if isempty(obj.design)
    obj.regressor_names_direct = out.variable_names;
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

% Compute the full diagnostic set (adds cVIF, condition number, collinearity
% report; uses canonical VIF/cVIF rather than regress's getvif)
obj = diagnostics(obj, 'noverbose');

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
