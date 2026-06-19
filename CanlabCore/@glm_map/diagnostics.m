function obj = diagnostics(obj, varargin)
% Compute and report design diagnostics for a glm_map object.
%
% Evaluates the conditioning of the design matrix X (and contrasts C):
% variance inflation factors (VIF) per regressor, contrast VIFs (cVIF),
% per-observation leverage and Cook's distance, condition number, rank
% deficiency, and a redundant/near-collinear column report. VIFs/cVIFs and the
% condition number are computed both for the full design and for the
% regressors of interest only (excluding nuisance covariates), so that the
% influence of nuisance covariates on the of-interest estimates is visible.
% By default a narrative report is printed with interpretation and warnings.
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
%        Regressor roles (of interest vs nuisance vs intercept) are read from
%        obj.wh_interest / obj.wh_nuisance / obj.wh_intercept.
%
% :Optional Inputs:
%
%   **'doverbose' / 'noverbose':**
%        Print the narrative report (default true).
%
%   **'vif_thresh', [value]:**
%        VIF/cVIF warning threshold (default 4).
%
%   **'cond_thresh', [value]:**
%        Condition-number warning threshold (default 30; > 100 = severe).
%
% :Outputs:
%
%   **obj:**
%        glm_map with obj.diagnostics populated. Fields include
%        Variance_inflation_factors, Contrast_variance_inflation_factors,
%        Leverages, Cooks_distance, condition_number, rank_deficient,
%        collinearity_report (all for the FULL design), plus the
%        *_interest_only counterparts, and warnings appended.
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
%    2026 - Operates on the design only (no fit required), so it can be used to
%    screen a design before fitting. Cook's distance additionally needs
%    residuals and is computed only when those are present.
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
doverbose = ~any(strcmpi(varargin, 'noverbose'));

vif_thresh = 4;
wh = find(strcmpi(varargin, 'vif_thresh'));
if ~isempty(wh), vif_thresh = varargin{wh(1) + 1}; end

cond_thresh = 30;
wh = find(strcmpi(varargin, 'cond_thresh'));
if ~isempty(wh), cond_thresh = varargin{wh(1) + 1}; end

X = obj.X;
if isempty(X)
    error('glm_map:NoDesign', 'No design matrix available (obj.X is empty). Build or supply a design first.');
end

mywarnings = {};

d = obj.diagnostics;
if ~isstruct(d), d = struct(); end
d.vif_threshold = vif_thresh;
d.cond_threshold = cond_thresh;

% Regressor roles
whI = obj.wh_interest(:)';
whN = obj.wh_nuisance(:)';
whB = obj.wh_intercept(:)';
wh_io = whI | whB;                      % interest + intercept (model without nuisance)
has_nuisance = any(whN);

% -------------------------------------------------------------------------
% Variance inflation factors -- full design
% -------------------------------------------------------------------------
d.Variance_inflation_factors = VIF(X);

if any(d.Variance_inflation_factors > vif_thresh)
    mywarnings{end + 1} = sprintf(['Design multicollinearity: %d regressor(s) have VIF > %g. ' ...
        'Check obj.diagnostics.Variance_inflation_factors and obj.regressor_names.'], ...
        sum(d.Variance_inflation_factors > vif_thresh), vif_thresh);
end

% Contrast VIFs -- full design
d.Contrast_variance_inflation_factors = [];
if ~isempty(obj.contrasts)
    if size(obj.contrasts, 1) ~= size(X, 2)
        mywarnings{end + 1} = sprintf(['Contrast matrix has %d rows but design has %d regressors; ' ...
            'skipping contrast VIFs.'], size(obj.contrasts, 1), size(X, 2));
    else
        d.Contrast_variance_inflation_factors = cVIF(X, obj.contrasts');
        if any(d.Contrast_variance_inflation_factors > vif_thresh)
            mywarnings{end + 1} = sprintf('%d contrast(s) have cVIF > %g (full design).', ...
                sum(d.Contrast_variance_inflation_factors > vif_thresh), vif_thresh);
        end
    end
end

% -------------------------------------------------------------------------
% Variance inflation factors -- regressors of interest only (no nuisance)
% Stored over the interest+intercept columns (see d.wh_interest_only_columns).
% -------------------------------------------------------------------------
d.Variance_inflation_factors_interest_only = [];
d.Contrast_variance_inflation_factors_interest_only = [];
d.wh_interest_only_columns = find(wh_io);

if has_nuisance && any(wh_io)
    Xio = X(:, wh_io);
    d.Variance_inflation_factors_interest_only = VIF(Xio);

    if ~isempty(obj.contrasts) && size(obj.contrasts, 1) == size(X, 2)
        Cio = obj.contrasts(wh_io, :);
        keepcon = any(Cio ~= 0, 1);                 % contrasts that survive restriction
        if any(keepcon)
            d.Contrast_variance_inflation_factors_interest_only = cVIF(Xio, Cio(:, keepcon)');
        end
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
% Cook's distance (per observation; requires residuals + sigma)
%   D_i = mean_v [ r_iv^2 / (k * s_v^2) ] * h_i / (1 - h_i)^2
% -------------------------------------------------------------------------
d.Cooks_distance = [];
have_resid = ~isempty(obj.residuals) && isa(obj.residuals, 'fmri_data') && ~isempty(obj.residuals.dat);
have_sigma = ~isempty(obj.sigma)     && isa(obj.sigma, 'fmri_data')     && ~isempty(obj.sigma.dat);

if have_resid && have_sigma
    h = d.Leverages(:);
    k = size(X, 2);
    r = obj.residuals.dat;
    s2 = obj.sigma.dat(:) .^ 2;
    if size(r, 2) == numel(h)
        good = isfinite(s2) & s2 > 0;
        if any(good)
            msr  = mean( (r(good, :) .^ 2) ./ s2(good), 1 );
            hfac = (h ./ (1 - h) .^ 2)';
            d.Cooks_distance = (msr / k) .* hfac;
            if any(d.Cooks_distance > 1)
                mywarnings{end + 1} = sprintf(['%d observation(s) have Cook''s distance > 1 ' ...
                    '(highly influential). Check obj.diagnostics.Cooks_distance.'], sum(d.Cooks_distance > 1));
            end
        end
    end
end

% -------------------------------------------------------------------------
% Conditioning / rank (full and interest-only)
%
% The condition number is computed after scaling each column to unit L2 norm
% (the Belsley scaled condition index). This makes it scale-invariant -- like
% VIF -- so it reflects true collinearity rather than differences in regressor
% magnitude. Without scaling, a constant intercept (large norm) alongside
% small-amplitude event regressors produces a large condition number that is a
% pure scaling artifact, even when the regressors are nearly orthogonal.
% -------------------------------------------------------------------------
d.condition_number = local_scaled_cond(X);
d.condition_number_interest_only = [];
if has_nuisance && any(wh_io), d.condition_number_interest_only = local_scaled_cond(X(:, wh_io)); end
d.rank_deficient = rank(X) < size(X, 2);

if d.rank_deficient
    mywarnings{end + 1} = 'Design matrix X is rank deficient (rank(X) < number of regressors).';
end
if d.condition_number > 100
    mywarnings{end + 1} = sprintf('Design is severely ill-conditioned (condition number %.0f > 100).', d.condition_number);
elseif d.condition_number > cond_thresh
    mywarnings{end + 1} = sprintf('Design is moderately ill-conditioned (condition number %.0f > %g).', d.condition_number, cond_thresh);
end

% Does adding nuisance covariates substantially inflate of-interest VIFs?
infl = local_inflation(d, whI, wh_io);
if ~isempty(infl) && infl.max_ratio >= 1.5
    mywarnings{end + 1} = sprintf(['Nuisance covariates inflate of-interest VIFs by up to %.1fx ' ...
        '(max VIF %.2f vs %.2f without nuisance); events of interest may correlate with nuisance variables.'], ...
        infl.max_ratio, infl.max_full, infl.max_io);
end

% -------------------------------------------------------------------------
% Design efficiency for contrasts (calcEfficiency). Uses the entered
% contrasts, or an orthogonal set spanning the regressors of interest when
% none are entered. Higher efficiency = lower contrast-estimate variance.
% -------------------------------------------------------------------------
[d.efficiency, d.efficiency_per_contrast, d.efficiency_contrast_names, ...
    d.efficiency_contrast_source, eff_note] = local_efficiency(obj, X, whI);

% -------------------------------------------------------------------------
% Redundant / near-collinear column report
% -------------------------------------------------------------------------
report = struct();
report.vif_threshold    = vif_thresh;
report.high_vif_columns = find(d.Variance_inflation_factors > vif_thresh);

ncol = size(X, 2);
dup_pairs = [];
for a = 1:ncol - 1
    for b = a + 1:ncol
        if isequal(X(:, a), X(:, b)), dup_pairs(end + 1, :) = [a b]; end %#ok<AGROW>
    end
end
report.duplicate_column_pairs = dup_pairs;
if ~isempty(dup_pairs)
    mywarnings{end + 1} = sprintf('%d pair(s) of identical design columns detected (see obj.diagnostics.collinearity_report).', size(dup_pairs, 1));
end

R = corrcoef(X);
R(logical(eye(ncol))) = 0;
[ia, ib] = find(triu(abs(R) > 0.95, 1));
report.high_correlation_pairs = [ia ib];

d.collinearity_report = report;

% -------------------------------------------------------------------------
% Store
% -------------------------------------------------------------------------
obj.diagnostics = d;
obj.warnings = [obj.warnings(:); mywarnings(:)]';
obj.history{end + 1} = 'diagnostics: VIF/cVIF (full + interest-only), leverage, Cook''s D, conditioning, collinearity';

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------
if doverbose
    local_report(obj, d, infl, vif_thresh, cond_thresh, mywarnings, eff_note);
end

end % diagnostics


% =========================================================================
% Local helpers
% =========================================================================
function [eff, eff_vec, cnames, src, note] = local_efficiency(obj, X, whI)
% Design efficiency for contrasts via calcEfficiency. Uses entered contrasts,
% else an orthogonal set spanning the regressors of interest. Returns [] for
% eff and an informative note when efficiency cannot be computed.
[eff, eff_vec, cnames, src, note] = deal([], [], {}, '', '');

Crows = [];
if ~isempty(obj.contrasts) && size(obj.contrasts, 1) == size(X, 2)
    Crows  = obj.contrasts';                      % [n_con x n_reg]
    cnames = obj.contrast_names;
    src    = 'entered contrasts';

elseif any(whI)
    whIidx = find(whI);
    if numel(whIidx) >= 2
        Cint = create_orthogonal_contrast_set(numel(whIidx));   % dispatches to the function
        Crows = zeros(size(Cint, 1), size(X, 2));
        Crows(:, whIidx) = Cint;
        cnames = arrayfun(@(i) sprintf('OrthC%d', i), 1:size(Crows, 1), 'UniformOutput', false);
        src = 'auto orthogonal set spanning the regressors of interest';
    else
        note = 'Efficiency not computed: need >= 2 regressors of interest, or enter contrasts.';
        return
    end
else
    note = ['Efficiency not computed: no regressors of interest are defined and no contrasts ' ...
        'entered. Define regressors of interest (build an event design, or set ' ...
        'obj.nuisance_columns for direct designs) or add contrasts.'];
    return
end

try
    [eff, eff_vec] = calcEfficiency(ones(1, size(Crows, 1)), Crows, pinv(X), []);
    eff_vec = eff_vec(:)';
catch ME
    [eff, eff_vec, cnames, src] = deal([], [], {}, '');
    note = sprintf('Efficiency not computed (calcEfficiency error): %s', ME.message);
end
end


function c = local_scaled_cond(X)
% Condition number after scaling each column to unit L2 norm (scale-invariant
% Belsley condition index).
nrm = vecnorm(X, 2, 1);
nrm(nrm == 0) = 1;
c = cond(X ./ nrm);
end


function infl = local_inflation(d, whI, wh_io)
% Compare full-design VIFs of the of-interest regressors against their
% interest-only VIFs. Returns a struct with the max inflation ratio, or [].
infl = [];
if isempty(d.Variance_inflation_factors_interest_only), return, end

vfull = d.Variance_inflation_factors;
vio   = d.Variance_inflation_factors_interest_only;          % over wh_io columns
io_cols = find(wh_io);

% map each interest column to its position within the interest-only design
isI_in_io = whI(io_cols);                  % which io columns are of-interest
full_on_I = vfull(io_cols(isI_in_io));
io_on_I   = vio(isI_in_io);

valid = io_on_I > 0 & isfinite(io_on_I) & isfinite(full_on_I);
if ~any(valid), return, end

ratios = full_on_I(valid) ./ io_on_I(valid);
[infl.max_ratio, w] = max(ratios);
fv = full_on_I(valid); iv = io_on_I(valid);
infl.max_full = fv(w);
infl.max_io   = iv(w);
end


function local_report(obj, d, infl, vif_thresh, cond_thresh, mywarnings, eff_note)

line = repmat('-', 1, 70);
fprintf('\n  glm_map design diagnostics\n  %s\n', line);
if ~isempty(obj.analysis_name), fprintf('  Analysis: %s\n', obj.analysis_name); end
fprintf('  %d images x %d regressors  (%d of interest, %d nuisance, %d intercept)\n', ...
    obj.num_images, obj.num_regressors, sum(obj.wh_interest), sum(obj.wh_nuisance), sum(obj.wh_intercept));

% ---- VIFs (full design) ----
rn = obj.regressor_names;
fprintf('\n  Variance inflation factors (VIF), full design\n');
fprintf('  VIF >= 1 (1 = orthogonal). ~1-2 low, 2-4 moderate, >%g high collinearity, Inf = exact dependence.\n', vif_thresh);
fprintf('  %-28s %8s  role\n', 'Regressor', 'VIF');
for i = 1:numel(d.Variance_inflation_factors)
    nm = local_name(rn, i, 'R');
    flag = local_flag(d.Variance_inflation_factors(i) > vif_thresh, '  <-- high');
    fprintf('  %-28s %8.2f  %s%s\n', nm, d.Variance_inflation_factors(i), local_role(obj, i), flag);
end

% ---- Interest-only comparison ----
if ~isempty(d.Variance_inflation_factors_interest_only)
    fprintf('\n  VIF, regressors of interest only (nuisance covariates removed)\n');
    io_cols = d.wh_interest_only_columns;
    vio = d.Variance_inflation_factors_interest_only;
    for j = 1:numel(io_cols)
        if obj.wh_intercept(io_cols(j)), continue, end   % skip intercept in this listing
        nm = local_name(rn, io_cols(j), 'R');
        fprintf('  %-28s %8.2f  (full: %.2f)\n', nm, vio(j), d.Variance_inflation_factors(io_cols(j)));
    end
    if ~isempty(infl)
        if infl.max_ratio >= 1.5
            fprintf('  => Nuisance covariates inflate of-interest VIFs up to %.1fx (%.2f vs %.2f).\n', ...
                infl.max_ratio, infl.max_full, infl.max_io);
            fprintf('     This suggests events of interest are correlated with nuisance variables.\n');
        else
            fprintf('  => Nuisance covariates have little effect on of-interest VIFs (max %.1fx).\n', infl.max_ratio);
        end
    end
end

% ---- Contrast VIFs ----
if ~isempty(d.Contrast_variance_inflation_factors)
    fprintf('\n  Contrast VIFs (cVIF), full design  [same scale as VIF]\n');
    for i = 1:numel(d.Contrast_variance_inflation_factors)
        nm = local_name(obj.contrast_names, i, 'Con');
        flag = local_flag(d.Contrast_variance_inflation_factors(i) > vif_thresh, '  <-- high');
        fprintf('  %-28s %8.2f%s\n', nm, d.Contrast_variance_inflation_factors(i), flag);
    end
    if ~isempty(d.Contrast_variance_inflation_factors_interest_only)
        fprintf('  Interest-only cVIF: %s\n', mat2str(round(d.Contrast_variance_inflation_factors_interest_only, 2)));
    end
end

% ---- Conditioning ----
fprintf('\n  Conditioning\n');
fprintf('  Scaled condition number (columns at unit norm; scale-invariant, like VIF).\n');
fprintf('  >= 1; <%g OK, %g-100 moderate, >100 severe multicollinearity.\n', cond_thresh, cond_thresh);
fprintf('  condition number (full)         : %.1f%s\n', d.condition_number, ...
    local_flag(d.condition_number > cond_thresh, '  <-- elevated'));
if ~isempty(d.condition_number_interest_only)
    fprintf('  condition number (interest only): %.1f\n', d.condition_number_interest_only);
end
fprintf('  rank deficient                  : %s\n', local_yn(d.rank_deficient));

% ---- Cook's distance ----
if ~isempty(d.Cooks_distance)
    [maxcd, wcd] = max(d.Cooks_distance);
    fprintf('\n  Influence: Cook''s distance per observation (image)\n');
    fprintf('  > 1 = highly influential (rule of thumb; 4/n = %.3f is a stricter cutoff).\n', 4 / obj.num_images);
    fprintf('  max Cook''s distance: %.3f (observation %d)%s\n', maxcd, wcd, local_flag(maxcd > 1, '  <-- influential'));
end

% ---- Design efficiency ----
fprintf('\n  Design efficiency for contrasts (calcEfficiency)\n');
fprintf('  Relative measure (no absolute scale); higher = lower contrast-estimate\n');
fprintf('  variance = more efficient. Compare designs/contrasts on the same data.\n');
if ~isempty(d.efficiency)
    fprintf('  Contrasts: %s\n', d.efficiency_contrast_source);
    for i = 1:numel(d.efficiency_per_contrast)
        nm = local_name(d.efficiency_contrast_names, i, 'Con');
        fprintf('  %-28s %12.4g\n', nm, d.efficiency_per_contrast(i));
    end
    fprintf('  %-28s %12.4g\n', 'overall (mean)', d.efficiency);
else
    fprintf('  %s\n', eff_note);
end

% ---- Warnings ----
if ~isempty(mywarnings)
    fprintf('\n  %d warning(s):\n', numel(mywarnings));
    for i = 1:numel(mywarnings), fprintf('    - %s\n', mywarnings{i}); end
else
    fprintf('\n  No diagnostic warnings.\n');
end
fprintf('  %s\n\n', line);

end


function r = local_role(obj, i)
if obj.wh_intercept(i), r = 'intercept';
elseif obj.wh_nuisance(i), r = 'nuisance';
else, r = 'of interest';
end
end

function nm = local_name(names, idx, prefix)
if iscell(names) && idx <= numel(names) && ~isempty(names{idx}), nm = names{idx};
else, nm = sprintf('%s%d', prefix, idx);
end
end

function s = local_flag(tf, str), if tf, s = str; else, s = ''; end, end
function s = local_yn(tf), if tf, s = 'yes'; else, s = 'no'; end, end
