function summary(obj)
% Print a narrative summary report for a glm_map object.
%
% Unlike disp/display (which list the object's properties), summary prints a
% human-readable report: the analysis name, a description of the model
% (level and input variables), the design diagnostics (VIF, conditioning,
% leverage, Cook's distance, collinearity), and -- if the model has been
% fitted -- a summary of the results, including how many result maps there
% are, the statistical threshold, and the number of significant voxels at
% that threshold for each predictor and contrast.
%
% :Usage:
% ::
%
%     summary(obj)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object (fitted or unfitted).
%
% :See also:
%   - glm_map, glm_map.diagnostics, glm_map.disp
%
% ..
%    Programmers' notes:
%    2026 - Rewritten as a narrative report (was a thin wrapper over disp).
% ..

line = repmat('=', 1, 68);
fprintf('\n%s\n  glm_map analysis summary\n%s\n', line, line);

% ---------------------------------------------------------------------
% Analysis name
% ---------------------------------------------------------------------
if ~isempty(obj.analysis_name)
    fprintf('  Analysis: %s\n\n', obj.analysis_name);
end

% ---------------------------------------------------------------------
% Model description (level and input variables)
% ---------------------------------------------------------------------
switch obj.level
    case 1, levelstr = 'first-level (within-run / single-subject)';
    case 2, levelstr = 'second-level (group)';
    otherwise, levelstr = num2str(obj.level);
end

if ~isempty(obj.design)
    modestr = 'event/1st-level (HRF-convolved design)';
else
    modestr = 'direct design matrix';
end

fprintf('  Model\n  %s\n', repmat('-', 1, 64));
fprintf('  Level           : %s\n', levelstr);
fprintf('  Design          : %s', modestr);
if obj.is_timeseries, fprintf('  [timeseries; AR errors allowed]'); end
fprintf('\n');
fprintf('  Observations    : %d images\n', obj.num_images);
fprintf('  Regressors      : %d (%d of interest, %d nuisance, %d intercept)\n', ...
    obj.num_regressors, sum(obj.wh_interest), sum(obj.wh_nuisance), sum(obj.wh_intercept));

rn = obj.regressor_names;
if ~isempty(rn)
    for i = 1:numel(rn)
        nm = rn{i}; if isempty(nm), nm = sprintf('R%d', i); end
        fprintf('      %2d. %-26s [%s]\n', i, nm, local_role(obj, i));
    end
end

if obj.num_contrasts > 0
    fprintf('  Contrasts       : %d\n', obj.num_contrasts);
    cn = obj.contrast_names;
    for i = 1:obj.num_contrasts
        if i <= numel(cn) && ~isempty(cn{i}), nm = cn{i}; else, nm = sprintf('Con%d', i); end
        fprintf('      %2d. %s\n', i, nm);
    end
end
fprintf('\n');

% ---------------------------------------------------------------------
% Design diagnostics
% ---------------------------------------------------------------------
d = obj.diagnostics;
if isstruct(d) && ~isempty(fieldnames(d))
    fprintf('  Design diagnostics\n  %s\n', repmat('-', 1, 64));

    if isfield(d, 'Variance_inflation_factors') && ~isempty(d.Variance_inflation_factors)
        [maxvif, wv] = max(d.Variance_inflation_factors);
        nm = local_name(rn, wv, 'R');
        fprintf('  Max VIF         : %.2f (%s)%s\n', maxvif, nm, local_flag(maxvif > 4, '  high'));
    end

    if isfield(d, 'Contrast_variance_inflation_factors') && ~isempty(d.Contrast_variance_inflation_factors)
        [maxcv, wc] = max(d.Contrast_variance_inflation_factors);
        nm = local_name(obj.contrast_names, wc, 'Con');
        fprintf('  Max contrast VIF: %.2f (%s)%s\n', maxcv, nm, local_flag(maxcv > 4, '  high'));
    end

    if isfield(d, 'condition_number') && ~isempty(d.condition_number)
        fprintf('  Condition number: %.2f\n', d.condition_number);
    end

    if isfield(d, 'rank_deficient') && ~isempty(d.rank_deficient)
        fprintf('  Rank deficient  : %s\n', local_yn(d.rank_deficient));
    end

    if isfield(d, 'Leverages') && ~isempty(d.Leverages)
        nhi = sum(abs(zscore(d.Leverages(:))) >= 3);
        fprintf('  High leverage   : %d observation(s) with abs(z) >= 3\n', nhi);
    end

    if isfield(d, 'Cooks_distance') && ~isempty(d.Cooks_distance)
        [maxcd, wcd] = max(d.Cooks_distance);
        fprintf('  Max Cook''s dist : %.3f (observation %d)%s\n', maxcd, wcd, local_flag(maxcd > 1, '  influential'));
    end

    if isfield(d, 'collinearity_report') && isstruct(d.collinearity_report)
        cr = d.collinearity_report;
        ndup = 0; if isfield(cr, 'duplicate_column_pairs'), ndup = size(cr.duplicate_column_pairs, 1); end
        ncor = 0; if isfield(cr, 'high_correlation_pairs'),  ncor = size(cr.high_correlation_pairs, 1); end
        if ndup > 0 || ncor > 0
            fprintf('  Collinearity    : %d duplicate column pair(s), %d high-correlation pair(s)\n', ndup, ncor);
        end
    end
    fprintf('\n');
end

% ---------------------------------------------------------------------
% Results (if fitted)
% ---------------------------------------------------------------------
if obj.is_fitted
    fprintf('  Results\n  %s\n', repmat('-', 1, 64));

    % How many maps
    fprintf('  Maps            : betas [%d], t [%d]', local_nimg(obj.betas), local_nimg(obj.t));
    if ~isempty(obj.contrast_estimates)
        fprintf(', contrast estimates [%d], contrast t [%d]', ...
            local_nimg(obj.contrast_estimates), local_nimg(obj.contrast_t));
    end
    fprintf('\n');
    fprintf('  Error df        : %g\n', obj.dfe);

    % Statistical threshold (from the thresholded t map)
    fprintf('  Threshold       : %s\n', local_threshold_str(obj.t, obj.fit_parameters));

    % Significant voxels per predictor (t map)
    fprintf('  Significant voxels at threshold (per regressor):\n');
    local_print_sig(obj.t, rn, 'R');

    % Significant voxels per contrast (contrast t map)
    if ~isempty(obj.contrast_t)
        fprintf('  Significant voxels at threshold (per contrast):\n');
        local_print_sig(obj.contrast_t, obj.contrast_names, 'Con');
    end
    fprintf('\n');
else
    fprintf('  Results         : not fitted (run fit(obj, fmri_data_obj))\n\n');
end

% ---------------------------------------------------------------------
% Warnings
% ---------------------------------------------------------------------
if ~isempty(obj.warnings)
    fprintf('  Warnings (%d)\n  %s\n', numel(obj.warnings), repmat('-', 1, 64));
    for i = 1:numel(obj.warnings)
        fprintf('    - %s\n', obj.warnings{i});
    end
    fprintf('\n');
end

fprintf('%s\n\n', line);

end % summary


% =====================================================================
% Local helpers
% =====================================================================
function r = local_role(obj, i)
if obj.wh_intercept(i), r = 'intercept';
elseif obj.wh_nuisance(i), r = 'nuisance';
else, r = 'of interest';
end
end

function nm = local_name(names, idx, prefix)
if iscell(names) && idx <= numel(names) && ~isempty(names{idx})
    nm = names{idx};
else
    nm = sprintf('%s%d', prefix, idx);
end
end


function s = local_flag(tf, str)
if tf, s = str; else, s = ''; end
end


function s = local_yn(tf)
if tf, s = 'yes'; else, s = 'no'; end
end


function n = local_nimg(map)
if isempty(map) || ~isprop(map, 'dat') || isempty(map.dat), n = 0; else, n = size(map.dat, 2); end
end


function s = local_threshold_str(map, fit_parameters)
% Build a readable threshold description. Prefer the glm_map fit_parameters
% (clean p-value/type), then fall back to the statistic_image's own record.
s = '(unthresholded)';

if isstruct(fit_parameters) && isfield(fit_parameters, 'pthresh') && ~isempty(fit_parameters.pthresh)
    ttype = ''; if isfield(fit_parameters, 'thresh_type'), ttype = char(fit_parameters.thresh_type); end
    s = strtrim(sprintf('p < %g %s', fit_parameters.pthresh, ttype));

elseif ~isempty(map) && isprop(map, 'threshold') && ~isempty(map.threshold)
    thr = unique(map.threshold(:)');          % collapse duplicated pos/neg thresholds
    ttype = '';
    if isprop(map, 'thr_type') && ~isempty(map.thr_type)
        if iscell(map.thr_type), ttype = strjoin(cellfun(@char, map.thr_type(:)', 'UniformOutput', false), '/');
        else, ttype = char(map.thr_type); end
    end
    if isscalar(thr)
        s = strtrim(sprintf('p < %g %s', thr, ttype));
    else
        s = strtrim(sprintf('%s %s', mat2str(thr, 4), ttype));
    end
end
end


function local_print_sig(map, names, prefix)
% Print the number of significant voxels per image of a thresholded
% statistic_image map.
if isempty(map) || ~isprop(map, 'dat') || isempty(map.dat)
    fprintf('      (no map)\n');
    return
end

nimg = size(map.dat, 2);

% Prefer the .sig logical mask; fall back to nonzero thresholded values
if isprop(map, 'sig') && ~isempty(map.sig) && size(map.sig, 2) == nimg
    counts = sum(map.sig ~= 0, 1);
else
    counts = sum(map.dat ~= 0 & ~isnan(map.dat), 1);
end

for i = 1:nimg
    nm = local_name(names, i, prefix);
    fprintf('      %-26s %d voxels\n', nm, counts(i));
end
end
