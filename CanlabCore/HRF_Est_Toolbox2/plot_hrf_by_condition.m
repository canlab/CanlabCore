function ax = plot_hrf_by_condition(results, varargin)
%PLOT_HRF_BY_CONDITION Plot one subject/run HRF curves by condition.
%
% ax = plot_hrf_by_condition(results, 'Model', 'sfir')
% ax = plot_hrf_by_condition(results, 'Model', {'fir','sfir','canonical'})
% ax = plot_hrf_by_condition(results, 'Model', 'mapscore', 'Signature', 'sig_all_NPS')
%
% PlotType='fit' plots fitted HRF/model curves. If the fit contains .se and
% .p fields, subject/run-level uncertainty is shaded and significant bins can
% be marked. PlotType='trialmean' plots raw event-locked trial means with SEM
% from repeated events in the run.

p = inputParser;
p.addRequired('results', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('Conditions', [], @(x) isempty(x) || isnumeric(x) || iscell(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x));
p.addParameter('Signature', '', @(x) ischar(x) || isstring(x));
p.addParameter('PlotType', 'fit', @(x) ischar(x) || isstring(x));
p.addParameter('ShowSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('ShowP', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('RecomputeSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('SEAlpha', 0.18, @(x) isscalar(x) && x >= 0 && x <= 1);
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('BaselineSeconds', 0, @(x) isscalar(x) && x >= 0);
p.addParameter('LineWidth', 1.8, @(x) isscalar(x) && x > 0);
p.parse(results, varargin{:});
opts = p.Results;

model_names = local_model_names(opts.Model);
plot_type = lower(char(opts.PlotType));
condition_specs = local_condition_specs(results, opts);

figure; ax = axes; hold(ax, 'on');
switch plot_type
    case 'fit'
        local_plot_fit(ax, results, model_names, condition_specs, opts);
    case {'trialmean', 'trial_mean', 'trials'}
        local_plot_trial_mean(ax, results, condition_specs, opts);
    otherwise
        error('Unknown PlotType: %s. Use ''fit'' or ''trialmean''.', char(opts.PlotType));
end
hline(0, 'k-');
end

function local_plot_fit(ax, results, model_names, condition_specs, opts)
[fit_struct, source_label, sig_label] = local_fit_struct(results, opts.Signature);
colors = lines(numel(condition_specs));
line_styles = {'-', '--', ':', '-.'};
legend_labels = {};
plotted_se = false(1, numel(model_names));
used_models = {};
for m = 1:numel(model_names)
    model_name = model_names{m};
    if ~isfield(fit_struct, model_name)
        error('Model %s not available in selected fit structure.', model_name);
    end

    fit = local_fit_with_uncertainty(results, fit_struct.(model_name), model_name, opts);
    y_mat = fit.hrf;
    x = local_fit_time(fit, results, size(y_mat, 1));
    has_se = isfield(fit, 'se') && ~isempty(fit.se) && all(size(fit.se) == size(y_mat));
    has_p = isfield(fit, 'p') && ~isempty(fit.p) && all(size(fit.p) == size(y_mat));
    style = line_styles{mod(m - 1, numel(line_styles)) + 1};
    used_models{end + 1} = model_name; %#ok<AGROW>

    for k = 1:numel(condition_specs)
        cidx = condition_specs(k).indices;
        y = local_mean_omitnan(y_mat(:, cidx), 2);
        if logical(opts.ShowSE) && has_se
            se = local_combine_condition_se(y_mat(:, cidx), fit.se(:, cidx));
            if any(~isnan(se))
                fill(ax, [x; flipud(x)], [y + se; flipud(y - se)], colors(k, :), ...
                    'FaceAlpha', opts.SEAlpha ./ max(numel(model_names), 1), ...
                    'EdgeColor', 'none', 'HandleVisibility', 'off');
                plotted_se(m) = true;
            end
        end
        plot(ax, x, y, 'LineWidth', opts.LineWidth, 'Color', colors(k, :), ...
            'LineStyle', style);
        if logical(opts.ShowP) && has_p && isscalar(cidx)
            sig = fit.p(:, cidx) < opts.Alpha;
            if any(sig)
                yrange = max(y_mat(:)) - min(y_mat(:));
                if yrange == 0 || isnan(yrange), yrange = 1; end
                ymark = y(sig) + 0.04 .* yrange;
                plot(ax, x(sig), ymark, '.', 'Color', colors(k, :), ...
                    'MarkerSize', 12, 'HandleVisibility', 'off');
            end
        end
        if isscalar(model_names)
            legend_labels{end + 1} = condition_specs(k).display_label; %#ok<AGROW>
        else
            legend_labels{end + 1} = sprintf('%s | %s', condition_specs(k).display_label, model_name); %#ok<AGROW>
        end
    end
end

legend(ax, format_strings_for_legend(legend_labels), 'Interpreter', 'none');
title(ax, local_fit_title(used_models, source_label, sig_label, any(plotted_se)), 'Interpreter', 'none');
xlabel(ax, 'Seconds after event onset');
ylabel(ax, local_ylabel(used_models, source_label));
end

function local_plot_trial_mean(ax, results, condition_specs, opts)
if ~isfield(results, 'events') || isempty(results.events)
    error('PlotType=''trialmean'' requires results.events.');
end
if ~isfield(results, 'settings') || ~isfield(results.settings, 'TR') || ...
        ~isfield(results.settings, 'window_seconds')
    error('PlotType=''trialmean'' requires results.settings.TR and results.settings.window_seconds.');
end
[tc, source_label] = local_trial_timeseries(results, opts.Signature);

colors = lines(numel(condition_specs));
legend_labels = cell(1, numel(condition_specs));
for k = 1:numel(condition_specs)
    avg = hrf_average_condition_trials(tc, results.events, condition_specs(k).matched_conditions, ...
        results.settings.TR, results.settings.window_seconds, ...
        'BaselineSeconds', opts.BaselineSeconds);
    x = avg.time;
    y = avg.mean;
    if logical(opts.ShowSE)
        fill(ax, [x; flipud(x)], [y + avg.sem; flipud(y - avg.sem)], colors(k, :), ...
            'FaceAlpha', opts.SEAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    plot(ax, x, y, 'LineWidth', opts.LineWidth, 'Color', colors(k, :));
    legend_labels{k} = sprintf('%s (n=%d)', condition_specs(k).display_label, avg.n_trials);
end

legend(ax, format_strings_for_legend(legend_labels), 'Interpreter', 'none');
title(ax, sprintf('model=trialmean, source=%s, ribbon=within-run trial SEM', source_label), ...
    'Interpreter', 'none');
xlabel(ax, 'Seconds after event onset');
ylabel(ax, 'Observed signal');
end

function [tc, source_label] = local_trial_timeseries(results, signature)
source_label = local_source_label(results);
sig = char(signature);
if ~isempty(sig)
    sig_field = local_timeseries_signature_field(results, sig);
    if ~isempty(sig_field)
        tc = results.timeseries_by_signature.(sig_field);
        source_label = sprintf('%s, %s', source_label, sig);
        return
    end
    if isfield(results, 'signature_meta') && isfield(results.signature_meta, 'selected_signature') && ...
            strcmp(char(results.signature_meta.selected_signature), sig) && ...
            isfield(results, 'timeseries') && ~isempty(results.timeseries)
        tc = results.timeseries;
        source_label = sprintf('%s, %s', source_label, sig);
        return
    end
    error(['PlotType=''trialmean'' requested Signature %s, but matching time series were not stored. ' ...
        'Rerun run_hrf_pipeline so results.timeseries_by_signature is saved, or omit Signature for the selected time series.'], sig);
end

if isfield(results, 'timeseries') && ~isempty(results.timeseries)
    tc = results.timeseries;
else
    error(['PlotType=''trialmean'' requires results.timeseries or results.timeseries_by_signature. ' ...
        'Map-score-only studies rebuilt from CSVs do not contain event-level time series.']);
end
end

function condition_specs = local_condition_specs(results, opts)
if ~isfield(results, 'conditions') || isempty(results.conditions)
    error('results.conditions is required.');
end
if ~isempty(opts.Conditions)
    cond_spec = opts.Conditions;
elseif ~isempty(char(opts.Condition))
    cond_spec = {char(opts.Condition)};
else
    cond_spec = 1:numel(results.conditions);
end
condition_specs = hrf_resolve_condition_patterns(results.conditions, cond_spec, 'DefaultMode', 'each');
condition_specs = local_apply_condition_group_labels(results, condition_specs);
end

function condition_specs = local_apply_condition_group_labels(results, condition_specs)
if ~isfield(results, 'condition_groups') || isempty(results.condition_groups)
    return
end
groups = results.condition_groups;
for i = 1:numel(condition_specs)
    if numel(condition_specs(i).indices) ~= 1
        continue
    end
    idx = condition_specs(i).indices(1);
    if idx <= numel(groups) && isfield(groups, 'display_label')
        condition_specs(i).display_label = groups(idx).display_label;
        if isfield(groups, 'matched_conditions')
            condition_specs(i).matched_conditions = groups(idx).matched_conditions;
        end
    end
end
end

function [fit_struct, source_label, sig_label] = local_fit_struct(results, signature)
source_label = local_source_label(results);
sig_label = char(signature);
if ~isempty(sig_label) && isfield(results, 'fits_by_signature')
    sig_field = local_signature_field(results, sig_label);
    if isempty(sig_field)
        error('Signature %s not found in results.fits_by_signature.', sig_label);
    end
    fit_struct = results.fits_by_signature.(sig_field);
    source_label = sprintf('%s, %s', source_label, sig_label);
elseif isfield(results, 'fits')
    fit_struct = results.fits;
else
    error('results must contain .fits or selected .fits_by_signature.');
end
end

function fit = local_fit_with_uncertainty(results, fit, model_name, opts)
if local_has_fit_se(fit) || ~logical(opts.ShowSE) || ~logical(opts.RecomputeSE)
    return
end
if ~ismember(lower(model_name), {'fir', 'sfir', 'canonical', 'spline'})
    return
end
if ~isempty(char(opts.Signature)) && ~local_is_selected_signature(results, char(opts.Signature))
    return
end
if ~isfield(results, 'timeseries') || isempty(results.timeseries) || ...
        ~isfield(results, 'stick_functions') || isempty(results.stick_functions) || ...
        ~isfield(results, 'settings') || ~isfield(results.settings, 'TR') || ...
        ~isfield(results.settings, 'window_seconds')
    return
end

try
    refit = hrf_fit_all_models(results.timeseries, results.settings.TR, ...
        results.stick_functions, results.settings.window_seconds, {model_name}, ...
        'DependencyPolicy', 'skip');
    if isfield(refit, model_name) && local_has_fit_se(refit.(model_name))
        fit = refit.(model_name);
        fit.uncertainty_source = sprintf('%s; recomputed for plotting', fit.uncertainty_source);
    end
catch
end
end

function tf = local_has_fit_se(fit)
tf = isfield(fit, 'se') && ~isempty(fit.se) && isfield(fit, 'hrf') && ...
    all(size(fit.se) == size(fit.hrf));
end

function tf = local_is_selected_signature(results, signature)
tf = false;
if isempty(signature)
    tf = true;
    return
end
if isfield(results, 'signature_meta') && isfield(results.signature_meta, 'selected_signature')
    tf = strcmp(char(results.signature_meta.selected_signature), signature);
end
end

function se = local_combine_condition_se(Y, SE)
n_cond = size(Y, 2);
has_model_se = any(~isnan(SE(:)));
if has_model_se
    valid = ~isnan(SE);
    SE(~valid) = 0;
    n = max(sum(valid, 2), 1);
    se = sqrt(sum(SE .^ 2, 2)) ./ n;
else
    se = nan(size(Y, 1), 1);
end
if n_cond > 1
    condition_sem = local_sem_omitnan(Y, 2);
    if has_model_se
        se = sqrt(se .^ 2 + condition_sem .^ 2);
    else
        se = condition_sem;
    end
end
end

function m = local_mean_omitnan(X, dim)
if nargin < 2, dim = 1; end
valid = ~isnan(X);
den = sum(valid, dim);
X(~valid) = 0;
m = sum(X, dim) ./ den;
m(den == 0) = NaN;
end

function se = local_sem_omitnan(X, dim)
if nargin < 2, dim = 1; end
mu = local_mean_omitnan(X, dim);
if dim == 1
    centered = X - repmat(mu, size(X, 1), 1);
else
    centered = X - repmat(mu, 1, size(X, 2));
end
centered(isnan(centered)) = 0;
n = sum(~isnan(X), dim);
sd = sqrt(sum(centered .^ 2, dim) ./ max(n - 1, 1));
se = sd ./ sqrt(n);
se(n == 0 | n == 1) = NaN;
end

function x = local_fit_time(fit, results, n)
if isfield(fit, 'time') && numel(fit.time) == n
    x = fit.time(:);
elseif isfield(fit, 'lag_seconds') && numel(fit.lag_seconds) == n
    x = fit.lag_seconds(:);
elseif isfield(results, 'settings') && isfield(results.settings, 'TR')
    x = 1 + (0:n - 1)' .* results.settings.TR;
else
    x = (1:n)';
end
end

function source_label = local_source_label(results)
source_label = 'selected signal';
if isfield(results, 'settings') && isfield(results.settings, 'signal_source')
    source_label = char(results.settings.signal_source);
elseif isfield(results, 'signature_meta') && isfield(results.signature_meta, 'signal_source')
    source_label = char(results.signature_meta.signal_source);
end
if isfield(results, 'settings') && isfield(results.settings, 'mapscore_object')
    source_label = sprintf('%s %s map scores', char(results.settings.mapscore_object), source_label);
elseif isfield(results, 'signature_meta') && isfield(results.signature_meta, 'object')
    source_label = sprintf('%s map scores', char(results.signature_meta.object));
end
end

function ttl = local_fit_title(model_names, source_label, sig_label, has_se)
if isempty(sig_label)
    sig_part = '';
else
    sig_part = sprintf(', score/signature=%s', sig_label);
end
if has_se
    se_part = ', ribbon=within-run SE';
else
    se_part = ', ribbon=none (SE unavailable)';
end
ttl = sprintf('model=%s, source=%s%s%s', strjoin(model_names, ' + '), source_label, sig_part, se_part);
end

function ylab = local_ylabel(model_names, source_label)
if any(strcmpi(model_names, 'mapscore')) || contains(lower(source_label), 'map score')
    ylab = 'Pattern expression / map score';
else
    ylab = 'Fitted response amplitude';
end
end

function sig_field = local_signature_field(results, sig)
sig_field = '';
if isfield(results.fits_by_signature, sig)
    sig_field = sig;
    return
end

candidate = matlab.lang.makeValidName(sig);
if isfield(results.fits_by_signature, candidate)
    sig_field = candidate;
    return
end

if isfield(results, 'signature_meta') && isfield(results.signature_meta, 'selected_signatures') && ...
        isfield(results.signature_meta, 'selected_signature_fields')
    names = cellstr(string(results.signature_meta.selected_signatures));
    fields = cellstr(string(results.signature_meta.selected_signature_fields));
    idx = find(strcmp(names, sig), 1);
    if ~isempty(idx) && idx <= numel(fields) && isfield(results.fits_by_signature, fields{idx})
        sig_field = fields{idx};
    end
end
end

function model_names = local_model_names(model_input)
model_names = cellstr(string(model_input));
model_names = cellfun(@(s) lower(strtrim(s)), model_names, 'UniformOutput', false);
model_names = model_names(~cellfun(@isempty, model_names));
if isempty(model_names)
    error('Model must contain at least one model name.');
end
model_names = unique(model_names, 'stable');
end

function sig_field = local_timeseries_signature_field(results, sig)
sig_field = '';
if ~isfield(results, 'timeseries_by_signature') || isempty(results.timeseries_by_signature)
    return
end
if isfield(results.timeseries_by_signature, sig)
    sig_field = sig;
    return
end

candidate = matlab.lang.makeValidName(sig);
if isfield(results.timeseries_by_signature, candidate)
    sig_field = candidate;
    return
end

if isfield(results, 'signature_meta') && isfield(results.signature_meta, 'selected_signatures') && ...
        isfield(results.signature_meta, 'selected_signature_fields')
    names = cellstr(string(results.signature_meta.selected_signatures));
    fields = cellstr(string(results.signature_meta.selected_signature_fields));
    idx = find(strcmp(names, sig), 1);
    if ~isempty(idx) && idx <= numel(fields) && isfield(results.timeseries_by_signature, fields{idx})
        sig_field = fields{idx};
    end
end
end
