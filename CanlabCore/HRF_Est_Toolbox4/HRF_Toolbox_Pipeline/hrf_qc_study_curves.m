function qc = hrf_qc_study_curves(study, varargin)
%HRF_QC_STUDY_CURVES Inspect noisy/outlying HRF curves in a study struct.
%
% qc = hrf_qc_study_curves(study, 'Model', 'sfir', 'Condition', 'pain')
%
% This summarizes one condition/contrast curve per run or subject. For
% map-score studies, this is the right place to flag outlying runs/subjects;
% trial-level QC requires original run results with events/timeseries.

p = inputParser;
p.addRequired('study', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('SourceModel', '', @(x) ischar(x) || isstring(x));
p.addParameter('Signature', '', @(x) ischar(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConditionB', '', @(x) ischar(x) || isstring(x));
p.addParameter('Unit', 'run', @(x) ischar(x) || isstring(x));
p.addParameter('OutlierZThreshold', 4, @(x) isscalar(x) && x > 0);
p.addParameter('Weighting', 'huber', @(x) ischar(x) || isstring(x));
p.addParameter('Plot', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.parse(study, varargin{:});
opts = p.Results;

stats = hrf_time_unfolding_stats(study, ...
    'Model', opts.Model, ...
    'SourceModel', opts.SourceModel, ...
    'Signature', opts.Signature, ...
    'ConditionA', opts.Condition, ...
    'ConditionB', opts.ConditionB, ...
    'Unit', 'run', ...
    'MissingPolicy', opts.MissingPolicy);

Y_runs = stats.run_level_data;
run_subject_ids = stats.run_subject_ids(:);
if isfield(stats, 'run_labels') && numel(stats.run_labels) == numel(run_subject_ids)
    source_run_labels = cellstr(string(stats.run_labels(:)));
else
    source_run_labels = local_unique_run_labels(run_subject_ids);
end
unit = lower(char(opts.Unit));
switch unit
    case 'run'
        Y = Y_runs;
        labels = source_run_labels;
        subject_ids = run_subject_ids;
        run_count = ones(size(labels));
    case 'subject'
        subject_ids = unique(run_subject_ids, 'stable');
        Y = nan(numel(subject_ids), size(Y_runs, 2));
        run_count = zeros(numel(subject_ids), 1);
        for i = 1:numel(subject_ids)
            wh = strcmp(run_subject_ids, subject_ids{i});
            Y(i, :) = local_mean_omitnan(Y_runs(wh, :), 1);
            run_count(i) = sum(wh);
        end
        labels = subject_ids;
    otherwise
        error('Unknown Unit: %s. Use run or subject.', unit);
end

[max_abs_z, weights, is_outlier] = local_curve_weights(Y, opts.Weighting, opts.OutlierZThreshold);
rms_value = sqrt(local_mean_omitnan(Y .^ 2, 2));
peak_abs = local_max_omitnan(abs(Y), 2);
auc_abs = local_sum_omitnan(abs(Y), 2);
n_nan = sum(isnan(Y), 2);

qc = struct();
qc.table = table((1:size(Y, 1))', string(labels(:)), string(subject_ids(:)), run_count(:), ...
    max_abs_z(:), rms_value(:), peak_abs(:), auc_abs(:), n_nan(:), ...
    weights(:), is_outlier(:), ...
    'VariableNames', {'index', 'label', 'subject', 'n_runs', 'max_abs_robust_z', ...
    'rms', 'peak_abs', 'auc_abs', 'n_nan', 'weight', 'is_outlier'});
qc.curves = Y;
qc.time = stats.time(:);
qc.stats = stats;
qc.unit = unit;
qc.weighting = char(opts.Weighting);
qc.outlier_z_threshold = opts.OutlierZThreshold;
qc.skipped = stats.skipped;

if logical(opts.Plot)
    local_plot_qc(qc, opts);
end
end

function labels = local_unique_run_labels(subject_ids)
labels = cell(size(subject_ids));
seen = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:numel(subject_ids)
    sid = char(subject_ids{i});
    if isKey(seen, sid)
        seen(sid) = seen(sid) + 1;
    else
        seen(sid) = 1;
    end
    labels{i} = sprintf('%s_run%02d', sid, seen(sid));
end
end

function [max_abs_z, weights, is_outlier] = local_curve_weights(Y, weighting, threshold)
weighting = lower(strtrim(char(weighting)));
n = size(Y, 1);
center = local_nanmedian(Y, 1);
scale = 1.4826 .* local_nanmedian(abs(Y - repmat(center, n, 1)), 1);
scale_bad = scale == 0 | isnan(scale);
fallback_scale = local_nanstd(Y, 1);
scale(scale_bad) = fallback_scale(scale_bad);
scale(scale == 0 | isnan(scale)) = NaN;
Z = abs((Y - repmat(center, n, 1)) ./ repmat(scale, n, 1));
max_abs_z = zeros(n, 1);
for i = 1:n
    zi = Z(i, :);
    zi = zi(~isnan(zi));
    if ~isempty(zi)
        max_abs_z(i) = max(zi);
    end
end

is_outlier = max_abs_z > threshold;
weights = ones(n, 1);
switch weighting
    case 'none'
        return
    case {'exclude', 'omit'}
        weights(is_outlier) = 0;
    case {'huber', 'downweight'}
        weights = min(1, threshold ./ max(max_abs_z, eps));
    case {'bisquare', 'tukey'}
        u = max_abs_z ./ threshold;
        weights = (1 - u .^ 2) .^ 2;
        weights(u >= 1) = 0;
    otherwise
        error('Unknown Weighting: %s. Use none, exclude, huber, or bisquare.', weighting);
end
end

function local_plot_qc(qc, opts)
figure;
subplot(2, 1, 1);
imagesc(qc.time, 1:size(qc.curves, 1), qc.curves);
colorbar;
set(gca, 'YTick', 1:size(qc.curves, 1), 'YTickLabel', cellstr(qc.table.label), 'TickLabelInterpreter', 'none');
xlabel('Time / lag');
ylabel(char(opts.Unit));
title(sprintf('Curve QC heatmap: model=%s, condition=%s', char(opts.Model), char(opts.Condition)), ...
    'Interpreter', 'none');

subplot(2, 1, 2);
hold on;
plot(qc.time, qc.curves', 'Color', [0.65 0.65 0.65]);
weighted_mean = local_weighted_mean(qc.curves, qc.table.weight);
plot(qc.time, weighted_mean, 'k-', 'LineWidth', 2.5);
xlabel('Time / lag');
ylabel('Curve value');
title(sprintf('black=weighted mean, weighting=%s, effective n=%0.3g', ...
    qc.weighting, sum(qc.table.weight)), 'Interpreter', 'none');
end

function m = local_weighted_mean(Y, weights)
m = nan(1, size(Y, 2));
weights = double(weights(:));
for j = 1:size(Y, 2)
    y = Y(:, j);
    valid = ~isnan(y) & weights > 0;
    if any(valid)
        w = weights(valid);
        m(j) = sum(w .* y(valid)) ./ sum(w);
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

function m = local_max_omitnan(X, dim)
if dim ~= 2
    error('local_max_omitnan supports dim 2.');
end
m = nan(size(X, 1), 1);
for i = 1:size(X, 1)
    y = X(i, :);
    y = y(~isnan(y));
    if ~isempty(y), m(i) = max(y); end
end
end

function s = local_sum_omitnan(X, dim)
if dim ~= 2
    error('local_sum_omitnan supports dim 2.');
end
s = nan(size(X, 1), 1);
for i = 1:size(X, 1)
    y = X(i, :);
    y = y(~isnan(y));
    if ~isempty(y), s(i) = sum(y); end
end
end

function m = local_nanmedian(X, dim)
if nargin < 2, dim = 1; end
if dim == 1
    m = nan(1, size(X, 2));
    for j = 1:size(X, 2)
        y = X(:, j);
        y = y(~isnan(y));
        if ~isempty(y), m(j) = median(y); end
    end
else
    m = nan(size(X, 1), 1);
    for i = 1:size(X, 1)
        y = X(i, :);
        y = y(~isnan(y));
        if ~isempty(y), m(i) = median(y); end
    end
end
end

function s = local_nanstd(X, dim)
mu = local_nanmean(X, dim);
if dim == 1
    centered = X - repmat(mu, size(X, 1), 1);
    n = sum(~isnan(X), 1);
    centered(isnan(centered)) = 0;
    s = sqrt(sum(centered .^ 2, 1) ./ max(n - 1, 1));
elseif dim == 2
    centered = X - repmat(mu, 1, size(X, 2));
    n = sum(~isnan(X), 2);
    centered(isnan(centered)) = 0;
    s = sqrt(sum(centered .^ 2, 2) ./ max(n - 1, 1));
else
    error('local_nanstd supports dim 1 or 2.');
end
s(n < 2) = NaN;
end

function m = local_nanmean(X, dim)
valid = ~isnan(X);
den = sum(valid, dim);
X(~valid) = 0;
m = sum(X, dim) ./ den;
m(den == 0) = NaN;
end
