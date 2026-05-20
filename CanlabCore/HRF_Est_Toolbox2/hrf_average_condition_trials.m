function out = hrf_average_condition_trials(tc, E, cond_name, TR, window_seconds, varargin)
%HRF_AVERAGE_CONDITION_TRIALS Epoch and average trials for one condition.
% out = hrf_average_condition_trials(tc, E, 'pain', 0.8, 20, ...)
%
% Optional name/value
%   'BaselineSeconds' : seconds before onset used for baseline correction (default 0)
%   'OutlierPolicy'   : none, exclude, huber, or bisquare trial weighting
%   'OutlierZThreshold' : robust max-|z| threshold for trial weighting

p = inputParser;
p.addRequired('tc', @(x) isnumeric(x) && isvector(x));
p.addRequired('E', @istable);
p.addRequired('cond_name', @(x) ischar(x) || isstring(x) || iscell(x));
p.addRequired('TR', @(x) isscalar(x) && x > 0);
p.addRequired('window_seconds', @(x) isscalar(x) && x > 0);
p.addParameter('BaselineSeconds', 0, @(x) isscalar(x) && x >= 0);
p.addParameter('OutlierPolicy', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('OutlierZThreshold', 4, @(x) isscalar(x) && x > 0);
p.parse(tc, E, cond_name, TR, window_seconds, varargin{:});
opts = p.Results;

tc = tc(:);
available_conditions = unique(cellstr(string(E.trial_type)), 'stable');
condition_specs = hrf_resolve_condition_patterns(available_conditions, cond_name, 'DefaultMode', 'all');
matched_conditions = unique([condition_specs.matched_conditions], 'stable');
condition_label = local_condition_label(condition_specs);

if ~ismember('onset', E.Properties.VariableNames) || ~ismember('trial_type', E.Properties.VariableNames)
    error('E must contain onset and trial_type columns.');
end

n_tp = numel(tc);
win_tp = round(window_seconds / TR) + 1;
base_tp = round(opts.BaselineSeconds / TR);
time = (0:win_tp-1)' * TR;

idx = ismember(cellstr(string(E.trial_type)), matched_conditions);
onsets = E.onset(idx);
onset_idx = round(onsets ./ TR) + 1;

valid = onset_idx >= 1 - base_tp & (onset_idx + win_tp - 1) <= n_tp;
onset_idx = onset_idx(valid);

if isempty(onset_idx)
    error('No valid trials found for condition "%s" in selected time window.', condition_label);
end

trials = nan(numel(onset_idx), win_tp);
for i = 1:numel(onset_idx)
    s = onset_idx(i);
    segment = tc(s:(s + win_tp - 1));

    if base_tp > 0
        bstart = max(1, s - base_tp);
        bseg = tc(bstart:(s - 1));
        segment = segment - mean(bseg);
    end

    trials(i, :) = segment(:)';
end

[trial_weights, trial_is_outlier, trial_max_abs_z] = local_trial_weights(trials, ...
    opts.OutlierPolicy, opts.OutlierZThreshold);

out = struct();
out.condition = condition_label;
out.matched_conditions = matched_conditions;
out.time = time;
out.trials = trials;
if strcmpi(char(opts.OutlierPolicy), 'none')
    out.mean = mean(trials, 1)';
    out.sem = std(trials, 0, 1)' ./ sqrt(size(trials, 1));
else
    out.mean = local_weighted_mean(trials, trial_weights)';
    out.sem = local_weighted_sem(trials, trial_weights)';
end
out.n_trials = size(trials, 1);
out.n_trials_used = sum(trial_weights > 0);
out.trial_weights = trial_weights;
out.trial_is_outlier = trial_is_outlier;
out.trial_max_abs_z = trial_max_abs_z;
out.onset_idx = onset_idx;
end

function [weights, is_outlier, max_abs_z] = local_trial_weights(trials, policy, threshold)
policy = lower(strtrim(char(policy)));
n = size(trials, 1);
max_abs_z = zeros(n, 1);
is_outlier = false(n, 1);
weights = ones(n, 1);
if strcmp(policy, 'none') || n < 3
    return
end

center = local_nanmedian(trials, 1);
scale = 1.4826 .* local_nanmedian(abs(trials - repmat(center, n, 1)), 1);
scale_bad = scale == 0 | isnan(scale);
fallback_scale = local_nanstd(trials, 1);
scale(scale_bad) = fallback_scale(scale_bad);
scale(scale == 0 | isnan(scale)) = NaN;
Z = abs((trials - repmat(center, n, 1)) ./ repmat(scale, n, 1));
for i = 1:n
    z = Z(i, :);
    z = z(~isnan(z));
    if isempty(z)
        max_abs_z(i) = 0;
    else
        max_abs_z(i) = max(z);
    end
end
is_outlier = max_abs_z > threshold;

switch policy
    case {'exclude', 'omit'}
        weights(is_outlier) = 0;
    case {'huber', 'downweight'}
        weights = min(1, threshold ./ max(max_abs_z, eps));
    case {'bisquare', 'tukey'}
        u = max_abs_z ./ threshold;
        weights = (1 - u .^ 2) .^ 2;
        weights(u >= 1) = 0;
    otherwise
        error('Unknown OutlierPolicy: %s. Use none, exclude, huber, or bisquare.', policy);
end
end

function m = local_weighted_mean(X, w)
m = nan(1, size(X, 2));
for j = 1:size(X, 2)
    y = X(:, j);
    valid = ~isnan(y) & w > 0;
    if any(valid)
        ww = w(valid);
        m(j) = sum(ww .* y(valid)) ./ sum(ww);
    end
end
end

function se = local_weighted_sem(X, w)
se = nan(1, size(X, 2));
for j = 1:size(X, 2)
    y = X(:, j);
    valid = ~isnan(y) & w > 0;
    if sum(valid) < 2
        continue
    end
    yy = y(valid);
    ww = w(valid);
    mu = sum(ww .* yy) ./ sum(ww);
    n_eff = (sum(ww) .^ 2) ./ sum(ww .^ 2);
    var_w = sum(ww .* (yy - mu) .^ 2) ./ sum(ww);
    se(j) = sqrt(var_w ./ max(n_eff, eps));
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

function label = local_condition_label(condition_specs)
if isscalar(condition_specs)
    label = condition_specs.display_label;
else
    labels = cell(1, numel(condition_specs));
    for i = 1:numel(condition_specs)
        labels{i} = condition_specs(i).display_label;
    end
    label = strjoin(labels, ' + ');
end
end
