function out = hrf_average_condition_trials(tc, E, cond_name, TR, window_seconds, varargin)
%HRF_AVERAGE_CONDITION_TRIALS Epoch and average trials for one condition.
% out = hrf_average_condition_trials(tc, E, 'pain', 0.8, 20, ...)
%
% Optional name/value
%   'BaselineSeconds' : seconds before onset used for baseline correction (default 0)

p = inputParser;
p.addRequired('tc', @(x) isnumeric(x) && isvector(x));
p.addRequired('E', @istable);
p.addRequired('cond_name', @(x) ischar(x) || isstring(x) || iscell(x));
p.addRequired('TR', @(x) isscalar(x) && x > 0);
p.addRequired('window_seconds', @(x) isscalar(x) && x > 0);
p.addParameter('BaselineSeconds', 0, @(x) isscalar(x) && x >= 0);
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

out = struct();
out.condition = condition_label;
out.matched_conditions = matched_conditions;
out.time = time;
out.trials = trials;
out.mean = mean(trials, 1)';
out.sem = std(trials, 0, 1)' ./ sqrt(size(trials, 1));
out.n_trials = size(trials, 1);
out.onset_idx = onset_idx;
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
