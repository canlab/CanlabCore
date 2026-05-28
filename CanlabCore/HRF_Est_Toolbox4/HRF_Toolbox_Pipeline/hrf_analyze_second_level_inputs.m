function analysis = hrf_analyze_second_level_inputs(second_level_inputs, varargin)
%HRF_ANALYZE_SECOND_LEVEL_INPUTS Analyze HRF map-score outputs across subjects.
%
% analysis = hrf_analyze_second_level_inputs(second_level_inputs, ...)
%
% second_level_inputs can be the table returned by
% hrf_collect_wholebrain_outputs or the corresponding CSV filename. This
% function reads *_beta_map_scores.csv or *_t_map_scores.csv files, builds a
% long subject-level table, optionally computes a condition contrast, and
% runs one-sample tests across subjects for each lag and signature/map score.

p = inputParser;
p.addRequired('second_level_inputs', @(x) istable(x) || ischar(x) || isstring(x));
p.addParameter('Object', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('ConditionA', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConditionB', '', @(x) ischar(x) || isstring(x));
p.addParameter('LagSeconds', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('ScoreColumns', {}, @(x) iscell(x) || isstring(x));
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.addParameter('OutputSummaryCsv', '', @(x) ischar(x) || isstring(x));
p.addParameter('OutputLongCsv', '', @(x) ischar(x) || isstring(x));
p.parse(second_level_inputs, varargin{:});
opts = p.Results;

inputs = local_read_inputs(second_level_inputs);
score_file_var = local_score_file_var(char(opts.Object));
if ~any(strcmp(score_file_var, inputs.Properties.VariableNames))
    error('second_level_inputs is missing column %s.', score_file_var);
end

[long_table, skipped] = local_build_long_table(inputs, score_file_var, opts);
if isempty(long_table)
    error('No usable score rows found. Skipped %d input row(s).', numel(skipped));
end

subject_table = local_subject_average(long_table);
summary = local_summarize(subject_table, opts.Alpha);

analysis = struct();
analysis.long_table = long_table;
analysis.subject_table = subject_table;
analysis.summary = summary;
analysis.interpretation = local_interpret(summary, opts.Alpha);
analysis.skipped = skipped;
analysis.object = char(opts.Object);
analysis.conditionA = char(opts.ConditionA);
analysis.conditionB = char(opts.ConditionB);
analysis.alpha = opts.Alpha;

if ~isempty(opts.OutputSummaryCsv)
    writetable(summary, char(opts.OutputSummaryCsv));
end
if ~isempty(opts.OutputLongCsv)
    writetable(subject_table, char(opts.OutputLongCsv));
end
end

function inputs = local_read_inputs(second_level_inputs)
if istable(second_level_inputs)
    inputs = second_level_inputs;
else
    inputs = readtable(char(second_level_inputs), 'TextType', 'string');
end
end

function varname = local_score_file_var(object_name)
switch lower(object_name)
    case {'beta', 'b'}
        varname = 'beta_scores_file';
    case {'t', 'tmap', 'tmaps'}
        varname = 't_scores_file';
    otherwise
        error('Unknown Object: %s. Use ''beta'' or ''t''.', object_name);
end
end

function [T, skipped] = local_build_long_table(inputs, score_file_var, opts)
rows = {};
skipped = struct('index', {}, 'subject', {}, 'reason', {});
missing_policy = lower(char(opts.MissingPolicy));

for i = 1:height(inputs)
    subject = local_table_value(inputs, i, 'subject');
    score_file = local_table_value(inputs, i, score_file_var);
    if isempty(score_file) || exist(score_file, 'file') ~= 2
        skipped = local_skip(skipped, i, subject, sprintf('missing score file %s', score_file), missing_policy);
        continue
    end

    S = readtable(score_file, 'TextType', 'string');
    score_cols = local_score_columns(S, opts.ScoreColumns);
    if isempty(score_cols)
        skipped = local_skip(skipped, i, subject, 'no numeric score columns', missing_policy);
        continue
    end
    S = local_filter_lags(S, opts.LagSeconds);

    if ~isempty(opts.ConditionA) && ~isempty(opts.ConditionB)
        rows = [rows; local_contrast_rows(S, score_cols, subject, char(opts.ConditionA), char(opts.ConditionB))]; %#ok<AGROW>
    else
        rows = [rows; local_condition_rows(S, score_cols, subject, char(opts.ConditionA))]; %#ok<AGROW>
    end
end

var_names = {'subject', 'condition', 'lag_index', 'lag_seconds', 'score_name', 'value'};
if isempty(rows)
    T = cell2table(cell(0, numel(var_names)), 'VariableNames', var_names);
else
    T = cell2table(rows, 'VariableNames', var_names);
    T.value = local_numeric_column(T.value);
    T.lag_index = local_numeric_column(T.lag_index);
    T.lag_seconds = local_numeric_column(T.lag_seconds);
end
end

function rows = local_condition_rows(S, score_cols, subject, condition)
rows = {};
if ~isempty(condition)
    if ~any(strcmp('condition', S.Properties.VariableNames))
        error('Score table is missing condition column.');
    end
    spec = local_condition_spec(S, condition);
    rows = local_averaged_condition_rows(S, score_cols, subject, spec);
else
    for r = 1:height(S)
        for c = 1:numel(score_cols)
            rows(end + 1, :) = {subject, local_condition_name(S, r), local_lag_index(S, r), ...
                local_lag_seconds(S, r), score_cols{c}, S.(score_cols{c})(r)}; %#ok<AGROW>
        end
    end
end
end

function rows = local_contrast_rows(S, score_cols, subject, conditionA, conditionB)
if ~any(strcmp('condition', S.Properties.VariableNames))
    error('Score table is missing condition column.');
end

specA = local_condition_spec(S, conditionA);
specB = local_condition_spec(S, conditionB);
cond = cellstr(string(S.condition));
A = S(ismember(cond, specA.matched_conditions), :);
B = S(ismember(cond, specB.matched_conditions), :);
rows = {};
if isempty(A) || isempty(B)
    return
end

lags = intersect(local_to_numeric(A.lag_index), local_to_numeric(B.lag_index), 'stable');
for li = 1:numel(lags)
    a = A(local_to_numeric(A.lag_index) == lags(li), :);
    b = B(local_to_numeric(B.lag_index) == lags(li), :);
    if isempty(a) || isempty(b), continue; end
    for c = 1:numel(score_cols)
        rows(end + 1, :) = {subject, sprintf('%s_minus_%s', specA.display_label, specB.display_label), ...
            local_lag_index(a, 1), local_lag_seconds(a, 1), score_cols{c}, ...
            local_mean_omitnan(a.(score_cols{c})) - local_mean_omitnan(b.(score_cols{c}))}; %#ok<AGROW>
    end
end
end

function rows = local_averaged_condition_rows(S, score_cols, subject, spec)
rows = {};
cond = cellstr(string(S.condition));
S = S(ismember(cond, spec.matched_conditions), :);
if isempty(S), return; end
lags = unique(local_to_numeric(S.lag_index), 'stable');
for li = 1:numel(lags)
    one_lag = S(local_to_numeric(S.lag_index) == lags(li), :);
    for c = 1:numel(score_cols)
        rows(end + 1, :) = {subject, spec.display_label, local_lag_index(one_lag, 1), ...
            local_lag_seconds(one_lag, 1), score_cols{c}, local_mean_omitnan(one_lag.(score_cols{c}))}; %#ok<AGROW>
    end
end
end

function spec = local_condition_spec(S, condition)
available = unique(cellstr(string(S.condition)), 'stable');
specs = hrf_resolve_condition_patterns(available, condition, 'DefaultMode', 'first');
spec = specs(1);
end

function S = local_filter_lags(S, lag_seconds)
if isempty(lag_seconds) || ~any(strcmp('lag_seconds', S.Properties.VariableNames))
    return
end

keep = false(height(S), 1);
all_lags = S.lag_seconds;
for i = 1:numel(lag_seconds)
    [~, idx] = min(abs(all_lags - lag_seconds(i)));
    keep = keep | all_lags == all_lags(idx);
end
S = S(keep, :);
end

function cols = local_score_columns(S, requested)
if ~isempty(requested)
    cols = cellstr(string(requested));
    missing = setdiff(cols, S.Properties.VariableNames);
    if ~isempty(missing)
        error('Requested score columns not found: %s', strjoin(missing, ', '));
    end
    return
end

metadata_cols = {'volume_index', 'condition', 'condition_index', 'lag_index', ...
    'lag_seconds', 'image_label', 'N', 'dfe', 'TR', 'mode'};
cols = {};
for i = 1:numel(S.Properties.VariableNames)
    name = S.Properties.VariableNames{i};
    if ismember(name, metadata_cols) || local_is_uncertainty_column(name)
        continue
    end
    if isnumeric(S.(name))
        cols{end + 1} = name; %#ok<AGROW>
    end
end
end

function tf = local_is_uncertainty_column(name)
tf = endsWith(char(name), '_se');
end

function subject_table = local_subject_average(T)
keys = local_group_keys(T);
rows = {};
for i = 1:size(keys, 1)
    wh = strcmp(T.subject, keys{i, 1}) & strcmp(T.condition, keys{i, 2}) & ...
        T.lag_index == keys{i, 3} & strcmp(T.score_name, keys{i, 5});
    rows(end + 1, :) = {keys{i, 1}, keys{i, 2}, keys{i, 3}, keys{i, 4}, keys{i, 5}, local_mean_omitnan(T.value(wh))}; %#ok<AGROW>
end
subject_table = cell2table(rows, 'VariableNames', T.Properties.VariableNames);
subject_table.value = local_numeric_column(subject_table.value);
subject_table.lag_index = local_numeric_column(subject_table.lag_index);
subject_table.lag_seconds = local_numeric_column(subject_table.lag_seconds);
end

function summary = local_summarize(T, alpha)
keys = local_summary_keys(T);
rows = {};
for i = 1:size(keys, 1)
    wh = strcmp(T.condition, keys{i, 1}) & T.lag_index == keys{i, 2} & strcmp(T.score_name, keys{i, 4});
    y = T.value(wh);
    y = y(~isnan(y));
    n = numel(y);
    if n > 1
        [~, pval, ~, st] = ttest(y);
        tval = st.tstat;
        sem = std(y) ./ sqrt(n);
    else
        pval = NaN;
        tval = NaN;
        sem = NaN;
    end
    rows(end + 1, :) = {keys{i, 1}, keys{i, 2}, keys{i, 3}, keys{i, 4}, ...
        n, local_mean_omitnan(y), sem, tval, pval, pval < alpha}; %#ok<AGROW>
end

summary = cell2table(rows, 'VariableNames', {'condition', 'lag_index', 'lag_seconds', ...
    'score_name', 'n', 'mean', 'sem', 't_value', 'p_value', 'significant'});
summary.lag_index = local_numeric_column(summary.lag_index);
summary.lag_seconds = local_numeric_column(summary.lag_seconds);
summary.n = local_numeric_column(summary.n);
summary.mean = local_numeric_column(summary.mean);
summary.sem = local_numeric_column(summary.sem);
summary.t_value = local_numeric_column(summary.t_value);
summary.p_value = local_numeric_column(summary.p_value);
summary.significant = logical(local_numeric_column(summary.significant));
end

function keys = local_summary_keys(T)
key_table = table(cellstr(string(T.condition)), T.lag_index, T.lag_seconds, cellstr(string(T.score_name)), ...
    'VariableNames', {'condition', 'lag_index', 'lag_seconds', 'score_name'});
[~, ia] = unique(key_table, 'rows', 'stable');
K = key_table(ia, :);
keys = [K.condition, num2cell(K.lag_index), num2cell(K.lag_seconds), K.score_name];
end

function keys = local_group_keys(T)
if any(strcmp('subject', T.Properties.VariableNames))
    subject = cellstr(string(T.subject));
else
    subject = repmat({''}, height(T), 1);
end
condition = cellstr(string(T.condition));
score_name = cellstr(string(T.score_name));
key_table = table(subject(:), condition(:), T.lag_index(:), T.lag_seconds(:), score_name(:), ...
    'VariableNames', {'subject', 'condition', 'lag_index', 'lag_seconds', 'score_name'});
[~, ia] = unique(key_table, 'rows', 'stable');
K = key_table(ia, :);
keys = [K.subject, K.condition, num2cell(K.lag_index), num2cell(K.lag_seconds), K.score_name];
end

function m = local_mean_omitnan(y)
y = y(~isnan(y));
if isempty(y)
    m = NaN;
else
    m = mean(y);
end
end

function x = local_numeric_column(x)
if iscell(x)
    x = cell2mat(x);
end
end

function x = local_to_numeric(x)
if isnumeric(x)
    x = double(x);
else
    x = str2double(string(x));
end
end

function interpretation = local_interpret(summary, alpha)
S = summary(~isnan(summary.p_value) & summary.p_value < alpha, :);
if isempty(S)
    interpretation = summary([], :);
    return
end
[~, idx] = sortrows([S.p_value, -abs(S.t_value)], [1 2]);
S = S(idx, :);
interpretation = S(1:min(20, height(S)), :);
end

function value = local_table_value(T, row, varname)
if ~any(strcmp(varname, T.Properties.VariableNames))
    value = '';
    return
end
value = T.(varname)(row);
if iscell(value), value = value{1}; end
value = char(string(value));
end

function name = local_condition_name(S, row)
if any(strcmp('condition', S.Properties.VariableNames))
    name = char(string(S.condition(row)));
else
    name = '';
end
end

function idx = local_lag_index(S, row)
if any(strcmp('lag_index', S.Properties.VariableNames))
    idx = S.lag_index(row);
else
    idx = row;
end
end

function seconds = local_lag_seconds(S, row)
if any(strcmp('lag_seconds', S.Properties.VariableNames))
    seconds = S.lag_seconds(row);
else
    seconds = NaN;
end
end

function skipped = local_skip(skipped, idx, subject, reason, missing_policy)
if strcmp(missing_policy, 'error')
    error('Input row %d (%s): %s', idx, subject, reason);
elseif ~strcmp(missing_policy, 'warn') && ~strcmp(missing_policy, 'silent')
    error('Unknown MissingPolicy: %s. Use ''warn'', ''silent'', or ''error''.', missing_policy);
end
skipped(end + 1) = struct('index', idx, 'subject', subject, 'reason', reason);
if strcmp(missing_policy, 'warn')
    warning('hrf_analyze_second_level_inputs:SkippingInput', 'Skipping input row %d (%s): %s', idx, subject, reason);
end
end
