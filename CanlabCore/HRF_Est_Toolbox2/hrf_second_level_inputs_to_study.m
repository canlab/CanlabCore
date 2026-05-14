function study = hrf_second_level_inputs_to_study(second_level_inputs, varargin)
%HRF_SECOND_LEVEL_INPUTS_TO_STUDY Convert map-score CSVs to study structure.
%
% study = hrf_second_level_inputs_to_study(second_level_inputs, ...)
%
% second_level_inputs can be the table returned by
% hrf_collect_wholebrain_outputs or the corresponding CSV filename. This
% helper reads *_beta_map_scores.csv or *_t_map_scores.csv files and creates
% a study.results cell array compatible with hrf_time_unfolding_stats.
%
% The converted "fits" are map-score curves over condition x lag, not
% subject-level HRF model fits. The default model name is therefore
% "mapscore". For beta map-score CSVs, matching t map-score CSVs are used
% by default to derive approximate score SE as abs(beta_score / t_score).

p = inputParser;
p.addRequired('second_level_inputs', @(x) istable(x) || ischar(x) || isstring(x));
p.addParameter('Object', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('ScoreColumns', {}, @(x) iscell(x) || isstring(x));
p.addParameter('ModelName', 'mapscore', @(x) ischar(x) || isstring(x));
p.addParameter('ApproxSEFromT', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.parse(second_level_inputs, varargin{:});
opts = p.Results;

inputs = local_read_inputs(second_level_inputs);
score_file_var = local_score_file_var(char(opts.Object));
if ~any(strcmp(score_file_var, inputs.Properties.VariableNames))
    error('second_level_inputs is missing column %s.', score_file_var);
end

n = height(inputs);
results = cell(n, 1);
subject_ids = cell(n, 1);
success = false(n, 1);
errors = cell(n, 1);
skipped = struct('index', {}, 'subject', {}, 'reason', {});
all_score_names = {};
missing_policy = lower(char(opts.MissingPolicy));
model_name = char(opts.ModelName);

for i = 1:n
    subject_ids{i} = local_table_value(inputs, i, 'subject');
    score_file = local_table_value(inputs, i, score_file_var);

    if isempty(score_file) || exist(score_file, 'file') ~= 2
        reason = sprintf('missing score file %s', score_file);
        skipped = local_skip(skipped, i, subject_ids{i}, reason, missing_policy);
        errors{i} = reason;
        continue
    end

    try
        S = readtable(score_file, 'TextType', 'string');
        Tscore = local_read_t_score_table(inputs, i, opts);
        score_cols = local_score_columns(S, opts.ScoreColumns);
        if isempty(score_cols)
            reason = 'no numeric score columns';
            skipped = local_skip(skipped, i, subject_ids{i}, reason, missing_policy);
            errors{i} = reason;
            continue
        end

        results{i} = local_score_table_to_result(S, Tscore, score_cols, model_name, ...
            char(opts.Object), score_file);
        success(i) = true;
        all_score_names = [all_score_names, score_cols(:)']; %#ok<AGROW>
    catch err
        skipped = local_skip(skipped, i, subject_ids{i}, err.message, missing_policy);
        errors{i} = err.message;
    end
end

study = struct();
study.results = results;
study.subject_ids = subject_ids;
study.success = success;
study.errors = errors;
study.skipped = skipped;
study.score_names = unique(all_score_names, 'stable');
study.object = char(opts.Object);
study.model_name = model_name;
study.source = 'second_level_inputs_map_scores';
study.second_level_inputs = inputs;
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

function Tscore = local_read_t_score_table(inputs, row, opts)
Tscore = table();
if ~logical(opts.ApproxSEFromT) || ~strcmpi(char(opts.Object), 'beta')
    return
end
t_score_file = local_table_value(inputs, row, 't_scores_file');
if isempty(t_score_file) || exist(t_score_file, 'file') ~= 2
    return
end
Tscore = readtable(t_score_file, 'TextType', 'string');
end

function result = local_score_table_to_result(S, Tscore, score_cols, model_name, object_name, score_file)
conditions = local_condition_values(S);
condition_names = unique(conditions, 'stable');
lag_index = local_lag_index_values(S, conditions);
lag_names = unique(lag_index, 'stable');
lag_seconds = local_lag_seconds_by_index(S, lag_index, lag_names);
t_conditions = {};
t_lag_index = [];
if ~isempty(Tscore)
    t_conditions = local_condition_values(Tscore);
    t_lag_index = local_lag_index_values(Tscore, t_conditions);
end

score_fields = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(score_cols));
fits_by_signature = struct();

for s = 1:numel(score_cols)
    H = nan(numel(lag_names), numel(condition_names));
    T = nan(numel(lag_names), numel(condition_names));
    values = S.(score_cols{s});
    t_values = local_t_values(Tscore, score_cols{s});
    for c = 1:numel(condition_names)
        for l = 1:numel(lag_names)
            wh = strcmp(conditions, condition_names{c}) & lag_index == lag_names(l);
            H(l, c) = local_mean_omitnan(values(wh));
            if ~isempty(t_values)
                T(l, c) = local_t_mean_for_bin(t_values, t_conditions, t_lag_index, ...
                    conditions, lag_index, condition_names{c}, lag_names(l), wh);
            end
        end
    end

    [SE, P, p_type, uncertainty_source] = local_approx_uncertainty(H, T, S);

    fit = struct();
    fit.(model_name) = struct();
    fit.(model_name).hrf = H;
    fit.(model_name).lag_index = lag_names(:);
    fit.(model_name).lag_seconds = lag_seconds(:);
    fit.(model_name).source_score = score_cols{s};
    fit.(model_name).source_file = score_file;
    fit.(model_name).se = SE;
    fit.(model_name).t = T;
    fit.(model_name).p = P;
    fit.(model_name).p_type = p_type;
    fit.(model_name).uncertainty_source = uncertainty_source;
    fits_by_signature.(score_fields{s}) = fit;
end

result = struct();
result.conditions = condition_names(:)';
result.fits_by_signature = fits_by_signature;
result.signature_meta = struct();
result.signature_meta.signal_source = 'wholebrain_map_score';
result.signature_meta.selected_signature = score_cols{1};
result.signature_meta.selected_signatures = score_cols(:)';
result.signature_meta.selected_signature_fields = score_fields(:)';
result.signature_meta.available_signatures = score_cols(:)';
result.signature_meta.n_signatures = numel(score_cols);
result.signature_meta.object = object_name;
result.signature_meta.source_file = score_file;
result.fits = fits_by_signature.(score_fields{1});
result.settings = struct('signal_source', 'second_level_map_scores', ...
    'model_name', model_name, 'object', object_name);
end

function t_values = local_t_values(Tscore, score_col)
t_values = [];
if isempty(Tscore) || ~istable(Tscore) || ~any(strcmp(score_col, Tscore.Properties.VariableNames))
    return
end
t_values = Tscore.(score_col);
end

function t_mean = local_t_mean_for_bin(t_values, t_conditions, t_lag_index, conditions, lag_index, condition_name, lag_name, fallback_wh)
t_mean = NaN;
if numel(t_values) == numel(t_conditions) && numel(t_values) == numel(t_lag_index)
    wh = strcmp(t_conditions, condition_name) & t_lag_index == lag_name;
    t_mean = local_mean_omitnan(t_values(wh));
elseif numel(t_values) == numel(conditions) && numel(t_values) == numel(lag_index)
    t_mean = local_mean_omitnan(t_values(fallback_wh));
end
end

function [SE, P, p_type, uncertainty_source] = local_approx_uncertainty(H, T, S)
SE = [];
P = [];
p_type = '';
uncertainty_source = 'not stored in map-score CSV; use group-level stats across subjects or refit source time series';
if isempty(T) || all(isnan(T(:)))
    return
end

SE = abs(H ./ T);
SE(~isfinite(SE)) = NaN;
P = [];
if any(strcmp('dfe', S.Properties.VariableNames))
    dfe = local_mean_omitnan(local_to_numeric(S.dfe));
    if ~isnan(dfe)
        P = 2 * (1 - tcdf(abs(T), dfe));
        P(P == 0) = eps;
        p_type = sprintf('Approximate two-tailed P-values from beta-map score / t-map score ratio, dfe = %.3f', dfe);
    end
end
uncertainty_source = 'approximate score SE from matching beta and t map-score CSVs';
end

function conditions = local_condition_values(S)
if any(strcmp('condition', S.Properties.VariableNames))
    conditions_string = string(S.condition);
else
    conditions_string = repmat("all", height(S), 1);
end
conditions_string(ismissing(conditions_string) | strlength(conditions_string) == 0) = "all";
conditions = cellstr(conditions_string);
conditions = conditions(:);
end

function lag_index = local_lag_index_values(S, conditions)
if any(strcmp('lag_index', S.Properties.VariableNames))
    lag_index = local_to_numeric(S.lag_index(:));
else
    lag_index = zeros(height(S), 1);
    condition_names = unique(conditions, 'stable');
    for c = 1:numel(condition_names)
        wh = find(strcmp(conditions, condition_names{c}));
        lag_index(wh) = (1:numel(wh))';
    end
end
lag_index = double(lag_index);
end

function lag_seconds = local_lag_seconds_by_index(S, lag_index, lag_names)
lag_seconds = nan(numel(lag_names), 1);
if ~any(strcmp('lag_seconds', S.Properties.VariableNames))
    return
end

for i = 1:numel(lag_names)
    wh = lag_index == lag_names(i);
    lag_seconds(i) = local_mean_omitnan(local_to_numeric(S.lag_seconds(wh)));
end
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
    if ismember(name, metadata_cols)
        continue
    end
    if isnumeric(S.(name))
        cols{end + 1} = name; %#ok<AGROW>
    end
end
end

function m = local_mean_omitnan(y)
y = y(:);
y = y(~isnan(y));
if isempty(y)
    m = NaN;
else
    m = mean(y);
end
end

function value = local_table_value(T, row, varname)
if ~any(strcmp(varname, T.Properties.VariableNames))
    value = '';
    return
end
value = T.(varname)(row);
if iscell(value), value = value{1}; end
value = string(value);
if ismissing(value) || strlength(value) == 0
    value = '';
else
    value = char(value);
end
end

function x = local_to_numeric(x)
if isnumeric(x)
    x = double(x);
else
    x = str2double(string(x));
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
    warning('hrf_second_level_inputs_to_study:SkippingInput', ...
        'Skipping input row %d (%s): %s', idx, subject, reason);
end
end
