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
% by default as a fallback to derive approximate score SE as
% abs(beta_score / t_score). If beta score CSVs contain propagated
% *_se columns from HRF beta SE maps, those are preferred.
% If multiple score columns are available, result.fits is the average score
% curve by default; individual signatures/maps remain in fits_by_signature.

p = inputParser;
p.addRequired('second_level_inputs', @(x) istable(x) || ischar(x) || isstring(x));
p.addParameter('Object', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('ScoreColumns', {}, @(x) iscell(x) || isstring(x));
p.addParameter('ModelName', 'mapscore', @(x) ischar(x) || isstring(x));
p.addParameter('SourceModel', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('AddAverageScore', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('AverageScoreName', 'mean_mapscore', @(x) ischar(x) || isstring(x));
p.addParameter('ApproxSEFromT', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.parse(second_level_inputs, varargin{:});
opts = p.Results;

inputs = local_read_inputs(second_level_inputs);
model_name = lower(strtrim(char(opts.ModelName)));
source_model_input = opts.SourceModel;
if isempty(local_source_model_filter(source_model_input)) && local_is_wholebrain_model_name(model_name)
    source_model_input = model_name;
    model_name = 'mapscore';
end
source_model_filter = local_source_model_filter(source_model_input);
score_file_var = local_score_file_var(char(opts.Object));
if ~any(strcmp(score_file_var, inputs.Properties.VariableNames))
    error('second_level_inputs is missing column %s.', score_file_var);
end

n = height(inputs);
results = cell(n, 1);
subject_ids = cell(n, 1);
run_labels = cell(n, 1);
success = false(n, 1);
errors = cell(n, 1);
skipped = struct('index', {}, 'subject', {}, 'reason', {});
all_score_names = {};
missing_policy = lower(char(opts.MissingPolicy));

for i = 1:n
    subject_ids{i} = local_table_value(inputs, i, 'subject');
    run_labels{i} = local_run_label(inputs, i);
    score_file = local_score_file_for_row(inputs, i, score_file_var);
    source_model = local_source_model(inputs, i, score_file);

    if ~local_source_model_matches(source_model, source_model_filter)
        reason = sprintf('source model %s did not match requested SourceModel', source_model);
        skipped = local_skip(skipped, i, subject_ids{i}, reason, missing_policy);
        errors{i} = reason;
        continue
    end

    if isempty(score_file) || exist(score_file, 'file') ~= 2
        reason = sprintf('missing score file %s', score_file);
        skipped = local_skip(skipped, i, subject_ids{i}, reason, missing_policy);
        errors{i} = reason;
        continue
    end

    try
        S = readtable(score_file, 'TextType', 'string');
        local_validate_score_table(S, inputs, i, score_file);
        Tscore = local_read_t_score_table(inputs, i, opts);
        if ~isempty(Tscore)
            local_validate_score_table(Tscore, inputs, i, local_score_file_for_row(inputs, i, 't_scores_file'));
        end
        score_cols = local_score_columns(S, opts.ScoreColumns);
        if isempty(score_cols)
            reason = 'no numeric score columns';
            skipped = local_skip(skipped, i, subject_ids{i}, reason, missing_policy);
            errors{i} = reason;
            continue
        end

        results{i} = local_score_table_to_result(S, Tscore, score_cols, model_name, ...
            char(opts.Object), score_file, source_model, logical(opts.AddAverageScore), ...
            char(opts.AverageScoreName));
        results{i}.subject_id = subject_ids{i};
        results{i}.run_label = run_labels{i};
        results{i}.input_row = i;
        results{i}.input_prefix = local_table_value(inputs, i, 'prefix');
        success(i) = true;
        if logical(opts.AddAverageScore)
            all_score_names = [all_score_names, {char(opts.AverageScoreName)}, score_cols(:)']; %#ok<AGROW>
        else
            all_score_names = [all_score_names, score_cols(:)']; %#ok<AGROW>
        end
    catch err
        skipped = local_skip(skipped, i, subject_ids{i}, err.message, missing_policy);
        errors{i} = err.message;
    end
end

study = struct();
study.results = results;
study.subject_ids = subject_ids;
study.run_labels = run_labels;
study.success = success;
study.errors = errors;
study.skipped = skipped;
study.score_names = unique(all_score_names, 'stable');
study.object = char(opts.Object);
study.model_name = model_name;
study.source_models = local_study_source_models(inputs);
study.source = 'second_level_inputs_map_scores';
study.second_level_inputs = inputs;
end

function run_label = local_run_label(inputs, row)
run_label = local_table_value(inputs, row, 'run_label');
if isempty(run_label)
    prefix = local_table_value(inputs, row, 'prefix');
    subject = local_table_value(inputs, row, 'subject');
    model_name = local_table_value(inputs, row, 'model');
    run_label = local_run_label_from_prefix(prefix, subject, model_name);
end
end

function run_label = local_run_label_from_prefix(prefix, subject, model_name)
run_label = '';
if isempty(prefix)
    return
end
[~, name] = fileparts(prefix);
run_label = regexprep(name, '_hrf.*$', '');
if ~isempty(subject)
    run_label = regexprep(run_label, ['^' regexptranslate('escape', subject) '_?'], '');
end
if ~isempty(model_name)
    run_label = regexprep(run_label, ['_' regexptranslate('escape', model_name) '$'], '');
end
if isempty(run_label)
    run_label = 'run-unknown';
end
end

function inputs = local_read_inputs(second_level_inputs)
if istable(second_level_inputs)
    inputs = second_level_inputs;
else
    inputs = readtable(char(second_level_inputs), 'TextType', 'string');
end
end

function models = local_source_model_filter(source_model)
if isempty(source_model)
    models = {};
else
    models = lower(cellstr(string(source_model)));
    models = cellfun(@strtrim, models, 'UniformOutput', false);
    models = models(~cellfun(@isempty, models));
end
end

function tf = local_is_wholebrain_model_name(model_name)
tf = ismember(lower(strtrim(char(model_name))), {'fir', 'sfir', 'canonical', 'spline'});
end

function tf = local_source_model_matches(source_model, requested)
if isempty(requested)
    tf = true;
else
    tf = ismember(lower(strtrim(source_model)), requested);
end
end

function source_model = local_source_model(inputs, row, score_file)
source_model = local_table_value(inputs, row, 'model');
if isempty(source_model)
    source_model = local_model_from_filename(score_file);
end
source_model = lower(strtrim(source_model));
end

function model_name = local_model_from_filename(file_name)
model_name = '';
if isempty(file_name)
    return
end
[~, name] = fileparts(file_name);
tok = regexp(name, '_([A-Za-z]+)_(beta|t)_map_scores$', 'tokens', 'once');
if ~isempty(tok) && ismember(lower(tok{1}), {'fir', 'sfir', 'canonical', 'spline'})
    model_name = lower(tok{1});
end
end

function models = local_study_source_models(inputs)
if any(strcmp('model', inputs.Properties.VariableNames))
    vals = cellstr(string(inputs.model));
    vals = vals(~cellfun(@(s) ismissing(string(s)) || strlength(string(s)) == 0, vals));
    models = unique(lower(vals), 'stable');
else
    models = {};
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

function score_file = local_score_file_for_row(inputs, row, score_file_var)
score_file = local_table_value(inputs, row, score_file_var);
if ~isempty(score_file)
    return
end
prefix = local_table_value(inputs, row, 'prefix');
if isempty(prefix)
    return
end
switch score_file_var
    case 'beta_scores_file'
        score_file = [prefix '_beta_map_scores.csv'];
    case 't_scores_file'
        score_file = [prefix '_t_map_scores.csv'];
end
end

function Tscore = local_read_t_score_table(inputs, row, opts)
Tscore = table();
if ~logical(opts.ApproxSEFromT) || ~strcmpi(char(opts.Object), 'beta')
    return
end
t_score_file = local_score_file_for_row(inputs, row, 't_scores_file');
if isempty(t_score_file) || exist(t_score_file, 'file') ~= 2
    return
end
Tscore = readtable(t_score_file, 'TextType', 'string');
end

function local_validate_score_table(S, inputs, row, score_file)
metadata_file = local_table_value(inputs, row, 'metadata_file');
if isempty(metadata_file) || exist(metadata_file, 'file') ~= 2
    return
end
M = readtable(metadata_file, 'TextType', 'string');
if height(S) ~= height(M)
    error('Score file %s has %d rows but metadata %s has %d rows. Regenerate map-score CSVs from the matching whole-brain maps.', ...
        score_file, height(S), metadata_file, height(M));
end
has_score_condition = any(strcmp('condition', S.Properties.VariableNames));
has_metadata_condition = any(strcmp('condition', M.Properties.VariableNames));
if has_metadata_condition && ~has_score_condition
    error('Score file %s is missing condition labels from metadata %s. Regenerate map-score CSVs.', ...
        score_file, metadata_file);
elseif has_score_condition && has_metadata_condition
    if ~isequal(string(S.condition), string(M.condition))
        error('Score file %s condition labels do not match metadata %s. Regenerate map-score CSVs.', ...
            score_file, metadata_file);
    end
end
has_score_lag = any(strcmp('lag_index', S.Properties.VariableNames));
has_metadata_lag = any(strcmp('lag_index', M.Properties.VariableNames));
if has_metadata_lag && ~has_score_lag
    error('Score file %s is missing lag indices from metadata %s. Regenerate map-score CSVs.', ...
        score_file, metadata_file);
elseif has_score_lag && has_metadata_lag
    if ~isequaln(local_to_numeric(S.lag_index), local_to_numeric(M.lag_index))
        error('Score file %s lag indices do not match metadata %s. Regenerate map-score CSVs.', ...
            score_file, metadata_file);
    end
end
end

function result = local_score_table_to_result(S, Tscore, score_cols, model_name, object_name, score_file, source_model, add_average, average_name)
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

if add_average
    signature_names = [{average_name}, score_cols(:)'];
    signature_fields = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(signature_names));
    average_field = signature_fields{1};
    score_fields = signature_fields(2:end);
else
    signature_names = score_cols(:)';
    signature_fields = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(signature_names));
    average_field = '';
    score_fields = signature_fields;
end
fits_by_signature = struct();
fit_data_by_score = cell(1, numel(score_cols));

for s = 1:numel(score_cols)
    H = nan(numel(lag_names), numel(condition_names));
    T = nan(numel(lag_names), numel(condition_names));
    propagated_SE = nan(numel(lag_names), numel(condition_names));
    values = S.(score_cols{s});
    t_values = local_t_values(Tscore, score_cols{s});
    se_values = local_score_se_values(S, score_cols{s});
    for c = 1:numel(condition_names)
        for l = 1:numel(lag_names)
            wh = strcmp(conditions, condition_names{c}) & lag_index == lag_names(l);
            H(l, c) = local_mean_omitnan(values(wh));
            if ~isempty(se_values)
                propagated_SE(l, c) = local_combine_independent_se(se_values(wh));
            end
            if ~isempty(t_values)
                T(l, c) = local_t_mean_for_bin(t_values, t_conditions, t_lag_index, ...
                    conditions, lag_index, condition_names{c}, lag_names(l), wh);
            end
        end
    end

    [SE, T, P, p_type, uncertainty_source] = local_score_uncertainty(H, T, propagated_SE, S);

    fit_data = local_score_fit_data(H, lag_names, lag_seconds, score_cols{s}, ...
        score_file, source_model, SE, T, P, p_type, uncertainty_source);
    fit_data_by_score{s} = fit_data;

    fit = struct();
    fit.(model_name) = fit_data;
    fits_by_signature.(score_fields{s}) = fit;
end

if add_average
    average_fit = struct();
    average_fit.(model_name) = local_average_score_fit(fit_data_by_score, average_name, score_cols, score_file, source_model);
    fits_by_signature.(average_field) = average_fit;
    selected_signature = average_name;
    selected_fit = average_fit;
else
    selected_signature = score_cols{1};
    selected_fit = fits_by_signature.(score_fields{1});
end

result = struct();
result.conditions = condition_names(:)';
result.fits_by_signature = fits_by_signature;
result.signature_meta = struct();
result.signature_meta.signal_source = 'wholebrain_map_score';
result.signature_meta.selected_signature = selected_signature;
result.signature_meta.selected_signatures = signature_names(:)';
result.signature_meta.selected_signature_fields = signature_fields(:)';
result.signature_meta.available_signatures = signature_names(:)';
result.signature_meta.n_signatures = numel(signature_names);
result.signature_meta.object = object_name;
result.signature_meta.source_file = score_file;
result.signature_meta.source_model = source_model;
result.fits = selected_fit;
result.settings = struct('signal_source', 'second_level_map_scores', ...
    'model_name', model_name, 'object', object_name, 'source_model', source_model);
end

function fit_data = local_score_fit_data(H, lag_names, lag_seconds, source_score, score_file, source_model, SE, T, P, p_type, uncertainty_source)
fit_data = struct();
fit_data.hrf = H;
fit_data.lag_index = lag_names(:);
fit_data.lag_seconds = lag_seconds(:);
fit_data.source_score = source_score;
fit_data.source_file = score_file;
fit_data.source_model = source_model;
fit_data.se = SE;
fit_data.t = T;
fit_data.p = P;
fit_data.p_type = p_type;
fit_data.uncertainty_source = uncertainty_source;
end

function fit_data = local_average_score_fit(fit_data_by_score, average_name, score_cols, score_file, source_model)
H = local_stack_field(fit_data_by_score, 'hrf');
SE = local_stack_field(fit_data_by_score, 'se');
T = local_stack_field(fit_data_by_score, 't');

fit_data = fit_data_by_score{1};
fit_data.hrf = local_mean_omitnan_stack(H);
fit_data.source_score = average_name;
fit_data.source_scores = score_cols(:)';
fit_data.source_file = score_file;
fit_data.source_model = source_model;

if ~isempty(SE)
    fit_data.se = local_combine_score_se(H, SE);
elseif size(H, 3) > 1
    fit_data.se = local_sem_omitnan_stack(H);
else
    fit_data.se = [];
end

if ~isempty(fit_data.se)
    fit_data.t = fit_data.hrf ./ fit_data.se;
    fit_data.t(~isfinite(fit_data.t)) = NaN;
else
    fit_data.t = local_mean_omitnan_stack(T);
end
fit_data.p = [];
fit_data.p_type = '';
fit_data.uncertainty_source = 'average across map-score columns';
end

function t_values = local_t_values(Tscore, score_col)
t_values = [];
if isempty(Tscore) || ~istable(Tscore) || ~any(strcmp(score_col, Tscore.Properties.VariableNames))
    return
end
t_values = Tscore.(score_col);
end

function se_values = local_score_se_values(S, score_col)
se_values = [];
se_col = [char(score_col) '_se'];
if any(strcmp(se_col, S.Properties.VariableNames))
    se_values = S.(se_col);
end
end

function se = local_combine_independent_se(vals)
vals = vals(:);
vals = vals(isfinite(vals));
if isempty(vals)
    se = NaN;
else
    se = sqrt(sum(vals .^ 2)) ./ numel(vals);
end
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

function [SE, T, P, p_type, uncertainty_source] = local_score_uncertainty(H, T, propagated_SE, S)
SE = [];
P = [];
p_type = '';
uncertainty_source = 'not stored in map-score CSV; use group-level stats across subjects or refit source time series';

if ~isempty(propagated_SE) && any(isfinite(propagated_SE(:)))
    SE = propagated_SE;
    T = H ./ SE;
    T(~isfinite(T)) = NaN;
    if any(strcmp('dfe', S.Properties.VariableNames))
        dfe = local_mean_omitnan(local_to_numeric(S.dfe));
        if ~isnan(dfe)
            P = 2 * (1 - tcdf(abs(T), dfe));
            P(P == 0) = eps;
            p_type = sprintf('Approximate two-tailed P-values from beta map-score / propagated beta-SE map-score, dfe = %.3f', dfe);
        end
    end
    uncertainty_source = 'propagated score SE from matching HRF beta SE maps, assuming diagonal voxel covariance';
    return
end

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

function X = local_stack_field(fit_data_by_score, field_name)
first_idx = [];
for i = 1:numel(fit_data_by_score)
    if isfield(fit_data_by_score{i}, field_name) && ~isempty(fit_data_by_score{i}.(field_name))
        first_idx = i;
        break
    end
end
if isempty(first_idx)
    X = [];
    return
end

first_x = fit_data_by_score{first_idx}.(field_name);
X = nan([size(first_x), numel(fit_data_by_score)]);
X(:, :, first_idx) = first_x;
for i = (first_idx + 1):numel(fit_data_by_score)
    if ~isfield(fit_data_by_score{i}, field_name) || isempty(fit_data_by_score{i}.(field_name))
        continue
    end
    this_x = fit_data_by_score{i}.(field_name);
    if isequal(size(this_x), size(first_x))
        X(:, :, i) = this_x;
    end
end
if ~isempty(X) && all(isnan(X(:)))
    X = [];
end
end

function M = local_mean_omitnan_stack(X)
if isempty(X)
    M = [];
    return
end
M = nan(size(X, 1), size(X, 2));
for r = 1:size(X, 1)
    for c = 1:size(X, 2)
        M(r, c) = local_mean_omitnan(squeeze(X(r, c, :)));
    end
end
end

function SE = local_sem_omitnan_stack(X)
if isempty(X)
    SE = [];
    return
end
SE = nan(size(X, 1), size(X, 2));
for r = 1:size(X, 1)
    for c = 1:size(X, 2)
        vals = squeeze(X(r, c, :));
        vals = vals(~isnan(vals));
        if numel(vals) > 1
            SE(r, c) = std(vals) ./ sqrt(numel(vals));
        end
    end
end
end

function SE = local_combine_score_se(H, score_SE)
SE = nan(size(H, 1), size(H, 2));
score_sem = local_sem_omitnan_stack(H);
for r = 1:size(score_SE, 1)
    for c = 1:size(score_SE, 2)
        vals = squeeze(score_SE(r, c, :));
        vals = vals(~isnan(vals));
        if isempty(vals)
            model_se = NaN;
        else
            model_se = sqrt(sum(vals .^ 2)) ./ numel(vals);
        end
        if isnan(model_se)
            SE(r, c) = score_sem(r, c);
        elseif isnan(score_sem(r, c))
            SE(r, c) = model_se;
        else
            SE(r, c) = sqrt(model_se .^ 2 + score_sem(r, c) .^ 2);
        end
    end
end
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
    if ismember(name, metadata_cols) || local_is_uncertainty_column(name)
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

function tf = local_is_uncertainty_column(name)
tf = endsWith(char(name), '_se');
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
