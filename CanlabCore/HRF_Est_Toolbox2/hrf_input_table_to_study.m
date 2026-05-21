function study = hrf_input_table_to_study(second_level_inputs, varargin)
%HRF_INPUT_TABLE_TO_STUDY Rebuild study/results structs from collected files.
%
% study = hrf_input_table_to_study(second_level_inputs, ...)
%
% second_level_inputs can be the table returned by
% hrf_collect_wholebrain_outputs or the corresponding CSV filename. This
% function rebuilds a study.results cell array whose entries can contain:
%   .wholebrain          - beta/T statistic_image objects from NIfTI sidecars
%   .fits_by_signature  - optional map-score curves from *_map_scores.csv
%
% Use the .wholebrain field with hrf_apply_maps_to_wholebrain and
% hrf_animate_wholebrain_stats. Use the mapscore fits with
% plot_hrf_study_by_subject and hrf_time_unfolding_stats.

p = inputParser;
p.addRequired('second_level_inputs', @(x) istable(x) || ischar(x) || isstring(x));
p.addParameter('LoadWholeBrain', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('LoadResultMat', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('IncludeMapScores', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Object', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('ScoreColumns', {}, @(x) iscell(x) || isstring(x));
p.addParameter('ModelName', 'mapscore', @(x) ischar(x) || isstring(x));
p.addParameter('SourceModel', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('AddAverageScore', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('AverageScoreName', 'mean_mapscore', @(x) ischar(x) || isstring(x));
p.addParameter('ApproxSEFromT', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.addParameter('NoVerbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(second_level_inputs, varargin{:});
opts = p.Results;

inputs = local_read_inputs(second_level_inputs);
[model_name, source_model] = local_resolve_model_and_source(opts.ModelName, opts.SourceModel);
inputs = local_filter_inputs_by_source_model(inputs, source_model);
n = height(inputs);
results = cell(n, 1);
subject_ids = cell(n, 1);
run_labels = cell(n, 1);
success = false(n, 1);
wholebrain_success = false(n, 1);
resultmat_success = false(n, 1);
mapscore_success = false(n, 1);
errors = cell(n, 1);
skipped = struct('index', {}, 'subject', {}, 'reason', {});
missing_policy = lower(char(opts.MissingPolicy));

score_study = struct();
if logical(opts.IncludeMapScores)
    score_study = local_score_study(inputs, opts, missing_policy, model_name, source_model);
end

for i = 1:n
    subject_ids{i} = local_table_value(inputs, i, 'subject');
    run_labels{i} = local_run_label(inputs, i);
    r = struct();
    row_errors = {};

    if logical(opts.LoadResultMat)
        [mat_result, ok, reason] = local_load_result_mat_row(inputs, i);
        if ok
            r = mat_result;
            resultmat_success(i) = true;
        else
            skipped = local_skip(skipped, i, subject_ids{i}, reason, missing_policy);
            row_errors{end + 1} = reason; %#ok<AGROW>
        end
    end

    if logical(opts.LoadWholeBrain)
        [wholebrain, ok, reason] = local_load_wholebrain_row(inputs, i, opts);
        if ok
            r.wholebrain = wholebrain;
            r.wholebrain_paths = wholebrain.paths;
            r.conditions = wholebrain.conditions;
            r.settings = struct('signal_source', 'wholebrain_hrf_maps', ...
                'source_prefix', local_table_value(inputs, i, 'prefix'));
            wholebrain_success(i) = true;
        else
            skipped = local_skip(skipped, i, subject_ids{i}, reason, missing_policy);
            row_errors{end + 1} = reason; %#ok<AGROW>
        end
    end

    if logical(opts.IncludeMapScores)
        if isfield(score_study, 'results') && numel(score_study.results) >= i && ~isempty(score_study.results{i})
            r = local_merge_mapscore_result(r, score_study.results{i});
            mapscore_success(i) = true;
        else
            score_reason = local_score_error(score_study, i);
            if ~isempty(score_reason)
                row_errors{end + 1} = score_reason; %#ok<AGROW>
            end
        end
    end

    if ~isempty(fieldnames(r))
        r.subject_id = subject_ids{i};
        r.run_label = run_labels{i};
        r.input_row = i;
        r.input_prefix = local_table_value(inputs, i, 'prefix');
        results{i} = r;
        success(i) = resultmat_success(i) || wholebrain_success(i) || mapscore_success(i);
    end
    errors{i} = strjoin(row_errors, ' | ');
end

study = struct();
study.results = results;
study.subject_ids = subject_ids;
study.run_labels = run_labels;
study.success = success;
study.resultmat_success = resultmat_success;
study.wholebrain_success = wholebrain_success;
study.mapscore_success = mapscore_success;
study.errors = errors;
study.skipped = skipped;
study.object = char(opts.Object);
if isfield(score_study, 'model_name')
    study.model_name = score_study.model_name;
else
    study.model_name = model_name;
end
study.source = 'second_level_inputs';
study.second_level_inputs = inputs;

if isfield(score_study, 'score_names')
    study.score_names = score_study.score_names;
else
    study.score_names = {};
end
if isfield(score_study, 'source_models')
    study.source_models = score_study.source_models;
elseif any(strcmp('model', inputs.Properties.VariableNames))
    study.source_models = unique(lower(cellstr(string(inputs.model))), 'stable');
else
    study.source_models = {};
end
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

function [model_name, source_model] = local_resolve_model_and_source(model_name_in, source_model_in)
model_name = lower(strtrim(char(model_name_in)));
source_model = local_source_model_filter(source_model_in);
if isempty(source_model) && local_is_wholebrain_model_name(model_name)
    source_model = {model_name};
    model_name = 'mapscore';
end
end

function inputs = local_filter_inputs_by_source_model(inputs, source_model)
if isempty(source_model) || ~any(strcmp('model', inputs.Properties.VariableNames))
    return
end
row_models = lower(cellstr(string(inputs.model)));
keep = ismember(row_models, source_model);
inputs = inputs(keep, :);
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

function score_study = local_score_study(inputs, opts, missing_policy, model_name, source_model)
score_study = struct();
score_var = local_score_file_var(char(opts.Object));
if ~any(strcmp(score_var, inputs.Properties.VariableNames))
    return
end

try
    score_study = hrf_second_level_inputs_to_study(inputs, ...
        'Object', opts.Object, ...
        'ScoreColumns', opts.ScoreColumns, ...
        'ModelName', model_name, ...
        'SourceModel', source_model, ...
        'AddAverageScore', opts.AddAverageScore, ...
        'AverageScoreName', opts.AverageScoreName, ...
        'ApproxSEFromT', opts.ApproxSEFromT, ...
        'MissingPolicy', local_nested_missing_policy(missing_policy));
catch err
    if strcmp(missing_policy, 'error')
        rethrow(err);
    elseif strcmp(missing_policy, 'warn')
        warning('hrf_input_table_to_study:MapScoreLoadFailed', ...
            'Could not load map-score curves: %s', err.message);
    end
end
end

function reason = local_score_error(score_study, row)
reason = '';
if isfield(score_study, 'errors') && numel(score_study.errors) >= row
    reason = char(string(score_study.errors{row}));
end
if isempty(reason) && isfield(score_study, 'skipped')
    idx = find([score_study.skipped.index] == row, 1, 'first');
    if ~isempty(idx)
        reason = score_study.skipped(idx).reason;
    end
end
end

function nested_policy = local_nested_missing_policy(missing_policy)
if strcmp(missing_policy, 'error')
    nested_policy = 'error';
else
    nested_policy = 'silent';
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

function [result, ok, reason] = local_load_result_mat_row(inputs, row)
result = struct();
ok = false;
reason = '';

mat_file = local_table_value(inputs, row, 'result_mat_file');
if isempty(mat_file)
    prefix = local_table_value(inputs, row, 'prefix');
    if ~isempty(prefix)
        mat_file = [prefix '_results.mat'];
    end
end
if isempty(mat_file) || exist(mat_file, 'file') ~= 2
    reason = 'missing result_mat_file for run-level time series';
    return
end

try
    S = load(mat_file, 'results');
    if isfield(S, 'results')
        result = S.results;
    else
        tmp = load(mat_file);
        names = fieldnames(tmp);
        if isempty(names)
            reason = sprintf('empty result MAT file %s', mat_file);
            return
        end
        result = tmp.(names{1});
    end
    result.result_mat_file = mat_file;
    ok = true;
catch err
    reason = sprintf('could not load result MAT file %s: %s', mat_file, err.message);
end
end

function [wholebrain, ok, reason] = local_load_wholebrain_row(inputs, row, opts)
wholebrain = struct();
ok = false;
reason = '';

prefix = local_table_value(inputs, row, 'prefix');
if isempty(prefix)
    prefix = local_prefix_from_beta_file(local_table_value(inputs, row, 'beta_file'));
end
if isempty(prefix)
    reason = 'missing prefix/beta_file for whole-brain load';
    return
end

try
    wholebrain = hrf_load_wholebrain_stats(prefix, 'NoVerbose', logical(opts.NoVerbose));
    ok = true;
catch err
    reason = err.message;
end
end

function prefix = local_prefix_from_beta_file(beta_file)
prefix = '';
if isempty(beta_file)
    return
end
suffix = '_beta.nii';
if endsWith(beta_file, suffix)
    prefix = extractBefore(beta_file, strlength(beta_file) - strlength(suffix) + 1);
end
end

function r = local_merge_mapscore_result(r, score_result)
if ~isfield(r, 'fits') || isempty(r.fits)
    r.fits = score_result.fits;
end
if ~isfield(r, 'fits_by_signature') || isempty(r.fits_by_signature)
    r.fits_by_signature = score_result.fits_by_signature;
else
    r.fits_by_signature = local_merge_fit_structs(r.fits_by_signature, score_result.fits_by_signature);
end
if ~isfield(r, 'signature_meta') || isempty(r.signature_meta)
    r.signature_meta = score_result.signature_meta;
else
    r.mapscore_signature_meta = score_result.signature_meta;
end
r.mapscore_source = score_result.signature_meta.source_file;
if isfield(score_result.signature_meta, 'source_model')
    r.mapscore_source_model = score_result.signature_meta.source_model;
end
if ~isfield(r, 'conditions') || isempty(r.conditions)
    r.conditions = score_result.conditions;
end
if isfield(r, 'settings')
    r.settings.mapscore_model_name = score_result.settings.model_name;
    r.settings.mapscore_object = score_result.settings.object;
    if isfield(score_result.settings, 'source_model')
        r.settings.mapscore_source_model = score_result.settings.source_model;
    end
else
    r.settings = score_result.settings;
end
end

function out = local_merge_fit_structs(out, incoming)
fields = fieldnames(incoming);
for i = 1:numel(fields)
    f = fields{i};
    if ~isfield(out, f) || isempty(out.(f))
        out.(f) = incoming.(f);
        continue
    end
    model_names = fieldnames(incoming.(f));
    for j = 1:numel(model_names)
        out.(f).(model_names{j}) = incoming.(f).(model_names{j});
    end
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

function skipped = local_skip(skipped, idx, subject, reason, missing_policy)
if strcmp(missing_policy, 'error')
    error('Input row %d (%s): %s', idx, subject, reason);
elseif ~strcmp(missing_policy, 'warn') && ~strcmp(missing_policy, 'silent')
    error('Unknown MissingPolicy: %s. Use ''warn'', ''silent'', or ''error''.', missing_policy);
end

skipped(end + 1) = struct('index', idx, 'subject', subject, 'reason', reason);
if strcmp(missing_policy, 'warn')
    warning('hrf_input_table_to_study:SkippingInput', ...
        'Skipping input row %d (%s): %s', idx, subject, reason);
end
end
