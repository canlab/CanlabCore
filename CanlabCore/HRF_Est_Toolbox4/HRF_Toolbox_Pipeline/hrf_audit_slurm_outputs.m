function [audit, summary] = hrf_audit_slurm_outputs(root_or_manifest, varargin)
%HRF_AUDIT_SLURM_OUTPUTS Reconcile an HRF SLURM manifest with written files.
%
% [audit, summary] = hrf_audit_slurm_outputs(output_dir)
% [audit, summary] = hrf_audit_slurm_outputs(manifest_file, ...)
%
% This audits expected per-task/per-model whole-brain outputs, result MATs,
% map-score CSVs, and SLURM error logs.  Use this before collecting
% second-level inputs, because failed jobs may not create any *_beta.nii file
% and therefore will be invisible to hrf_collect_wholebrain_outputs().

p = inputParser;
p.addRequired('root_or_manifest', @(x) ischar(x) || isstring(x));
p.addParameter('ManifestFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConfigMat', '', @(x) ischar(x) || isstring(x));
p.addParameter('LogDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('Models', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ScoreObjects', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('RequireScores', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.addParameter('CheckNiftiVolumes', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('OutputCsv', '', @(x) ischar(x) || isstring(x));
p.addParameter('RepairMissing', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('SignatureSets', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ImageSets', {}, @(x) ischar(x) || iscell(x) || isstring(x) || isa(x, 'image_vector'));
p.addParameter('SimilarityMetric', '', @(x) ischar(x) || isstring(x));
p.addParameter('PropagateSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('UseParallel', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Verbose', false, @(x) islogical(x) || isnumeric(x));
p.parse(root_or_manifest, varargin{:});
opts = p.Results;

[root_dir, manifest_file] = local_resolve_manifest(root_or_manifest, opts.ManifestFile);
config_mat = local_default_file(opts.ConfigMat, fullfile(root_dir, 'hrf_study_config.mat'));
log_dir = local_default_file(opts.LogDir, fullfile(root_dir, 'logs'));

manifest = local_read_table(manifest_file);
config = local_load_config(config_mat);
models = local_expected_models(config, opts.Models);
score_objects = local_expected_score_objects(config, opts.ScoreObjects);
require_scores = local_require_scores(config, score_objects, opts.RequireScores);
logs = local_read_error_logs(log_dir);

rows = {};
for i = 1:height(manifest)
    task_index = local_table_number(manifest, i, 'index', i);
    subject = local_table_value(manifest, i, 'subject', sprintf('task-%03d', task_index));
    run_label = local_table_value(manifest, i, 'run_label', '');
    fmri_file = local_normalize_path(local_table_value(manifest, i, 'fmri_file', ''));
    events_file = local_normalize_path(local_table_value(manifest, i, 'events_file', ''));
    output_prefix = local_normalize_path(local_table_value(manifest, i, 'output_prefix', ''));
    output_mat = local_normalize_path(local_table_value(manifest, i, 'output_mat', ''));
    if isempty(output_mat) && ~isempty(output_prefix)
        output_mat = [output_prefix '_results.mat'];
    end

    [error_type, error_message, error_log_file] = local_task_error(logs, task_index);
    result_mat_exists = exist(output_mat, 'file') == 2;
    fmri_exists = exist(fmri_file, 'file') == 2;
    events_exists = exist(events_file, 'file') == 2;

    for m = 1:numel(models)
        model_name = models{m};
        model_prefix = local_model_prefix(output_prefix, model_name, numel(models));
        beta_file = [model_prefix '_beta.nii'];
        t_file = [model_prefix '_t.nii'];
        se_file = [model_prefix '_se.nii'];
        p_file = [model_prefix '_p.nii'];
        metadata_file = [model_prefix '_metadata.csv'];

        beta_exists = exist(beta_file, 'file') == 2;
        t_exists = exist(t_file, 'file') == 2;
        se_exists = exist(se_file, 'file') == 2;
        p_exists = exist(p_file, 'file') == 2;
        metadata_exists = exist(metadata_file, 'file') == 2;

        [metadata_rows, metadata_ok, metadata_error] = local_metadata_rows(metadata_file);
        [beta_vols, beta_ok, beta_error] = local_nifti_volumes(beta_file, opts.CheckNiftiVolumes);
        [t_vols, t_ok, t_error] = local_nifti_volumes(t_file, opts.CheckNiftiVolumes);
        [se_vols, se_ok, se_error] = local_nifti_volumes(se_file, opts.CheckNiftiVolumes);
        [p_vols, p_ok, p_error] = local_nifti_volumes(p_file, opts.CheckNiftiVolumes);

        if logical(opts.CheckNiftiVolumes)
            volume_match = local_rows_match(metadata_rows, [beta_vols t_vols se_vols p_vols]);
        else
            volume_match = metadata_ok;
        end
        nifti_ok = local_existing_niftis_ok([beta_exists t_exists se_exists p_exists], [beta_ok t_ok se_ok p_ok]);

        [score_files_exist, score_files_missing, beta_scores_file, t_scores_file] = ...
            local_score_status(output_prefix, model_name, numel(models), score_objects);
        score_complete = ~require_scores || score_files_exist;

        core_complete = result_mat_exists && beta_exists && t_exists && se_exists && p_exists && ...
            metadata_exists && metadata_ok && nifti_ok && volume_match;
        complete = core_complete && score_complete;

        failed_reason = local_failed_reason(error_type, result_mat_exists, fmri_exists, events_exists, ...
            beta_exists, t_exists, se_exists, p_exists, metadata_exists, metadata_ok, ...
            metadata_error, volume_match, require_scores, score_files_missing, ...
            beta_error, t_error, se_error, p_error);

        rows(end + 1, :) = {task_index, subject, run_label, model_name, output_prefix, model_prefix, ...
            fmri_file, events_file, output_mat, result_mat_exists, fmri_exists, events_exists, ...
            beta_file, t_file, se_file, p_file, metadata_file, beta_exists, t_exists, se_exists, p_exists, metadata_exists, ...
            metadata_rows, beta_vols, t_vols, se_vols, p_vols, beta_scores_file, t_scores_file, ...
            score_complete, score_files_missing, core_complete, complete, ...
            error_type, error_message, error_log_file, failed_reason}; %#ok<AGROW>
    end
end

var_names = {'task_index', 'subject', 'run_label', 'model', 'output_prefix', 'model_prefix', ...
    'fmri_file', 'events_file', 'result_mat_file', 'result_mat_exists', 'fmri_exists', 'events_exists', ...
    'beta_file', 't_file', 'se_file', 'p_file', 'metadata_file', ...
    'beta_exists', 't_exists', 'se_exists', 'p_exists', 'metadata_exists', ...
    'metadata_rows', 'beta_vols', 't_vols', 'se_vols', 'p_vols', ...
    'beta_scores_file', 't_scores_file', 'score_complete', 'score_files_missing', ...
    'core_complete', 'complete', 'error_type', 'error_message', 'error_log_file', 'failed_reason'};

if isempty(rows)
    audit = cell2table(cell(0, numel(var_names)), 'VariableNames', var_names);
else
    audit = cell2table(rows, 'VariableNames', var_names);
end

summary = local_summary(audit, manifest_file, config_mat, log_dir, models, require_scores);

if logical(opts.RepairMissing)
    [audit, repair_summary] = local_repair_missing_scores(audit, config, ...
        models, score_objects, opts);
    summary = local_summary(audit, manifest_file, config_mat, log_dir, models, require_scores);
    summary.repair = repair_summary;
end

if ~isempty(char(opts.OutputCsv))
    writetable(audit, char(opts.OutputCsv));
end
end

% =========================================================================
% Score-file repair (RepairMissing mode)
%   For every audit row where core_complete is true but score_complete is
%   not, call hrf_score_one_prefix to backfill the missing CSVs. After
%   scoring, re-evaluate score-file presence and update the audit table.
% =========================================================================
function [audit, repair_summary] = local_repair_missing_scores(audit, config, ...
    models, score_objects, opts)

repair_summary = struct( ...
    'n_rows_attempted', 0, ...
    'n_rows_repaired', 0, ...
    'n_rows_failed', 0, ...
    'n_files_written', 0, ...
    'errors', {{}}, ...
    'skipped_reason', '');

signature_sets = local_repair_signature_sets(opts.SignatureSets, config);
image_sets = local_repair_image_sets(opts.ImageSets, config);

if isempty(signature_sets) && isempty(image_sets)
    repair_summary.skipped_reason = ['No SignatureSets or ImageSets provided to RepairMissing ' ...
        'and none recorded in the SLURM config. Nothing to score.'];
    warning('hrf_audit_slurm_outputs:NothingToRepair', '%s', repair_summary.skipped_reason);
    return
end

similarity_metric = local_repair_similarity_metric(opts.SimilarityMetric, config);
propagate_se = logical(opts.PropagateSE);
n_models = max(numel(models), 1);

needs_repair = audit.core_complete & ~audit.score_complete;
row_indices = find(needs_repair);
repair_summary.n_rows_attempted = numel(row_indices);

if isempty(row_indices)
    return
end

if logical(opts.Verbose)
    fprintf('hrf_audit_slurm_outputs: RepairMissing -- attempting %d row(s).\n', numel(row_indices));
end

% Pre-compute the per-row argument bundles (this preserves parfor safety:
% parfor cannot slice a table by row, so we marshal each row's relevant
% scalars into a struct first).
arg_bundles = cell(numel(row_indices), 1);
for k = 1:numel(row_indices)
    r = row_indices(k);
    arg_bundles{k} = struct( ...
        'row', r, ...
        'output_prefix', local_char(audit.output_prefix(r)), ...
        'model_name', local_char(audit.model(r)), ...
        'result_mat_file', local_char(audit.result_mat_file(r)), ...
        'metadata_file', local_char(audit.metadata_file(r)), ...
        'subject', local_char(audit.subject(r)), ...
        'task_index', audit.task_index(r));
end

helper_statuses = cell(numel(arg_bundles), 1);
helper_errors = cell(numel(arg_bundles), 1);

if logical(opts.UseParallel)
    parfor k = 1:numel(arg_bundles)
        [helper_statuses{k}, helper_errors{k}] = local_repair_one_row( ...
            arg_bundles{k}, score_objects, signature_sets, image_sets, ...
            similarity_metric, propagate_se, n_models);
    end
else
    for k = 1:numel(arg_bundles)
        [helper_statuses{k}, helper_errors{k}] = local_repair_one_row( ...
            arg_bundles{k}, score_objects, signature_sets, image_sets, ...
            similarity_metric, propagate_se, n_models);
    end
end

for k = 1:numel(arg_bundles)
    r = arg_bundles{k}.row;
    s = helper_statuses{k};
    e = helper_errors{k};

    if ~isempty(s)
        repair_summary.n_files_written = repair_summary.n_files_written + numel(s.wrote_files);
        for ei = 1:numel(s.errors)
            repair_summary.errors{end + 1} = sprintf( ...
                'task %d (%s) object=%s: %s', ...
                arg_bundles{k}.task_index, arg_bundles{k}.subject, ...
                s.errors(ei).object, s.errors(ei).message);
        end

        % Update beta_scores_file / t_scores_file in the audit table from
        % whatever the helper produced or verified.
        for fi = 1:numel(s.wrote_files)
            audit = local_set_audit_scores_path(audit, r, s.wrote_files{fi});
        end
        for fi = 1:numel(s.skipped_existing)
            audit = local_set_audit_scores_path(audit, r, s.skipped_existing{fi});
        end
    end
    if ~isempty(e)
        repair_summary.errors{end + 1} = sprintf( ...
            'task %d (%s) helper error: %s', ...
            arg_bundles{k}.task_index, arg_bundles{k}.subject, e);
    end

    % Re-evaluate score-file presence for this row.
    [score_files_exist, score_files_missing, beta_scores_file, t_scores_file] = ...
        local_score_status(arg_bundles{k}.output_prefix, arg_bundles{k}.model_name, n_models, score_objects);
    audit.beta_scores_file{r} = beta_scores_file;
    audit.t_scores_file{r} = t_scores_file;
    audit.score_complete(r) = score_files_exist;
    audit.score_files_missing{r} = score_files_missing;
    audit.complete(r) = audit.core_complete(r) && audit.score_complete(r);
    audit.failed_reason{r} = local_recompute_failed_reason(audit, r);

    if audit.score_complete(r) && audit.core_complete(r)
        repair_summary.n_rows_repaired = repair_summary.n_rows_repaired + 1;
    else
        repair_summary.n_rows_failed = repair_summary.n_rows_failed + 1;
    end
end
end

function [helper_status, helper_error] = local_repair_one_row(bundle, score_objects, ...
    signature_sets, image_sets, similarity_metric, propagate_se, n_models)
helper_status = [];
helper_error = '';
try
    helper_status = hrf_score_one_prefix(bundle.output_prefix, ...
        'ModelName', bundle.model_name, ...
        'NumModels', n_models, ...
        'ScoreObjects', score_objects, ...
        'SignatureSets', signature_sets, ...
        'ImageSets', image_sets, ...
        'SimilarityMetric', similarity_metric, ...
        'PropagateSE', propagate_se, ...
        'MetadataFile', bundle.metadata_file, ...
        'ResultMatFile', bundle.result_mat_file, ...
        'Overwrite', false, ...
        'OverwriteStale', true, ...
        'WarningContext', sprintf('task=%d; subject=%s; repair=true', ...
            bundle.task_index, bundle.subject));
catch err
    helper_error = err.message;
end
end

function audit = local_set_audit_scores_path(audit, row, file_path)
[~, base] = fileparts(char(file_path));
if endsWith(base, '_beta_map_scores')
    audit.beta_scores_file{row} = char(file_path);
elseif endsWith(base, '_t_map_scores')
    audit.t_scores_file{row} = char(file_path);
end
end

function metric = local_repair_similarity_metric(metric_override, config)
metric = char(metric_override);
if ~isempty(metric), return; end
if isstruct(config) && isfield(config, 'similarity_metric') && ~isempty(config.similarity_metric)
    metric = char(config.similarity_metric);
else
    metric = 'dotproduct';
end
end

function sigsets = local_repair_signature_sets(sig_override, config)
if ~isempty(sig_override)
    sigsets = sig_override;
    return
end
sigsets = local_config_field(config, 'signature_sets', {});
end

function image_sets = local_repair_image_sets(image_override, config)
if ~isempty(image_override)
    image_sets = image_override;
    return
end
image_sets = local_config_field(config, 'image_sets', {});
end

function s = local_char(value)
if iscell(value)
    if isempty(value), s = ''; return; end
    value = value{1};
end
if isstring(value)
    if isempty(value) || ismissing(value) || strlength(value) == 0, s = ''; return; end
    s = char(value);
elseif ischar(value)
    s = value;
elseif isnumeric(value)
    s = num2str(value);
else
    try, s = char(string(value)); catch, s = ''; end
end
end

function reason = local_recompute_failed_reason(audit, row)
parts = {};
existing = char(audit.error_type{row});
if ~isempty(existing)
    parts{end + 1} = ['log error: ' existing];
end
if ~audit.result_mat_exists(row), parts{end + 1} = 'missing result MAT'; end
if ~audit.fmri_exists(row), parts{end + 1} = 'missing input fMRI'; end
if ~audit.events_exists(row), parts{end + 1} = 'missing events file'; end
missing = {};
if ~audit.beta_exists(row), missing{end + 1} = 'beta'; end
if ~audit.t_exists(row), missing{end + 1} = 't'; end
if ~audit.se_exists(row), missing{end + 1} = 'se'; end
if ~audit.p_exists(row), missing{end + 1} = 'p'; end
if ~audit.metadata_exists(row), missing{end + 1} = 'metadata'; end
if ~isempty(missing)
    parts{end + 1} = ['missing outputs: ' strjoin(missing, ',')];
end
if ~audit.score_complete(row) && ~isempty(char(audit.score_files_missing{row}))
    parts{end + 1} = ['missing score CSVs: ' char(audit.score_files_missing{row})];
end
reason = strjoin(parts, '; ');
end

function [root_dir, manifest_file] = local_resolve_manifest(root_or_manifest, manifest_override)
root_or_manifest = char(root_or_manifest);
if ~isempty(char(manifest_override))
    manifest_file = char(manifest_override);
    root_dir = root_or_manifest;
elseif exist(root_or_manifest, 'dir') == 7
    root_dir = root_or_manifest;
    manifest_file = fullfile(root_dir, 'hrf_study_manifest.csv');
else
    manifest_file = root_or_manifest;
    root_dir = fileparts(manifest_file);
    if isempty(root_dir), root_dir = pwd; end
end

manifest_file = local_normalize_path(manifest_file);
root_dir = local_normalize_path(root_dir);
if exist(manifest_file, 'file') ~= 2
    error('Manifest file not found: %s', manifest_file);
end
end

function path_out = local_default_file(user_path, default_path)
if isempty(char(user_path))
    path_out = local_normalize_path(default_path);
else
    path_out = local_normalize_path(char(user_path));
end
end

function T = local_read_table(filename)
try
    T = readtable(filename, 'TextType', 'string', 'Delimiter', ',', 'VariableNamingRule', 'preserve');
catch
    T = readtable(filename, 'TextType', 'string', 'Delimiter', ',');
end
end

function config = local_load_config(config_mat)
config = struct();
if exist(config_mat, 'file') ~= 2
    return
end
S = load(config_mat, 'config');
if isfield(S, 'config')
    config = S.config;
end
end

function models = local_expected_models(config, override)
if ~isempty(override)
    requested = local_to_cell(override);
    model_args = requested;
else
    pipeline_args = local_config_field(config, 'pipeline_args', {});
    requested = local_to_cell(local_arg_value(pipeline_args, 'WholeBrainMode'));
    model_args = local_to_cell(local_arg_value(pipeline_args, 'Models'));
    if isempty(requested), requested = {'auto'}; end
    if isempty(model_args)
        model_args = {'logit', 'fir', 'sfir', 'canonical', 'spline', 'nlgamma'};
    end
end

if isscalar(requested) && strcmpi(requested{1}, 'auto')
    requested = model_args;
end

supported = {'fir', 'sfir', 'canonical', 'spline'};
models = {};
for i = 1:numel(requested)
    name = lower(strtrim(char(requested{i})));
    if strcmp(name, 'auto')
        nested = local_expected_models(struct('pipeline_args', {{'Models', model_args, 'WholeBrainMode', 'auto'}}), {});
        models = [models nested]; %#ok<AGROW>
    elseif ismember(name, supported)
        models{end + 1} = name; %#ok<AGROW>
    end
end
models = unique(models, 'stable');
if isempty(models)
    models = {'fir'};
end
end

function score_objects = local_expected_score_objects(config, override)
if ~isempty(override)
    score_objects = local_to_cell(override);
else
    score_objects = local_to_cell(local_config_field(config, 'score_objects', {'beta'}));
end
score_objects = cellfun(@(x) lower(char(x)), score_objects, 'UniformOutput', false);
end

function tf = local_require_scores(config, score_objects, override)
if ~isempty(override)
    tf = logical(override);
    return
end
has_signature_sets = ~isempty(local_to_cell(local_config_field(config, 'signature_sets', {})));
has_image_sets = ~isempty(local_to_cell(local_config_field(config, 'image_sets', {})));
tf = ~isempty(score_objects) && (has_signature_sets || has_image_sets);
end

function value = local_config_field(config, field_name, default_value)
if isstruct(config) && isfield(config, field_name)
    value = config.(field_name);
else
    value = default_value;
end
end

function value = local_arg_value(args, name)
value = [];
if isempty(args) || ~iscell(args)
    return
end
names = args(1:2:end);
idx = find(strcmpi(string(names), string(name)), 1, 'last');
if ~isempty(idx)
    value = args{idx * 2};
end
end

function c = local_to_cell(x)
if isempty(x)
    c = {};
elseif iscell(x)
    c = x;
elseif isstring(x)
    c = cellstr(x);
else
    c = {x};
end
end

function value = local_table_value(T, row_idx, name, default_value)
col = find(strcmpi(T.Properties.VariableNames, name), 1);
if isempty(col)
    value = default_value;
    return
end
raw = T{row_idx, col};
value = local_cell_to_char(raw);
if isempty(value)
    value = default_value;
end
end

function value = local_table_number(T, row_idx, name, default_value)
txt = local_table_value(T, row_idx, name, '');
value = str2double(txt);
if isnan(value)
    raw_col = find(strcmpi(T.Properties.VariableNames, name), 1);
    if ~isempty(raw_col)
        raw = T{row_idx, raw_col};
        if isnumeric(raw) && isscalar(raw)
            value = raw;
            return
        end
    end
    value = default_value;
end
end

function out = local_cell_to_char(raw)
if iscell(raw)
    if isempty(raw) || isempty(raw{1})
        out = '';
        return
    end
    raw = raw{1};
end
if isa(raw, 'missing')
    out = '';
    return
end
try
    s = string(raw);
    if isempty(s) || ismissing(s(1)) || strlength(s(1)) == 0
        out = '';
    else
        out = char(s(1));
    end
catch
    out = '';
end
end

function path_out = local_normalize_path(path_in)
path_out = char(path_in);
if isunix
    path_out = strrep(path_out, '\', filesep);
end
end

function prefix = local_model_prefix(output_prefix, model_name, n_models)
if n_models == 1
    prefix = output_prefix;
else
    prefix = sprintf('%s_%s', output_prefix, lower(char(model_name)));
end
end

function [metadata_rows, ok, message] = local_metadata_rows(metadata_file)
metadata_rows = NaN;
ok = false;
message = '';
if exist(metadata_file, 'file') ~= 2
    return
end
try
    T = local_read_table(metadata_file);
    metadata_rows = height(T);
    ok = true;
catch err
    message = err.message;
end
end

function [n_vols, ok, message] = local_nifti_volumes(filename, do_check)
n_vols = NaN;
ok = false;
message = '';
if exist(filename, 'file') ~= 2
    return
end
if ~logical(do_check)
    ok = true;
    return
end
try
    info = niftiinfo(filename);
    if numel(info.ImageSize) >= 4
        n_vols = info.ImageSize(4);
    else
        n_vols = 1;
    end
    ok = true;
catch err
    message = err.message;
end
end

function tf = local_existing_niftis_ok(exists_flags, ok_flags)
tf = all(~logical(exists_flags) | logical(ok_flags));
end

function tf = local_rows_match(metadata_rows, image_vols)
if isnan(metadata_rows)
    tf = false;
    return
end
present_vols = image_vols(~isnan(image_vols));
tf = ~isempty(present_vols) && all(present_vols == metadata_rows);
end

function [all_exist, missing, beta_scores_file, t_scores_file] = local_score_status(output_prefix, model_name, n_models, score_objects)
all_exist = true;
missing_names = {};
beta_scores_file = '';
t_scores_file = '';
for i = 1:numel(score_objects)
    object_name = lower(char(score_objects{i}));
    score_file = local_score_file(output_prefix, model_name, n_models, object_name);
    if strcmp(object_name, 'beta'), beta_scores_file = score_file; end
    if strcmp(object_name, 't'), t_scores_file = score_file; end
    if exist(score_file, 'file') ~= 2
        all_exist = false;
        missing_names{end + 1} = object_name; %#ok<AGROW>
    end
end
missing = strjoin(missing_names, ',');
end

function score_file = local_score_file(output_prefix, model_name, n_models, object_name)
if n_models == 1
    score_file = sprintf('%s_%s_map_scores.csv', output_prefix, object_name);
else
    score_file = sprintf('%s_%s_%s_map_scores.csv', output_prefix, lower(char(model_name)), object_name);
end
end

function logs = local_read_error_logs(log_dir)
logs = struct('task_index', {}, 'file', {}, 'text', {});
if exist(log_dir, 'dir') ~= 7
    return
end
files = dir(fullfile(log_dir, '*.err'));
if ~isempty(files)
    [~, order] = sort([files.datenum], 'descend');
    files = files(order);
end
for i = 1:numel(files)
    tok = regexp(files(i).name, '_(\d+)\.err$', 'tokens', 'once');
    if isempty(tok)
        continue
    end
    task_index = str2double(tok{1});
    path_in = fullfile(files(i).folder, files(i).name);
    try
        txt = fileread(path_in);
    catch
        txt = '';
    end
    logs(end + 1) = struct('task_index', task_index, 'file', path_in, 'text', txt); %#ok<AGROW>
end
end

function [error_type, error_message, error_log_file] = local_task_error(logs, task_index)
error_type = '';
error_message = '';
error_log_file = '';
for i = 1:numel(logs)
    if logs(i).task_index ~= task_index
        continue
    end
    [candidate_type, candidate_message] = local_classify_error(logs(i).text);
    if ~isempty(candidate_type)
        error_type = candidate_type;
        error_message = candidate_message;
        error_log_file = logs(i).file;
        return
    end
end
end

function [error_type, error_message] = local_classify_error(txt)
error_type = '';
error_message = '';
if isempty(txt)
    return
end
if contains(txt, 'Image size does not match header description')
    error_type = 'nifti_image_size_mismatch';
elseif contains(txt, 'zeroinsert') && contains(txt, 'replace_empty')
    error_type = 'image_write_removed_images';
elseif contains(txt, 'Unable to write header')
    error_type = 'nifti_write_failed';
elseif contains(txt, 'Map/signature output length')
    error_type = 'mapscore_length_mismatch';
elseif contains(txt, 'missing metadata')
    error_type = 'missing_metadata';
elseif contains(txt, 'Error using')
    error_type = 'matlab_error';
end
if ~isempty(error_type)
    error_message = local_first_error_block(txt);
end
end

function message = local_first_error_block(txt)
lines = regexp(txt, '\r\n|\n|\r', 'split');
lines = lines(~contains(lines, 'Unable to load ApplicationService'));
start_idx = find(contains(lines, 'Error using'), 1);
if isempty(start_idx)
    patterns = {'Unable to perform assignment', 'Image size does not match header description', ...
        'Unable to write header', 'Out of memory', 'File not found or permission denied'};
    start_idx = [];
    for p = 1:numel(patterns)
        start_idx = find(contains(lines, patterns{p}), 1);
        if ~isempty(start_idx), break; end
    end
    if isempty(start_idx)
        message = strtrim(strjoin(lines(~cellfun('isempty', strtrim(lines))), ' '));
        return
    end
end
parts = {};
for i = start_idx:min(numel(lines), start_idx + 4)
    line = strtrim(lines{i});
    if isempty(line)
        continue
    end
    if i > start_idx && startsWith(line, 'Error in')
        break
    end
    parts{end + 1} = line; %#ok<AGROW>
end
message = strjoin(parts, ' ');
end

function reason = local_failed_reason(error_type, result_mat_exists, fmri_exists, events_exists, ...
    beta_exists, t_exists, se_exists, p_exists, metadata_exists, metadata_ok, ...
    metadata_error, volume_match, require_scores, score_files_missing, ...
    beta_error, t_error, se_error, p_error)
parts = {};
if ~isempty(error_type)
    parts{end + 1} = ['log error: ' error_type];
end
if ~result_mat_exists, parts{end + 1} = 'missing result MAT'; end
if ~fmri_exists, parts{end + 1} = 'missing input fMRI'; end
if ~events_exists, parts{end + 1} = 'missing events file'; end
missing = {};
if ~beta_exists, missing{end + 1} = 'beta'; end
if ~t_exists, missing{end + 1} = 't'; end
if ~se_exists, missing{end + 1} = 'se'; end
if ~p_exists, missing{end + 1} = 'p'; end
if ~metadata_exists, missing{end + 1} = 'metadata'; end
if ~isempty(missing)
    parts{end + 1} = ['missing outputs: ' strjoin(missing, ',')];
end
if metadata_exists && ~metadata_ok
    parts{end + 1} = ['cannot read metadata: ' metadata_error];
end
if metadata_ok && ~volume_match
    parts{end + 1} = 'metadata/image volume mismatch';
end
nifti_errors = local_join_nonempty({beta_error, t_error, se_error, p_error});
if ~isempty(nifti_errors)
    parts{end + 1} = ['cannot read output NIfTI header: ' nifti_errors];
end
if require_scores && ~isempty(score_files_missing)
    parts{end + 1} = ['missing score CSVs: ' score_files_missing];
end
reason = strjoin(parts, '; ');
end

function out = local_join_nonempty(values)
keep = false(size(values));
for i = 1:numel(values)
    keep(i) = ~isempty(values{i});
end
out = strjoin(values(keep), ' | ');
end

function summary = local_summary(audit, manifest_file, config_mat, log_dir, models, require_scores)
summary = struct();
summary.manifest_file = manifest_file;
summary.config_mat = config_mat;
summary.log_dir = log_dir;
summary.expected_models = models;
summary.require_scores = require_scores;
summary.n_manifest_tasks = numel(unique(audit.task_index));
summary.n_expected_model_outputs = height(audit);
summary.n_core_complete = sum(audit.core_complete);
summary.n_complete = sum(audit.complete);
summary.n_missing_or_failed = sum(~audit.complete);
summary.n_result_mats = numel(unique(audit.result_mat_file(audit.result_mat_exists)));

task_ids = unique(audit.task_index);
task_complete = false(size(task_ids));
for i = 1:numel(task_ids)
    task_complete(i) = all(audit.complete(audit.task_index == task_ids(i)));
end
summary.n_complete_tasks = sum(task_complete);
summary.n_missing_or_failed_tasks = sum(~task_complete);

error_types = audit.error_type(~cellfun('isempty', audit.error_type));
summary.error_types = unique(error_types);
end
