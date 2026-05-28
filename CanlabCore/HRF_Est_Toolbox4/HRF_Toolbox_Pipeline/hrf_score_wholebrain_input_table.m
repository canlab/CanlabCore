function [input_table, status] = hrf_score_wholebrain_input_table(second_level_inputs, varargin)
%HRF_SCORE_WHOLEBRAIN_INPUT_TABLE Backfill map-score CSVs for collected HRF maps.
%
% [input_table, status] = hrf_score_wholebrain_input_table(input_table, ...)
%
% Iterates the table returned by hrf_collect_wholebrain_outputs and applies
% signatures/image sets to rows whose *_<object>_map_scores.csv files are
% missing or invalid. Updates beta_scores_file / t_scores_file in the
% returned table.
%
% This is now a thin orchestrator over HRF_SCORE_ONE_PREFIX, which is also
% called from the per-task SLURM worker and from HRF_AUDIT_SLURM_OUTPUTS
% (in RepairMissing mode). The per-prefix metadata resolution, validity
% checks, and statistic_image loading all live in HRF_SCORE_ONE_PREFIX so
% that the three call sites stay in sync.
%
% See also: hrf_score_one_prefix, hrf_audit_slurm_outputs,
%           hrf_apply_maps_to_wholebrain.

p = inputParser;
p.addRequired('second_level_inputs', @(x) istable(x) || ischar(x) || isstring(x));
p.addParameter('SourceModel', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ScoreObjects', {'beta', 't'}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('SignatureSets', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ImageSets', {}, @(x) ischar(x) || iscell(x) || isstring(x) || isa(x, 'image_vector'));
p.addParameter('AtlasObj', [], @(x) isempty(x) || isa(x, 'atlas') || isa(x, 'image_vector'));
p.addParameter('AtlasName', '', @(x) ischar(x) || isstring(x));
p.addParameter('Regions', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.addParameter('PropagateSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Overwrite', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('OverwriteStale', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Append', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('RequireMetadata', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('OutputCsv', '', @(x) ischar(x) || isstring(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.addParameter('NoVerbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(second_level_inputs, varargin{:});
opts = p.Results;

input_table = local_read_inputs(second_level_inputs);
source_models = local_source_model_filter(opts.SourceModel);
missing_policy = lower(char(opts.MissingPolicy));

status = table('Size', [0 6], ...
    'VariableTypes', {'double', 'string', 'string', 'string', 'logical', 'string'}, ...
    'VariableNames', {'row', 'subject', 'model', 'object', 'wrote_file', 'message'});

for i = 1:height(input_table)
    model_name = local_table_value(input_table, i, 'model');
    if ~local_source_model_matches(model_name, source_models)
        continue
    end

    prefix = local_table_value(input_table, i, 'prefix');
    subject = local_table_value(input_table, i, 'subject');
    if isempty(prefix)
        status = local_add_status(status, i, subject, model_name, '', false, 'missing prefix');
        local_handle_missing(missing_policy, i, subject, 'missing prefix');
        continue
    end

    result_mat_file = local_table_value(input_table, i, 'result_mat_file');
    metadata_file = local_table_value(input_table, i, 'metadata_file');

    helper_status = hrf_score_one_prefix(prefix, ...
        'ModelName', model_name, ...
        'NumModels', 1, ...
        'ScoreObjects', opts.ScoreObjects, ...
        'SignatureSets', opts.SignatureSets, ...
        'ImageSets', opts.ImageSets, ...
        'AtlasObj', opts.AtlasObj, ...
        'AtlasName', opts.AtlasName, ...
        'Regions', opts.Regions, ...
        'SimilarityMetric', opts.SimilarityMetric, ...
        'PropagateSE', opts.PropagateSE, ...
        'MetadataFile', metadata_file, ...
        'ResultMatFile', result_mat_file, ...
        'Overwrite', opts.Overwrite, ...
        'OverwriteStale', opts.OverwriteStale, ...
        'Append', opts.Append, ...
        'RequireMetadata', opts.RequireMetadata, ...
        'NoVerbose', opts.NoVerbose, ...
        'WarningContext', sprintf('row=%d; subject=%s', i, subject));

    if ~isempty(helper_status.metadata_file)
        input_table = local_set_file(input_table, i, 'metadata_file', helper_status.metadata_file);
    end

    objects_requested = local_resolve_score_objects_arg(opts.ScoreObjects);

    % Existing CSVs that the helper verified as valid
    for k = 1:numel(helper_status.skipped_existing)
        file_path = helper_status.skipped_existing{k};
        object_name = local_object_from_score_file(file_path, prefix, model_name);
        input_table = local_set_score_file(input_table, i, object_name, file_path);
        status = local_add_status(status, i, subject, model_name, object_name, false, 'exists');
    end

    % Newly written CSVs
    for k = 1:numel(helper_status.wrote_files)
        file_path = helper_status.wrote_files{k};
        object_name = local_object_from_score_file(file_path, prefix, model_name);
        input_table = local_set_score_file(input_table, i, object_name, file_path);
        status = local_add_status(status, i, subject, model_name, object_name, true, file_path);
    end

    % CSVs left stale (Overwrite=false, OverwriteStale=false)
    for k = 1:numel(helper_status.skipped_stale)
        file_path = helper_status.skipped_stale{k};
        object_name = local_object_from_score_file(file_path, prefix, model_name);
        msg = sprintf('existing score file is stale: %s', file_path);
        status = local_add_status(status, i, subject, model_name, object_name, false, msg);
        local_handle_missing(missing_policy, i, subject, msg);
    end

    % Errors
    for k = 1:numel(helper_status.errors)
        err = helper_status.errors(k);
        status = local_add_status(status, i, subject, model_name, err.object, false, err.message);
        local_handle_missing(missing_policy, i, subject, err.message);
    end

    % Catch-all: report any object that was requested but produced no
    % helper output (e.g. because RequireMetadata=false but metadata was
    % missing, so the helper bailed before iterating objects).
    reported_objects = string(status.object(status.row == i));
    for k = 1:numel(objects_requested)
        if ~any(reported_objects == string(objects_requested{k}))
            status = local_add_status(status, i, subject, model_name, ...
                objects_requested{k}, false, 'no action taken (no metadata or core inputs)');
        end
    end
end

if ~isempty(opts.OutputCsv)
    writetable(input_table, char(opts.OutputCsv));
end
end

% =========================================================================
% Local helpers
% =========================================================================
function input_table = local_read_inputs(second_level_inputs)
if istable(second_level_inputs)
    input_table = second_level_inputs;
else
    input_table = readtable(char(second_level_inputs), 'TextType', 'string');
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

function tf = local_source_model_matches(model_name, requested)
if isempty(requested)
    tf = true;
else
    tf = ismember(lower(strtrim(char(model_name))), requested);
end
end

function objects = local_resolve_score_objects_arg(score_objects)
objects = lower(cellstr(string(score_objects)));
objects = cellfun(@strtrim, objects, 'UniformOutput', false);
if any(strcmp(objects, 'both')) || any(strcmp(objects, 'all'))
    objects = {'beta', 't'};
end
objects(strcmp(objects, 'b')) = {'beta'};
objects(ismember(objects, {'tmap', 'tmaps'})) = {'t'};
objects = unique(objects(~cellfun(@isempty, objects)), 'stable');
end

function object_name = local_object_from_score_file(file_path, prefix, model_name)
[~, base] = fileparts(char(file_path));
prefix_base = char(prefix);
[~, prefix_name] = fileparts(prefix_base);
suffix = base;
if startsWith(suffix, prefix_name)
    suffix = suffix(numel(prefix_name) + 1:end);
end
if ~isempty(model_name)
    model_tag = ['_' lower(char(model_name))];
    if startsWith(suffix, model_tag)
        suffix = suffix(numel(model_tag) + 1:end);
    end
end
suffix = regexprep(suffix, '_map_scores$', '');
suffix = regexprep(suffix, '^_', '');
object_name = suffix;
if isempty(object_name) || ~ismember(object_name, {'beta', 't'})
    if endsWith(base, '_beta_map_scores')
        object_name = 'beta';
    elseif endsWith(base, '_t_map_scores')
        object_name = 't';
    else
        object_name = '';
    end
end
end

function input_table = local_set_score_file(input_table, row, object_name, score_file)
varname = local_score_varname(object_name);
if isempty(varname), return; end
input_table = local_set_file(input_table, row, varname, score_file);
end

function input_table = local_set_file(input_table, row, varname, file_name)
if ~any(strcmp(varname, input_table.Properties.VariableNames))
    input_table.(varname) = strings(height(input_table), 1);
end
if iscell(input_table.(varname))
    input_table.(varname){row} = char(file_name);
else
    input_table.(varname)(row) = string(file_name);
end
end

function varname = local_score_varname(object_name)
switch lower(char(object_name))
    case 'beta'
        varname = 'beta_scores_file';
    case 't'
        varname = 't_scores_file';
    otherwise
        varname = '';
end
end

function status = local_add_status(status, row, subject, model_name, object_name, wrote_file, message)
new_row = table(row, string(subject), string(model_name), string(object_name), ...
    logical(wrote_file), string(message), ...
    'VariableNames', status.Properties.VariableNames);
status = [status; new_row];
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

function local_handle_missing(missing_policy, row, subject, reason)
switch missing_policy
    case 'error'
        error('hrf_score_wholebrain_input_table:Row', ...
            'Input row %d (%s): %s', row, subject, reason);
    case 'warn'
        warning('hrf_score_wholebrain_input_table:SkippingInput', ...
            'Skipping input row %d (%s): %s', row, subject, reason);
    case 'silent'
        return
    otherwise
        error('hrf_score_wholebrain_input_table:UnknownMissingPolicy', ...
            'Unknown MissingPolicy: %s. Use warn, silent, or error.', missing_policy);
end
end
