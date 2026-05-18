function [input_table, status] = hrf_score_wholebrain_input_table(second_level_inputs, varargin)
%HRF_SCORE_WHOLEBRAIN_INPUT_TABLE Backfill map-score CSVs for collected HRF maps.
%
% [input_table, status] = hrf_score_wholebrain_input_table(input_table, ...)
%
% This helper takes the table returned by hrf_collect_wholebrain_outputs and
% applies signatures/image sets to rows whose *_map_scores.csv files are
% missing. It updates beta_scores_file/t_scores_file in the returned table.

p = inputParser;
p.addRequired('second_level_inputs', @(x) istable(x) || ischar(x) || isstring(x));
p.addParameter('SourceModel', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ScoreObjects', {'beta', 't'}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('SignatureSets', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ImageSets', {}, @(x) ischar(x) || iscell(x) || isstring(x) || isa(x, 'image_vector'));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.addParameter('Overwrite', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('OutputCsv', '', @(x) ischar(x) || isstring(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.addParameter('NoVerbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(second_level_inputs, varargin{:});
opts = p.Results;

input_table = local_read_inputs(second_level_inputs);
objects = local_score_objects(opts.ScoreObjects);
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

    [metadata_table, metadata_file] = local_metadata_for_row(input_table, i, prefix, model_name);
    if ~isempty(metadata_table) && ~isempty(metadata_file)
        input_table = local_set_file(input_table, i, 'metadata_file', metadata_file);
    end

    for j = 1:numel(objects)
        object_name = objects{j};
        score_file = local_score_file(input_table, i, object_name, prefix);
        if exist(score_file, 'file') == 2 && ~logical(opts.Overwrite)
            input_table = local_set_score_file(input_table, i, object_name, score_file);
            status = local_add_status(status, i, subject, model_name, object_name, false, 'exists');
            continue
        end

        try
            wholebrain = hrf_load_wholebrain_stats(prefix, 'NoVerbose', logical(opts.NoVerbose));
            hrf_apply_maps_to_wholebrain(wholebrain, ...
                'Object', object_name, ...
                'SignatureSets', opts.SignatureSets, ...
                'ImageSets', opts.ImageSets, ...
                'SimilarityMetric', opts.SimilarityMetric, ...
                'MetadataTable', metadata_table, ...
                'OutputCsv', score_file);
            input_table = local_set_score_file(input_table, i, object_name, score_file);
            status = local_add_status(status, i, subject, model_name, object_name, true, score_file);
        catch err
            status = local_add_status(status, i, subject, model_name, object_name, false, err.message);
            local_handle_missing(missing_policy, i, subject, err.message);
        end
    end
end

if ~isempty(opts.OutputCsv)
    writetable(input_table, char(opts.OutputCsv));
end
end

function input_table = local_read_inputs(second_level_inputs)
if istable(second_level_inputs)
    input_table = second_level_inputs;
else
    input_table = readtable(char(second_level_inputs), 'TextType', 'string');
end
end

function objects = local_score_objects(score_objects)
objects = lower(cellstr(string(score_objects)));
objects = cellfun(@strtrim, objects, 'UniformOutput', false);
if any(strcmp(objects, 'both')) || any(strcmp(objects, 'all'))
    objects = {'beta', 't'};
end
objects = unique(objects(~cellfun(@isempty, objects)), 'stable');
valid = {'beta', 'b', 't', 'tmap', 'tmaps'};
bad = setdiff(objects, valid);
if ~isempty(bad)
    error('Unknown ScoreObjects: %s. Use beta, t, or both.', strjoin(bad, ', '));
end
objects(strcmp(objects, 'b')) = {'beta'};
objects(ismember(objects, {'tmap', 'tmaps'})) = {'t'};
objects = unique(objects, 'stable');
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

function score_file = local_score_file(input_table, row, object_name, prefix)
varname = local_score_varname(object_name);
score_file = local_table_value(input_table, row, varname);
if ~isempty(score_file)
    return
end
score_file = [prefix '_' object_name '_map_scores.csv'];
end

function input_table = local_set_score_file(input_table, row, object_name, score_file)
varname = local_score_varname(object_name);
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
        error('Unknown object: %s. Use beta or t.', char(object_name));
end
end

function [metadata_table, metadata_file] = local_metadata_for_row(input_table, row, prefix, model_name)
metadata_file = local_table_value(input_table, row, 'metadata_file');
if isempty(metadata_file)
    metadata_file = [prefix '_metadata.csv'];
end
if exist(metadata_file, 'file') == 2
    metadata_table = readtable(metadata_file, 'TextType', 'string');
    return
end

metadata_table = local_metadata_from_result_mat(input_table, row, model_name);
if ~isempty(metadata_table)
    try
        writetable(metadata_table, metadata_file);
    catch
        metadata_file = '';
    end
end
end

function metadata_table = local_metadata_from_result_mat(input_table, row, model_name)
metadata_table = table();
mat_file = local_table_value(input_table, row, 'result_mat_file');
if isempty(mat_file) || exist(mat_file, 'file') ~= 2
    return
end

try
    S = load(mat_file, 'results');
catch
    return
end
if ~isfield(S, 'results')
    return
end
R = S.results;
model_field = matlab.lang.makeValidName(lower(char(model_name)));

if isfield(R, 'wholebrain_metadata_by_model') && isfield(R.wholebrain_metadata_by_model, model_field)
    metadata_table = R.wholebrain_metadata_by_model.(model_field);
elseif isfield(R, 'wholebrain_by_model') && isfield(R.wholebrain_by_model, model_field) && ...
        isfield(R.wholebrain_by_model.(model_field), 'metadata_table')
    metadata_table = R.wholebrain_by_model.(model_field).metadata_table;
elseif isfield(R, 'wholebrain_metadata_table')
    metadata_table = R.wholebrain_metadata_table;
elseif isfield(R, 'wholebrain') && isfield(R.wholebrain, 'metadata_table')
    metadata_table = R.wholebrain.metadata_table;
end
end

function status = local_add_status(status, row, subject, model_name, object_name, wrote_file, message)
new_row = table(row, string(subject), string(model_name), string(object_name), logical(wrote_file), string(message), ...
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
        error('Input row %d (%s): %s', row, subject, reason);
    case 'warn'
        warning('hrf_score_wholebrain_input_table:SkippingInput', ...
            'Skipping input row %d (%s): %s', row, subject, reason);
    case 'silent'
        return
    otherwise
        error('Unknown MissingPolicy: %s. Use warn, silent, or error.', missing_policy);
end
end
