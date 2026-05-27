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
p.addParameter('PropagateSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Overwrite', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('OverwriteStale', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('RequireMetadata', true, @(x) islogical(x) || isnumeric(x));
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
    elseif logical(opts.RequireMetadata)
        msg = sprintf('missing metadata for %s; cannot write condition-labeled map scores', prefix);
        status = local_add_status(status, i, subject, model_name, '', false, msg);
        local_handle_missing(missing_policy, i, subject, msg);
        continue
    end

    for j = 1:numel(objects)
        object_name = objects{j};
        score_file = local_score_file(input_table, i, object_name, prefix);
        require_uncertainty = local_requires_uncertainty(prefix, object_name, opts);
        if exist(score_file, 'file') == 2 && ~logical(opts.Overwrite) && ...
                local_existing_score_is_valid(score_file, metadata_table, require_uncertainty, opts)
            input_table = local_set_score_file(input_table, i, object_name, score_file);
            status = local_add_status(status, i, subject, model_name, object_name, false, 'exists');
            continue
        elseif exist(score_file, 'file') == 2 && ~logical(opts.Overwrite) && ~logical(opts.OverwriteStale)
            msg = sprintf('existing score file is stale: %s', score_file);
            status = local_add_status(status, i, subject, model_name, object_name, false, msg);
            local_handle_missing(missing_policy, i, subject, msg);
            continue
        end

        try
            score_obj = local_load_score_object(prefix, object_name, metadata_table, logical(opts.NoVerbose));
            local_validate_score_object_metadata(score_obj, metadata_table, object_name, prefix);
            se_obj = local_uncertainty_object(prefix, object_name, metadata_table, logical(opts.NoVerbose), logical(opts.PropagateSE));
            hrf_apply_maps_to_wholebrain(score_obj, ...
                'Object', object_name, ...
                'SignatureSets', opts.SignatureSets, ...
                'ImageSets', opts.ImageSets, ...
                'SimilarityMetric', opts.SimilarityMetric, ...
                'SEInput', se_obj, ...
                'PropagateSE', opts.PropagateSE, ...
                'MetadataTable', metadata_table, ...
                'OutputCsv', score_file, ...
                'WarningContext', local_warning_context(i, subject, model_name, object_name, prefix));
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
if isempty(metadata_table)
    metadata_table = local_metadata_from_sibling(prefix, model_name);
end
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

function metadata_table = local_metadata_from_sibling(prefix, model_name)
metadata_table = table();
base_prefix = local_base_prefix(prefix, model_name);
candidate_files = { ...
    [base_prefix '_metadata.csv'], ...
    [base_prefix '_sfir_metadata.csv'], ...
    [base_prefix '_canonical_metadata.csv'], ...
    [base_prefix '_spline_metadata.csv']};
candidate_files = unique(candidate_files, 'stable');
this_file = [prefix '_metadata.csv'];

for i = 1:numel(candidate_files)
    fname = candidate_files{i};
    if strcmp(fname, this_file) || exist(fname, 'file') ~= 2
        continue
    end
    try
        T = readtable(fname, 'TextType', 'string');
    catch
        continue
    end
    if any(strcmp('condition', T.Properties.VariableNames)) && ...
            any(strcmp('lag_index', T.Properties.VariableNames))
        metadata_table = local_relabel_metadata(T, model_name);
        return
    end
end
end

function base_prefix = local_base_prefix(prefix, model_name)
base_prefix = char(prefix);
model_suffix = ['_' lower(char(model_name))];
if endsWith(lower(base_prefix), model_suffix)
    base_prefix = base_prefix(1:end - numel(model_suffix));
end
end

function T = local_relabel_metadata(T, model_name)
model_name = lower(char(model_name));
T.volume_index = (1:height(T))';
if any(strcmp('mode', T.Properties.VariableNames))
    T.mode = repmat(string(upper(model_name)), height(T), 1);
end
if any(strcmp('image_label', T.Properties.VariableNames)) && ...
        any(strcmp('condition', T.Properties.VariableNames)) && ...
        any(strcmp('lag_index', T.Properties.VariableNames))
    T.image_label = local_image_labels(T, model_name);
end
end

function labels = local_image_labels(T, model_name)
conditions = cellstr(string(T.condition));
lag_index = local_to_numeric(T.lag_index);
if any(strcmp('lag_seconds', T.Properties.VariableNames))
    lag_seconds = local_to_numeric(T.lag_seconds);
else
    lag_seconds = lag_index;
end
labels = cell(height(T), 1);
for i = 1:height(T)
    if ismember(model_name, {'canonical', 'spline'})
        labels{i} = sprintf('%s_%s_lag%03d_%0.3gs', ...
            local_safe_label(model_name), local_safe_label(conditions{i}), ...
            lag_index(i), lag_seconds(i));
    else
        labels{i} = sprintf('%s_lag%03d_%0.3gs', ...
            local_safe_label(conditions{i}), lag_index(i), lag_seconds(i));
    end
end
end

function s = local_safe_label(s)
s = matlab.lang.makeValidName(char(s));
end

function tf = local_existing_score_is_valid(score_file, metadata_table, require_uncertainty, opts)
tf = true;
try
    S = readtable(score_file, 'TextType', 'string');
catch
    tf = false;
    return
end
if logical(require_uncertainty) && ~local_score_table_has_uncertainty(S)
    tf = false;
    return
end
if ~local_has_requested_score_sets(S, opts)
    tf = false;
    return
end
if isempty(metadata_table)
    return
end
if height(S) ~= height(metadata_table)
    tf = false;
    return
end
has_score_condition = any(strcmp('condition', S.Properties.VariableNames));
has_metadata_condition = any(strcmp('condition', metadata_table.Properties.VariableNames));
if has_metadata_condition && ~has_score_condition
    tf = false;
    return
elseif has_score_condition && has_metadata_condition && ...
        ~isequal(string(S.condition), string(metadata_table.condition))
    tf = false;
    return
end
has_score_lag = any(strcmp('lag_index', S.Properties.VariableNames));
has_metadata_lag = any(strcmp('lag_index', metadata_table.Properties.VariableNames));
if has_metadata_lag && ~has_score_lag
    tf = false;
    return
elseif has_score_lag && has_metadata_lag && ...
        ~isequaln(local_to_numeric(S.lag_index), local_to_numeric(metadata_table.lag_index))
    tf = false;
end
end

function tf = local_has_requested_score_sets(S, opts)
tf = true;
sigsets = local_to_cell(opts.SignatureSets);
for i = 1:numel(sigsets)
    if ~local_has_numeric_prefix(S, local_var_prefix({'sig', sigsets{i}}))
        tf = false;
        return
    end
end

image_sets = local_to_cell(opts.ImageSets);
for i = 1:numel(image_sets)
    if isa(image_sets{i}, 'image_vector')
        set_name = 'imageset';
    else
        set_name = char(image_sets{i});
    end
    if ~local_has_numeric_prefix(S, local_var_prefix({'map', set_name}))
        tf = false;
        return
    end
end
end

function tf = local_has_numeric_prefix(S, prefix)
tf = false;
names = S.Properties.VariableNames;
for i = 1:numel(names)
    name = names{i};
    if startsWith(name, prefix) && isnumeric(S.(name)) && ~local_is_uncertainty_column(name)
        tf = true;
        return
    end
end
end

function prefix = local_var_prefix(parts)
prefix = matlab.lang.makeValidName(strjoin(cellfun(@char, parts, 'UniformOutput', false), '_'));
prefix = [prefix '_'];
end

function c = local_to_cell(x)
if isempty(x)
    c = {};
elseif isa(x, 'image_vector')
    c = {x};
elseif ischar(x) || isstring(x)
    c = cellstr(string(x));
else
    c = x;
end
end

function tf = local_requires_uncertainty(prefix, object_name, opts)
tf = strcmpi(char(object_name), 'beta') && logical(opts.PropagateSE) && ...
    local_is_linear_metric(opts.SimilarityMetric) && ...
    exist([prefix '_beta.nii'], 'file') == 2 && exist([prefix '_se.nii'], 'file') == 2;
end

function tf = local_is_linear_metric(metric)
metric = lower(strrep(strtrim(char(metric)), '_', ''));
tf = ismember(metric, {'dotproduct', 'dot'});
end

function tf = local_score_table_has_uncertainty(S)
score_cols = local_score_columns_for_validation(S);
tf = true;
for i = 1:numel(score_cols)
    if ~any(strcmp([score_cols{i} '_se'], S.Properties.VariableNames))
        tf = false;
        return
    end
end
end

function cols = local_score_columns_for_validation(S)
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

function x = local_to_numeric(x)
if isnumeric(x)
    x = double(x);
else
    x = str2double(string(x));
end
end

function obj = local_load_score_object(prefix, object_name, metadata_table, noverbose)
image_file = local_score_image_file(prefix, object_name);
if exist(image_file, 'file') ~= 2
    error('Missing %s image: %s', object_name, image_file);
end

load_args = {};
if noverbose
    load_args = {'noverbose'};
end
obj = statistic_image(fmri_data(image_file, load_args{:}));
switch lower(char(object_name))
    case 'beta'
        obj.type = 'FIR HRF beta';
    case 't'
        obj.type = 'T';
end

if ~isempty(metadata_table) && height(metadata_table) == size(obj.dat, 2)
    if any(strcmp('image_label', metadata_table.Properties.VariableNames))
        obj.image_labels = cellstr(string(metadata_table.image_label));
    end
    if any(strcmp('N', metadata_table.Properties.VariableNames))
        obj.N = metadata_table.N(1);
    end
    if any(strcmp('dfe', metadata_table.Properties.VariableNames))
        obj.dfe = metadata_table.dfe(1);
    end
end
if isempty(obj.sig)
    obj.sig = true(size(obj.dat));
end
end

function obj = local_uncertainty_object(prefix, object_name, metadata_table, noverbose, propagate_se)
obj = [];
if ~propagate_se || ~strcmpi(char(object_name), 'beta')
    return
end
image_file = [prefix '_se.nii'];
if exist(image_file, 'file') ~= 2
    return
end

load_args = {};
if noverbose
    load_args = {'noverbose'};
end
obj = statistic_image(fmri_data(image_file, load_args{:}));
obj.type = 'HRF beta standard error';
if ~isempty(metadata_table) && height(metadata_table) == size(obj.dat, 2)
    if any(strcmp('image_label', metadata_table.Properties.VariableNames))
        obj.image_labels = cellstr(string(metadata_table.image_label));
    end
    if any(strcmp('N', metadata_table.Properties.VariableNames))
        obj.N = metadata_table.N(1);
    end
    if any(strcmp('dfe', metadata_table.Properties.VariableNames))
        obj.dfe = metadata_table.dfe(1);
    end
end
if isempty(obj.sig)
    obj.sig = true(size(obj.dat));
end
end

function image_file = local_score_image_file(prefix, object_name)
switch lower(char(object_name))
    case 'beta'
        image_file = [prefix '_beta.nii'];
    case 't'
        image_file = [prefix '_t.nii'];
    otherwise
        error('Unknown object: %s. Use beta or t.', char(object_name));
end
end

function local_validate_score_object_metadata(obj, metadata_table, object_name, prefix)
if isempty(metadata_table)
    return
end
n_images = size(obj.dat, 2);
if height(metadata_table) ~= n_images
    error('Cannot score %s: metadata has %d rows but %s image has %d volumes. Regenerate the whole-brain maps and metadata together.', ...
        prefix, height(metadata_table), object_name, n_images);
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

function context = local_warning_context(row, subject, model_name, object_name, prefix)
context = sprintf('row=%d; subject=%s; model=%s; object=%s; prefix=%s', ...
    row, char(subject), char(model_name), char(object_name), char(prefix));
end
