function scores = hrf_apply_maps_to_wholebrain(stats_input, varargin)
%HRF_APPLY_MAPS_TO_WHOLEBRAIN Apply signatures/imagesets to 4D HRF maps.
%
% scores = hrf_apply_maps_to_wholebrain(stats_input, ...)
%
% stats_input may be:
%   - output from hrf_fit_wholebrain_stats
%   - a statistic_image/fmri_data/image_vector object
%   - a 4D NIfTI filename
%
% The output is a table with one row per 4D map volume and one column per
% signature/map score. This is the lightweight route for summaries after
% writing whole-brain beta/T maps.

p = inputParser;
p.addRequired('stats_input');
p.addParameter('Object', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('SignatureSets', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('ImageSets', {}, @(x) ischar(x) || iscell(x) || isstring(x) || isa(x, 'image_vector'));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.addParameter('SEInput', [], @(x) isempty(x) || isstruct(x) || isa(x, 'image_vector') || ischar(x) || isstring(x));
p.addParameter('PropagateSE', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('SEScoreSuffix', '_se', @(x) ischar(x) || isstring(x));
p.addParameter('MetadataTable', table(), @(x) isempty(x) || istable(x));
p.addParameter('OutputCsv', '', @(x) ischar(x) || isstring(x));
p.addParameter('WarningContext', '', @(x) ischar(x) || isstring(x));
p.parse(stats_input, varargin{:});
opts = p.Results;

[obj, metadata_table] = local_get_object(stats_input, opts.Object, opts.MetadataTable);
n_images = size(obj.dat, 2);
scores = local_base_table(obj, metadata_table, n_images);
empty_insertions = struct('score', {}, 'n_inserted', {});
se_obj = local_get_se_object(stats_input, opts.SEInput, opts.Object, obj);
propagate_se = logical(opts.PropagateSE) && ~isempty(se_obj) && local_is_linear_metric(opts.SimilarityMetric);
if logical(opts.PropagateSE) && ~isempty(se_obj) && ~local_is_linear_metric(opts.SimilarityMetric)
    warning('hrf_apply_maps_to_wholebrain:CannotPropagateSEForMetric', ...
        'Pattern-score SE propagation is only implemented for dotproduct/dot_product metrics; skipping SE columns for metric %s.', ...
        char(opts.SimilarityMetric));
end

signature_sets = local_to_cell(opts.SignatureSets);
for s = 1:numel(signature_sets)
    sigset = signature_sets{s};
    [signature_obj, signature_names] = local_load_signature_set(sigset, obj);
    if isempty(signature_names)
        warning('No signatures returned for image_set %s.', sigset);
        continue
    end
    for i = 1:numel(signature_names)
        name = signature_names{i};
        this_map = get_wh_image(signature_obj, i);
        v = apply_mask(obj, this_map, 'pattern_expression', 'ignore_missing', char(opts.SimilarityMetric));
        varname = local_unique_varname(scores, local_varname({'sig', sigset, name}));
        [v, n_inserted] = local_match_length(v, n_images, obj);
        scores.(varname) = v;
        empty_insertions = local_record_empty_insertion(empty_insertions, varname, n_inserted);
        if propagate_se
            scores = local_add_propagated_se(scores, se_obj, this_map, n_images, varname, opts.SEScoreSuffix);
        end
    end
end

image_sets = local_to_cell(opts.ImageSets);
for s = 1:numel(image_sets)
    image_set = image_sets{s};
    if isa(image_set, 'image_vector')
        maps = image_set;
        set_name = 'imageset';
        map_names = local_map_names(maps);
    else
        set_name = char(image_set);
        [maps, map_names] = local_load_named_image_set(set_name, obj);
    end

    for i = 1:numel(map_names)
        this_map = get_wh_image(maps, i);
        v = apply_mask(obj, this_map, 'pattern_expression', 'ignore_missing', char(opts.SimilarityMetric));
        varname = local_unique_varname(scores, local_varname({'map', set_name, map_names{i}}));
        [v, n_inserted] = local_match_length(v, n_images, obj);
        scores.(varname) = v;
        empty_insertions = local_record_empty_insertion(empty_insertions, varname, n_inserted);
        if propagate_se
            scores = local_add_propagated_se(scores, se_obj, this_map, n_images, varname, opts.SEScoreSuffix);
        end
    end
end

local_warn_empty_insertions(empty_insertions, opts.WarningContext, opts.OutputCsv);

if ~isempty(opts.OutputCsv)
    writetable(scores, char(opts.OutputCsv));
end
end

function se_obj = local_get_se_object(stats_input, se_input, which_obj, ref_obj)
se_obj = [];
if ~ismember(lower(char(which_obj)), {'beta', 'b'})
    return
end

if ~isempty(se_input)
    se_obj = local_as_image_vector(se_input);
elseif isstruct(stats_input) && isfield(stats_input, 'b') && local_has_ste(stats_input.b)
    se_obj = stats_input.b;
    se_obj.dat = stats_input.b.ste;
    se_obj.type = 'HRF beta standard error';
elseif isa(stats_input, 'statistic_image') && local_has_ste(stats_input)
    se_obj = stats_input;
    se_obj.dat = stats_input.ste;
    se_obj.type = 'HRF beta standard error';
end

if isempty(se_obj)
    return
end
local_validate_se_object(se_obj, ref_obj);
end

function tf = local_has_ste(obj)
tf = isprop(obj, 'ste') && ~isempty(obj.ste);
end

function obj = local_as_image_vector(input_obj)
if isa(input_obj, 'image_vector')
    obj = input_obj;
elseif ischar(input_obj) || isstring(input_obj)
    obj = statistic_image(fmri_data(char(input_obj), 'noverbose'));
else
    error('Unsupported SEInput type.');
end
end

function local_validate_se_object(se_obj, ref_obj)
if size(se_obj.dat, 1) ~= size(ref_obj.dat, 1) || size(se_obj.dat, 2) ~= size(ref_obj.dat, 2)
    error('SEInput image data size (%d-by-%d) does not match scored image data size (%d-by-%d).', ...
        size(se_obj.dat, 1), size(se_obj.dat, 2), size(ref_obj.dat, 1), size(ref_obj.dat, 2));
end
end

function tf = local_is_linear_metric(metric)
metric = lower(strrep(strtrim(char(metric)), '_', ''));
tf = ismember(metric, {'dotproduct', 'dot'});
end

function scores = local_add_propagated_se(scores, se_obj, weight_map, n_images, score_varname, suffix)
se_values = local_pattern_score_se(se_obj, weight_map, n_images);
se_varname = local_unique_varname(scores, [char(score_varname) char(suffix)]);
scores.(se_varname) = se_values;
end

function score_se = local_pattern_score_se(se_obj, weight_map, n_images)
if size(weight_map.dat, 1) ~= size(se_obj.dat, 1)
    error('Weight map voxel count (%d) does not match SE image voxel count (%d).', ...
        size(weight_map.dat, 1), size(se_obj.dat, 1));
end

w = double(weight_map.dat(:, 1));
valid_w = isfinite(w);
w(~valid_w) = 0;
w2 = w .^ 2;
score_se = nan(n_images, 1);

chunk_size = 64;
for first_img = 1:chunk_size:n_images
    last_img = min(first_img + chunk_size - 1, n_images);
    idx = first_img:last_img;
    se_dat = double(se_obj.dat(:, idx));
    valid = isfinite(se_dat) & repmat(valid_w, 1, numel(idx));
    se_dat(~valid) = 0;
    score_se(idx) = sqrt(w2' * (se_dat .^ 2))';
    score_se(idx(sum(valid, 1) == 0)) = NaN;
end
end

function [signature_obj, signature_names] = local_load_signature_set(sigset, ref_obj)
persistent signature_cache
if isempty(signature_cache)
    signature_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

key = local_cache_key('sig', sigset, ref_obj);
if isKey(signature_cache, key)
    entry = signature_cache(key);
    signature_obj = entry.obj;
    signature_names = entry.names;
    return
end

if ~exist('load_image_set', 'file')
    error('load_image_set not found on path. Add CanlabCore dependencies first.');
end
[signature_obj, signature_names] = load_image_set(char(sigset), 'noverbose');
signature_names = local_clean_signature_names(signature_names);
signature_obj = resample_space(signature_obj, ref_obj);
signature_cache(key) = struct('obj', signature_obj, 'names', {signature_names});
end

function [maps, map_names] = local_load_named_image_set(set_name, ref_obj)
persistent image_set_cache
if isempty(image_set_cache)
    image_set_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

key = local_cache_key('map', set_name, ref_obj);
if isKey(image_set_cache, key)
    entry = image_set_cache(key);
    maps = entry.obj;
    map_names = entry.names;
    return
end

if ~exist('load_image_set', 'file')
    error('load_image_set not found on path. Add CanlabCore dependencies first.');
end
[maps, map_names] = load_image_set(char(set_name), 'noverbose');
map_names = cellstr(string(map_names));
maps = resample_space(maps, ref_obj);
image_set_cache(key) = struct('obj', maps, 'names', {map_names});
end

function key = local_cache_key(kind, name, ref_obj)
parts = {char(kind), char(name), class(ref_obj)};
try
    vi = ref_obj.volInfo;
    parts{end + 1} = mat2str(local_volinfo_value(vi, 'dim'), 8);
    parts{end + 1} = mat2str(local_volinfo_value(vi, 'mat'), 8);
    parts{end + 1} = mat2str(local_volinfo_value(vi, 'nvox'), 8);
    parts{end + 1} = mat2str(local_volinfo_value(vi, 'n_inmask'), 8);
catch
    parts{end + 1} = sprintf('nvoxels-%d', size(ref_obj.dat, 1));
end
key = strjoin(parts, '|');
end

function value = local_volinfo_value(vi, field_name)
if isstruct(vi) && isfield(vi, field_name)
    value = vi.(field_name);
elseif isobject(vi) && isprop(vi, field_name)
    value = vi.(field_name);
else
    value = [];
end
value = double(value(:)');
end

function names = local_clean_signature_names(names)
names = cellstr(string(names));
names = strrep(names, '-', '_');
names = strrep(names, ' ', '_');
names = strrep(names, '.', '');
names = strrep(names, '(', '');
names = strrep(names, ')', '');
names = strrep(names, '^', '_');
end

function [obj, metadata_table] = local_get_object(stats_input, which_obj, metadata_table)
if isstruct(stats_input) && isfield(stats_input, 'b') && isfield(stats_input, 't')
    switch lower(char(which_obj))
        case {'beta', 'b'}
            obj = stats_input.b;
        case {'t', 'tmap', 'tmaps'}
            obj = stats_input.t;
        otherwise
            error('Unknown Object: %s. Use ''beta'' or ''t''.', char(which_obj));
    end
    if isempty(metadata_table) && isfield(stats_input, 'metadata_table')
        metadata_table = stats_input.metadata_table;
    end
elseif isa(stats_input, 'image_vector')
    obj = stats_input;
elseif ischar(stats_input) || isstring(stats_input)
    obj = statistic_image(char(stats_input), 'type', 'generic');
else
    error('Unsupported stats_input type.');
end
end

function T = local_base_table(obj, metadata_table, n_images)
if ~isempty(metadata_table)
    if height(metadata_table) ~= n_images
        error('Metadata row count (%d) does not match number of 4D images (%d). Regenerate matching whole-brain maps/metadata before scoring.', ...
            height(metadata_table), n_images);
    end
    T = metadata_table;
elseif isa(obj, 'statistic_image') && ~isempty(obj.image_labels) && numel(obj.image_labels) == n_images
    T = table((1:n_images)', obj.image_labels(:), 'VariableNames', {'volume_index', 'image_label'});
else
    T = table((1:n_images)', 'VariableNames', {'volume_index'});
end
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

function names = local_map_names(maps)
if isprop(maps, 'metadata_table') && ~isempty(maps.metadata_table) && ...
        any(strcmp('target', maps.metadata_table.Properties.VariableNames))
    names = cellstr(string(maps.metadata_table.target));
elseif isprop(maps, 'image_labels') && ~isempty(maps.image_labels)
    names = cellstr(string(maps.image_labels));
elseif ~isempty(maps.image_names)
    names = cellstr(string(maps.image_names));
else
    names = arrayfun(@(i) sprintf('map_%03d', i), 1:size(maps.dat, 2), 'UniformOutput', false);
end
end

function [v, n_inserted] = local_match_length(v, n_images, obj)
v = v(:);
n_inserted = 0;
if numel(v) == n_images
    return
end

empty_images = local_empty_images(obj, n_images);
if numel(v) + sum(empty_images) == n_images
    full_v = nan(n_images, 1);
    full_v(~empty_images) = v;
    v = full_v;
    n_inserted = sum(empty_images);
    return
end

error('Map/signature output length (%d) does not match number of 4D volumes (%d), and the mismatch cannot be explained by all-zero/all-NaN maps.', ...
    numel(v), n_images);
end

function empty_images = local_empty_images(obj, n_images)
empty_images = false(n_images, 1);
if ~isprop(obj, 'dat') || isempty(obj.dat) || size(obj.dat, 2) ~= n_images
    return
end

empty_images = all(obj.dat == 0 | isnan(obj.dat), 1)';

if isprop(obj, 'removed_images') && ~isempty(obj.removed_images) && ...
        numel(obj.removed_images) == n_images
    empty_images = empty_images | obj.removed_images(:);
end
end

function empty_insertions = local_record_empty_insertion(empty_insertions, score_name, n_inserted)
if n_inserted == 0
    return
end
empty_insertions(end + 1) = struct('score', char(score_name), 'n_inserted', n_inserted);
end

function local_warn_empty_insertions(empty_insertions, warning_context, output_csv)
if isempty(empty_insertions)
    return
end

n_each = [empty_insertions.n_inserted];
unique_n = unique(n_each);
score_names = {empty_insertions.score};
n_show = min(numel(score_names), 8);
shown_scores = strjoin(score_names(1:n_show), ', ');
if numel(score_names) > n_show
    shown_scores = sprintf('%s, ...', shown_scores);
end

context = char(warning_context);
if isempty(context)
    context = 'context not provided';
end
if ~isempty(output_csv)
    context = sprintf('%s; output_csv=%s', context, char(output_csv));
end

warning('hrf_apply_maps_to_wholebrain:ReinsertedEmptyImages', ...
    ['%s: Reinserted NaN scores for all-zero/all-NaN 4D volume(s) to preserve metadata alignment. ' ...
    '%d score column(s) affected; inserted counts={%s}; first affected columns: %s'], ...
    context, numel(empty_insertions), strjoin(cellstr(string(unique_n)), ', '), shown_scores);
end

function name = local_varname(parts)
parts = cellfun(@char, parts, 'UniformOutput', false);
name = matlab.lang.makeValidName(strjoin(parts, '_'));
end

function name = local_unique_varname(T, base)
name = base;
if istable(T)
    existing = T.Properties.VariableNames;
else
    existing = {};
end
if ~ismember(name, existing)
    return
end

k = 2;
while ismember(sprintf('%s_%d', base, k), existing)
    k = k + 1;
end
name = sprintf('%s_%d', base, k);
end
