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
p.addParameter('MetadataTable', table(), @(x) isempty(x) || istable(x));
p.addParameter('OutputCsv', '', @(x) ischar(x) || isstring(x));
p.addParameter('WarningContext', '', @(x) ischar(x) || isstring(x));
p.parse(stats_input, varargin{:});
opts = p.Results;

[obj, metadata_table] = local_get_object(stats_input, opts.Object, opts.MetadataTable);
n_images = size(obj.dat, 2);
scores = local_base_table(obj, metadata_table, n_images);
empty_insertions = struct('score', {}, 'n_inserted', {});

signature_sets = local_to_cell(opts.SignatureSets);
for s = 1:numel(signature_sets)
    if ~exist('apply_all_signatures', 'file')
        error('apply_all_signatures not found on path. Add CanlabCore dependencies first.');
    end
    sigset = signature_sets{s};
    S = apply_all_signatures(obj, 'similarity_metric', char(opts.SimilarityMetric), 'image_set', sigset);
    if ~isfield(S, 'signaturenames') || isempty(S.signaturenames)
        warning('No signatures returned for image_set %s.', sigset);
        continue
    end
    for i = 1:numel(S.signaturenames)
        name = S.signaturenames{i};
        v = local_get_signal(S.(name));
        varname = local_unique_varname(scores, local_varname({'sig', sigset, name}));
        [v, n_inserted] = local_match_length(v, n_images, obj);
        scores.(varname) = v;
        empty_insertions = local_record_empty_insertion(empty_insertions, varname, n_inserted);
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
        [maps, map_names] = load_image_set(set_name, 'noverbose');
        map_names = cellstr(string(map_names));
    end

    for i = 1:numel(map_names)
        this_map = get_wh_image(maps, i);
        v = apply_mask(obj, this_map, 'pattern_expression', 'ignore_missing', char(opts.SimilarityMetric));
        varname = local_unique_varname(scores, local_varname({'map', set_name, map_names{i}}));
        [v, n_inserted] = local_match_length(v, n_images, obj);
        scores.(varname) = v;
        empty_insertions = local_record_empty_insertion(empty_insertions, varname, n_inserted);
    end
end

local_warn_empty_insertions(empty_insertions, opts.WarningContext, opts.OutputCsv);

if ~isempty(opts.OutputCsv)
    writetable(scores, char(opts.OutputCsv));
end
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

function v = local_get_signal(sig_struct)
if istable(sig_struct)
    vn = sig_struct.Properties.VariableNames;
    v = sig_struct.(vn{1});
elseif isstruct(sig_struct)
    f = fieldnames(sig_struct);
    v = sig_struct.(f{1});
else
    v = sig_struct;
end
if istable(v)
    vn = v.Properties.VariableNames;
    v = v.(vn{1});
end
v = v(:);
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
