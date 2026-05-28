function [TC, meta] = hrf_extract_imageset_timeseries(fmri_nii, image_set, varargin)
%HRF_EXTRACT_IMAGESET_TIMESERIES Apply load_image_set maps to each fMRI volume.
% Supports Buckner networks, Margulies gradients, Hansen receptor maps, etc.

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x) || isa(x, 'fmri_data'));
p.addRequired('image_set', @(x) ischar(x) || isstring(x) || isa(x, 'image_vector'));
p.addParameter('MapNames', {}, @(x) iscell(x) || isstring(x));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.parse(fmri_nii, image_set, varargin{:});
opts = p.Results;

if isa(fmri_nii, 'fmri_data')
    dat = fmri_nii;
else
    dat = fmri_data(char(fmri_nii));
end

if isa(image_set, 'image_vector')
    maps = image_set;
    image_set_name = class(image_set);
    if isprop(maps, 'metadata_table') && ~isempty(maps.metadata_table) && ...
            any(strcmp('target', maps.metadata_table.Properties.VariableNames))
        map_names = maps.metadata_table.target;
    elseif isprop(maps, 'image_labels') && ~isempty(maps.image_labels)
        map_names = maps.image_labels;
    elseif ~isempty(maps.image_names)
        map_names = cellstr(maps.image_names);
    else
        map_names = arrayfun(@(i) sprintf('map_%03d', i), 1:size(maps.dat, 2), 'UniformOutput', false);
    end
else
    image_set_name = char(image_set);
    [maps, map_names] = load_image_set(char(image_set), 'noverbose');
end
map_names = cellstr(string(map_names));

if ~isempty(opts.MapNames)
    req = cellstr(string(opts.MapNames));
    [tf, idx] = ismember(req, map_names);
    idx = idx(tf);
    maps = get_wh_image(maps, idx);
    map_names = map_names(idx);
end
if isempty(map_names)
    error('No maps matched the requested ImageSet/MapNames.');
end

n_map = numel(map_names);
n_tp = size(dat.dat, 2);
TC = nan(n_tp, n_map);
for i = 1:n_map
    this_map = get_wh_image(maps, i);
    y = apply_mask(dat, this_map, 'pattern_expression', 'ignore_missing', char(opts.SimilarityMetric));
    y = y(:);
    TC(:, i) = local_zscore(y);
end

meta = struct();
meta.available_signatures = map_names;
meta.image_set = image_set_name;
meta.similarity_metric = char(opts.SimilarityMetric);
end

function y = local_zscore(y)
y = y(:);
s = std(y);
if s == 0 || isnan(s)
    y = zeros(size(y));
else
    y = (y - mean(y)) ./ s;
end
end
