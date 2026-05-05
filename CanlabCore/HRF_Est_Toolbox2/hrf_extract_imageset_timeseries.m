function [TC, meta] = hrf_extract_imageset_timeseries(fmri_nii, image_set, varargin)
%HRF_EXTRACT_IMAGESET_TIMESERIES Apply load_image_set maps to each fMRI volume.
% Supports Buckner networks, Margulies gradients, Hansen receptor maps, etc.

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x));
p.addRequired('image_set', @(x) ischar(x) || isstring(x));
p.addParameter('MapNames', {}, @(x) iscell(x) || isstring(x));
p.addParameter('SimilarityMetric', 'dotproduct', @(x) ischar(x) || isstring(x));
p.parse(fmri_nii, image_set, varargin{:});
opts = p.Results;

dat = fmri_data(char(fmri_nii));
[maps, map_names] = load_image_set(char(image_set), 'noverbose');
map_names = cellstr(string(map_names));

if ~isempty(opts.MapNames)
    req = cellstr(string(opts.MapNames));
    [tf, idx] = ismember(req, map_names);
    idx = idx(tf);
    maps = get_wh_image(maps, idx);
    map_names = map_names(idx);
end

n_map = numel(map_names);
n_tp = size(dat.dat, 2);
TC = nan(n_tp, n_map);
for i = 1:n_map
    this_map = get_wh_image(maps, i);
    y = apply_mask(dat, this_map, 'pattern_expression', 'ignore_missing', char(opts.SimilarityMetric));
    y = y(:);
    TC(:, i) = (y - mean(y)) ./ std(y);
end

meta = struct();
meta.available_signatures = map_names;
meta.image_set = char(image_set);
meta.similarity_metric = char(opts.SimilarityMetric);
end
