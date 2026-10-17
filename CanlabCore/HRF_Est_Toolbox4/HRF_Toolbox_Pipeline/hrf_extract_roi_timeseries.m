function [TC, meta] = hrf_extract_roi_timeseries(fmri_nii, atlas_obj, varargin)
%HRF_EXTRACT_ROI_TIMESERIES Extract ROI-mean time series from atlas regions.

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x));
p.addRequired('atlas_obj', @(x) isa(x, 'atlas'));
p.addParameter('Regions', {}, @(x) iscell(x) || isstring(x));
p.parse(fmri_nii, atlas_obj, varargin{:});
opts = p.Results;

dat = fmri_data(char(fmri_nii));

if isempty(opts.Regions)
    regions = atlas_obj.labels;
else
    regions = cellstr(string(opts.Regions));
end

n_reg = numel(regions);
n_tp = size(dat.dat, 2);
TC = nan(n_tp, n_reg);
for r = 1:n_reg
    at_sub = atlas_obj.select_atlas_subset(regions(r), 'exact');
    y = mean(apply_mask(dat, at_sub).dat, 1);
    y = y(:);
    TC(:, r) = (y - mean(y)) ./ std(y);
end

meta = struct();
meta.available_signatures = regions;
meta.image_set = 'atlas_rois';
meta.similarity_metric = 'mean_signal';
end
