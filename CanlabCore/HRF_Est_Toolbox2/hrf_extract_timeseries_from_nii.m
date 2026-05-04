function [tc, TR, n_tp] = hrf_extract_timeseries_from_nii(fmri_nii, mask_nii)
%HRF_EXTRACT_TIMESERIES_FROM_NII Mean fMRI timeseries from 4D NIfTI.
info = niftiinfo(fmri_nii);
V = double(niftiread(info));
if ndims(V) ~= 4
    error('Expected 4D fMRI image, got %d dimensions.', ndims(V));
end
n_tp = size(V, 4);

TR = [];
if isfield(info, 'PixelDimensions') && numel(info.PixelDimensions) >= 4
    TR = info.PixelDimensions(4);
end
if isempty(TR) || TR <= 0
    error('Could not infer TR from NIfTI header. Pass ''TR'' explicitly.');
end

if nargin < 2 || isempty(mask_nii)
    mask = true(size(V,1), size(V,2), size(V,3));
else
    M = niftiread(mask_nii);
    mask = M > 0;
    if ~isequal(size(mask), size(V(:,:,:,1)))
        error('Mask dimensions must match fMRI spatial dimensions.');
    end
end

vox_by_time = reshape(V, [], n_tp);
mask_lin = mask(:);
masked = vox_by_time(mask_lin, :);
tc = mean(masked, 1)';
tc = (tc - mean(tc)) ./ std(tc);
end
