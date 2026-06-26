function obj = remove_empty(obj, varargin)
% remove_empty (fmri_surface_data): no-op. Grayordinate data is already compact.
%
% :Usage:
% ::
%     obj = remove_empty(obj)
%
% Unlike image_vector/fmri_data -- which store sparse 3-D volumes and drop
% all-zero / NaN voxels to save space -- fmri_surface_data holds CIFTI
% grayordinate data, which is already compact (the medial wall is excluded by
% construction and there are no out-of-brain voxels). The .dat matrix is ALWAYS
% the full [nGrayordinates x nMaps] set, in 1:1 row correspondence with
% .brain_model, and is never squeezed.
%
% This override therefore returns the object UNCHANGED. It exists so that (a)
% inherited methods that call remove_empty internally still compose, and (b)
% user code that calls remove_empty by habit gets a valid object back. The
% .removed_voxels / .removed_images vectors stay all-false.
%
% See docs/fmri_surface_data_design_plan.md (decision D5b).
%
% :See also: replace_empty, fmri_surface_data, get_wh_image

% Intentionally a no-op. (varargin accepted and ignored for signature parity
% with image_vector/remove_empty.)
end
