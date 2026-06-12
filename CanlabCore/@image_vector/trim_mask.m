function obj = trim_mask(obj)
% trim_mask Exclude empty voxels from the mask in obj.volInfo and rebuild .volInfo.
%
% Removes voxels whose values across all images are zero or NaN, and
% updates the .volInfo bookkeeping (image_indx, wh_inmask, xyzlist,
% cluster, n_inmask) accordingly. For fmri_data objects, the embedded
% .mask is updated in parallel.
%
% :Usage:
% ::
%
%     obj = trim_mask(obj)
%
% :Inputs:
%
%   **obj:**
%        An image_vector / fmri_data object.
%
% :Outputs:
%
%   **obj:**
%        The input with all-zero / all-NaN voxels removed from .dat and
%        from .volInfo (and from .mask for fmri_data objects).
%
% :Examples:
% ::
%
%     obj = trim_mask(obj);
%
% :See also:
%   - remove_empty
%   - replace_empty
%   - rebuild_volinfo_from_dat
%
% ..
%    Tor Wager, 2013
% ..

obj = replace_empty(obj);
whomit = all(obj.dat == 0 | isnan(obj.dat), 2);

obj.dat(whomit, :) = [];
if ~isempty(obj.removed_voxels) && length(obj.removed_voxels) == obj.volInfo.n_inmask
    obj.removed_voxels(whomit) = [];
end
obj.volInfo.image_indx(obj.volInfo.wh_inmask(whomit)) = 0;
obj.volInfo.wh_inmask(whomit) = [];
obj.volInfo.xyzlist(whomit, :) = [];
obj.volInfo.cluster(whomit) = [];
obj.volInfo.n_inmask = length(obj.volInfo.wh_inmask);

if isa(obj, 'fmri_data')
    obj.mask.dat(whomit, :) = [];
    if ~isempty(obj.mask.removed_voxels) && length(obj.mask.removed_voxels) == obj.mask.volInfo.n_inmask
        obj.mask.removed_voxels(whomit) = [];
    end
    obj.mask.volInfo.image_indx(obj.mask.volInfo.wh_inmask(whomit)) = 0;
    obj.mask.volInfo.wh_inmask(whomit) = [];
    obj.mask.volInfo.xyzlist(whomit, :) = [];
    obj.mask.volInfo.cluster(whomit) = [];
    obj.mask.volInfo.n_inmask = length(obj.mask.volInfo.wh_inmask);
end

end % function
