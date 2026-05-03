function xyzmm = get_xyzmm_coordinates(obj)
% get_xyzmm_coordinates Return MNI/world (mm) coordinates for non-empty in-mask voxels.
%
% Extract non-empty, in-mask coordinates for all voxels in an image
% object and convert from voxel indices to mm using the affine matrix
% in obj.volInfo.mat.
%
% :Usage:
% ::
%
%     xyzmm = get_xyzmm_coordinates(obj)
%
% :Inputs:
%
%   **obj:**
%        An image_vector (or subclass) object with a populated
%        volInfo.xyzlist and volInfo.mat. Empty voxels are removed
%        before coordinate conversion.
%
% :Outputs:
%
%   **xyzmm:**
%        An [n x 3] matrix of mm coordinates (one row per non-empty
%        in-mask voxel).
%
% :Examples:
% ::
%
%     xyz = get_xyzmm_coordinates(dat);
%
% :See also:
%   - voxel2mm
%   - remove_empty
%   - reconstruct_image

xyzvox = obj.volInfo.xyzlist;

obj = remove_empty(obj);

xyzvox(obj.removed_voxels, :) = [];

xyzmm = voxel2mm(xyzvox', obj.volInfo.mat)';

end

