function xyzmm = get_xyzmm_coordinates(obj)
% Extract non-empty, inmask coordinates for all voxels in an image object

xyzvox = obj.volInfo.xyzlist;

obj = remove_empty(obj);

xyzvox(obj.removed_voxels, :) = [];

xyzmm = voxel2mm(xyzvox', obj.volInfo.mat)';

end

