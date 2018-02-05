function [vol_in_cubic_mm, voxcount] = get_region_volumes(atlas_obj)
% Get the volume (and raw voxel count) of each region in an atlas object
%
% [vol_in_cubic_mm, voxcount] = get_region_volumes(atlas_obj)
%


[imtx, ivals] = condf2indic(atlas_obj.dat);

wh = find(ivals == 0);
ivals(wh) = [];
imtx(wh, :) = [];

voxcount = double(sum(imtx));

% count per cubic mm
vox_volume_in_mm = prod(abs(diag(atlas_obj.volInfo.mat(1:3, 1:3))));

vol_in_cubic_mm = voxcount .* vox_volume_in_mm;


end % function


