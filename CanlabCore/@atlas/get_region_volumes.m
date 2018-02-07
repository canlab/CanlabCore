function [vol_in_cubic_mm, voxcount] = get_region_volumes(atlas_obj)
% Get the volume (and raw voxel count) of each region in an atlas object
%
% [vol_in_cubic_mm, voxcount] = get_region_volumes(atlas_obj)
%


dat = double(atlas_obj.dat);
[imtx, parcels_we_have] = condf2indic(dat, 'integers');
imtx(isnan(imtx)) = 0;

voxcount = double(sum(imtx));

% count per cubic mm
vox_volume_in_mm = prod(abs(diag(atlas_obj.volInfo.mat(1:3, 1:3))));

vol_in_cubic_mm = voxcount .* vox_volume_in_mm;


end % function


