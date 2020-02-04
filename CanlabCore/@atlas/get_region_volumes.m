function [vol_in_cubic_mm, voxcount, imtx, parcel_index_values] = get_region_volumes(atlas_obj)
% Get the volume (and raw voxel count) of each region in an atlas object
%
% [vol_in_cubic_mm, voxcount, imtx, parcel_index_values] = get_region_volumes(atlas_obj)
%
% imtx : indicator matrix, 1/0, for each region

[n_regions, n_regions_with_data, missing_regions] = num_regions(atlas_obj);

dat = double(atlas_obj.dat);
[imtx, parcel_index_values] = condf2indic(dat, 'integers');
imtx(isnan(imtx)) = 0;

% if last regions are missing, condf2indic will be missing last entries
if any(missing_regions == n_regions)
    [~,f] = size(imtx);
    imtx(:, (f+1):n_regions) = 0;
    parcel_index_values = [parcel_index_values; n_regions];
end

voxcount = double(sum(imtx));

% count per cubic mm
vox_volume_in_mm = prod(abs(diag(atlas_obj.volInfo.mat(1:3, 1:3))));

vol_in_cubic_mm = voxcount .* vox_volume_in_mm;


end % function


