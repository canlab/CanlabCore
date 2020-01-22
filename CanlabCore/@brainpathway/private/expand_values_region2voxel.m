function output_dat = expand_values_region2voxel(obj, input_dat)
% Map region-wise values to voxel values
% output_dat = expand_values_region2voxel(obj, input_dat)
%
% obj: A brainpathway object
% input_dat: regions x images values to expand
% output_dat: voxels x images
%
% Tor Wager, 12/2019

% Integer codes for which voxels belong to each region
voxel_indx = obj.region_atlas.dat;

[n_regions, n_images] = size(input_dat);

output_dat = zeros([size(voxel_indx, 1), n_images]);


for i = 1:n_regions
    
    wh = voxel_indx == i;
    
    for j = 1:n_images
        
        myvalue = input_dat(i, j);
        
        output_dat(wh, j) = myvalue;
        
    end
    
end

end % function
