function [voldata, vectorized_voldata, xyz_coord_struct] = reconstruct_image(obj)
% Reconstruct a 3-D or 4-D image from image_vector object obj
%
% voldata is and X x Y x Z x Images matrix
% vectorized_voldata is the same, with all voxels vectorized
%
% This output has one element for every voxel in THE ENTIRE IMAGE, and so
% can be very memory-intensive.  But it's useful for lining up voxels
% across images with different masks/in-mask voxels.
%
% This function returns output in memory;
% see image_vector.write for writing .img files to disk.
%
% :Outputs:
%
%   **voldata:**
%        3-D recon volume
%   **vectorized_voldata:**
%        volume in column vetor, iimg_xxx function format
%   **xyz_coord_struct:**
%        has fields with coordinate information in mm (world) space
%          - x, y, z : vectors of coordinates in mm for each of the 3
%            dimensions of the image
%          - X, Y, Z : output matrices from meshgrid with mm coordinates,
%            for volume visualization.
%            These can be passed to surf or isocaps functions for volume
%            visualization in world space (mm).
%
% ..
%    Copyright 2011 tor wager
%
%    Programmers' notes:
%
%    Aug 2012: This function does not flip the data based on the sign of x dimension.  
%    The flipping is applied in image writing / display in
%    iimg_reconstruct_vols, write method, spm_orthviews, etc.
%
%    July 2013 : Tor : Added xyz_coord_struct output
% ..

obj = replace_empty(obj);
voldata = iimg_reconstruct_vols(obj.dat, obj.volInfo);

if nargout > 1
    vectorized_voldata = voldata(:);
end

if nargout > 2
    disp('Returning coordinates in mm and meshgrid matrices.');
    
    % return xyz mm coordinates for full volume
    xyzmin = voxel2mm([1 1 1]', obj.volInfo.mat);
    xyzmax = voxel2mm(obj.volInfo.dim', obj.volInfo.mat); 
   
    %voxsize = diag(obj.volInfo.mat(1:3, 1:3));
    
    % rotate and flip for compatibility with addbrain.m direction

    for i = 1:size(voldata, 3)
        
        vout(:, :, i) = rot90(voldata(:, :, i));
        
        vout(:, :, i) = flipdim(vout(:, :, i), 1);
    end
    
    dims = size(voldata);

    % reverse dims 2 and 1 because dim 1 is y
    x = linspace(xyzmin(1), xyzmax(1), dims(1))';
    y = linspace(xyzmin(2), xyzmax(2), dims(2))';
    z = linspace(xyzmin(3), xyzmax(3), dims(3))';
    
    [X, Y, Z] = meshgrid(x, y, z);
    
    xyz_coord_struct = struct('voldata', vout, 'x', x, 'y', y, 'z', z, 'X', X, 'Y', Y, 'Z', Z);

end
