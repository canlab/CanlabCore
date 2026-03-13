function [voldata, vectorized_voldata, xyz_coord_struct] = reconstruct_image(obj)
% Reconstruct a 3-D or 4-D image from image_vector object obj
%
% [voldata, vectorized_voldata, xyz_coord_struct] = reconstruct_image(obj)
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
%        volume in column vector, iimg_xxx function format
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

if isa(obj, 'statistic_image')
    % Apply threshold
    obj.dat = obj.dat .* logical(obj.sig);
end

voldata = iimg_reconstruct_vols(obj.dat, obj.volInfo);

if nargout > 1
    vectorized_voldata = voldata(:);
end

if nargout > 2
    % disp('Returning coordinates in mm and meshgrid matrices.');
    
    % 3/12/2026 Tor Wager - best way is to use affine transformation matrix
    % M rather than bounding box method

    % -------------------------------------------------------------------------
    % Create grid in mm coordinates (world) for surface mapping
    % -------------------------------------------------------------------------
    xd = obj.volInfo.dim(1);
    yd = obj.volInfo.dim(2);
    zd = obj.volInfo.dim(3);
    [Xt, Yt, Zt] = meshgrid(1:xd, 1:yd, 1:zd); % length(y) × length(x) × length(z)
    ijk_vox = [Xt(:)'; Yt(:)'; Zt(:)'; ones(1,numel(Xt))];
    ijk_xyzmm = obj.volInfo.mat * ijk_vox;                  % where M is affine matrix
    X = reshape(ijk_xyzmm(1,  :)', [yd xd zd]);
    Y = reshape(ijk_xyzmm(2,  :)', [yd xd zd]);
    Z = reshape(ijk_xyzmm(3,  :)', [yd xd zd]);

    % reverse y and x to match meshgrid
    y = squeeze(X(:,1,1));
    x = squeeze(Y(1,:,1))';
    z = squeeze(Z(1,1,:));

    % rotate and flip for compatibility with addbrain.m direction
    vout = zeros(yd, xd, zd);
    if ndims(voldata) > 3, vtmp = voldata(:, :, :, 1); else, vtmp = voldata; end

    for i = 1:size(voldata, 3)

        vout(:, :, i) = rot90(vtmp(:, :, i));

        vout(:, :, i) = flipdim(vout(:, :, i), 1);
    end

    % render_on_surface uses
    % c = interp3(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, mesh_struct.voldata, ...
    % sh.Vertices(:,1), sh.Vertices(:,2), sh.Vertices(:,3), interp);

    xyz_coord_struct = struct('voldata', vout, 'x', x, 'y', y, 'z', z, 'X', X, 'Y', Y, 'Z', Z);

    % % -------------------------------------------------------------------------
    % % Old bounding box method
    % % -------------------------------------------------------------------------
    % 
    % % disp('Returning coordinates in mm and meshgrid matrices.');
    % 
    % % return xyz mm coordinates for full volume
    % % note: Nifti starts at [0 0 0] but SPM starts at [1 1 1]
    % xyzmin = voxel2mm([1 1 1]', obj.volInfo.mat);
    % xyzmax = voxel2mm(obj.volInfo.dim', obj.volInfo.mat); 
    % 
    % % rotate and flip for compatibility with addbrain.m direction
    % 
    % for i = 1:size(voldata, 3)
    % 
    %     vout(:, :, i) = rot90(voldata(:, :, i));
    % 
    %     vout(:, :, i) = flipdim(vout(:, :, i), 1);
    % end
    % 
    % dims = size(voldata);
    % 
    % % reverse dims 2 and 1 because dim 1 is y
    % % * NOTE: x and y probably need to be reversed here, but needs thorough testing to make sure everything works *
    % x = linspace(xyzmin(1), xyzmax(1), dims(1))';
    % y = linspace(xyzmin(2), xyzmax(2), dims(2))';
    % z = linspace(xyzmin(3), xyzmax(3), dims(3))';
    % 
    % [X, Y, Z] = meshgrid(x, y, z);
    % 
    % xyz_coord_struct = struct('voldata', vout, 'x', x, 'y', y, 'z', z, 'X', X, 'Y', Y, 'Z', Z);

end
