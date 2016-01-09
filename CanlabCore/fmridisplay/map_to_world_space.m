function SPACE = map_to_world_space(V)
% :Usage:
% ::
%
%     SPACE = map_to_world_space(V)
%
% :Inputs:
%
%   **V:**
%        spm-style .mat structure, e.g., from spm_vol
%
%   **V.mat:**
%        4 x 4 matrix of voxel sizes and mm coords for the bottom
%        back left vox
%
%   **V.dim:**
%        dimensions of image
%
% :Outputs: SPACE structure, with fields:
%
%   **Xmm, Ymm, Zmm:**
%        Meshgrid for voxel volume in mm space
%
%   **xcoords, ycoords, zcoords:**
%        mm coordinates for rows, cols, slices
%
% ..
%    Tor Wager, Feb 2011
% ..

bottomleft_mm = V.mat(1:3, 4);
topright_mm = voxel2mm([V.dim(1:3)]', V.mat);

xcoords = linspace(bottomleft_mm(1), topright_mm(1), V.dim(1));
ycoords = linspace(bottomleft_mm(2), topright_mm(2), V.dim(2));
zcoords = linspace(bottomleft_mm(3), topright_mm(3), V.dim(3));

[Xmm, Ymm, Zmm] = meshgrid(ycoords, xcoords, zcoords);



SPACE = struct('V', V, 'Xmm', Xmm, 'Ymm', Ymm, 'Zmm', Zmm, ...
    'xcoords', xcoords, 'ycoords', ycoords, 'zcoords', zcoords);


end

