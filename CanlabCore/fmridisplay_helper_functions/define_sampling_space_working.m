function SPACE = define_sampling_space(V, varargin)
% Build mm-coordinate sampling grids for interpolation of NIfTI volumes.
%
% This function is primarily used in:
% 1 - Resampling spaces from one coordinate system to
%   another so that voxels align.  primarily image_vector/resample_space
%   Uses .Xmm, .Ymm, .Zmm coordinates
%   resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, varargin{:});
%
% 2 - display, e.g., display_slice, fmridisplay/render_blobs
%   This uses SPACE.X and SPACE.Y fields.
%
% :Usage:
% ::
%
%     SPACE = define_sampling_space(V, [desired_res_mm])
%
% :Inputs:
%
%   **V:**
%        spm-style .mat structure, e.g., from spm_vol
%        SPM-style volume struct with fields:
%           V.mat : 4x4 affine voxel->mm transform
%           V.dim : [nx ny nz]
%
%   **V.mat:**
%        4 x 4 matrix of voxel sizes and mm coords for the bottom
%        back left vox
%
%   **V.dim:**
%        dimensions of image
%
% :Outputs:
%
% SPACE : struct with fields
%     .Xmm .Ymm .Zmm : mm coordinate grids for interpolation
%     .X .Y .Z : voxel index grids
%
%   **Xo, Yo:**
%        Meshgrid for original voxel space, Matlab indexing 1...i, 1...j
%
%   **X, Y, Z:**
%        Meshgrid for upsampled voxel space, 1 mm res
%
%   **Xmm, Ymm, Zmm:**
%        Meshgrid for upsampled space in mm, 1 mm res 
%
%   **xcoords, ycoords, zcoords:**
%        mm coordinates for rows and cols for slice locations
%
%   **new_voxSize:**
%        new voxel size in mm for upsampled space
%
%   **usfactor:**
%        Upsampling factor for new sampling space, x y z triplet
%
% :Examples:
% ::
%
%    overlay = which('SPM8_colin27T1_seg.img');  % spm8 seg cleaned up
%    V = spm_vol(overlay);
%    SPACE = define_sampling_space(V)
%
%    % Define mm sampling space in original voxel coord resolution
%    SPACE = define_sampling_space(V, 1)
% 
%    original (o) and new (X, Y) grid space
%    xcoords, ycoords: mm coords centered on origin


% Programmers' notes:
% 2026-3-8 Tor Wager: This version fixes a voxel-center convention bug that caused a
% +1 voxel shift in x and y when used with resample_space().
% Removed bounding-box calculation
% bottomleft_mm = V.mat(1:3, 4);
% topright_mm = voxel2mm([V.dim(1:3)]', V.mat);
%
% voxel2mm([dim]) implicitly treats the last voxel index as dim, while NIfTI voxel indices correspond to centers 0..dim-1 (or equivalently 1..dim but offset by 0.5 depending on convention).
% Using [dim] instead of [dim-1] shifts the grid when the mm sampling lattice is later constructed.
% 
% The safest solution is to compute mm coordinates directly from voxel indices using the affine, ensuring consistent voxel-center conventions.
% 
% bottomleft_mm = V.mat(1:3, 4);
% topright_mm = voxel2mm([V.dim(1:3)]', V.mat); % x_world = M * [x_i y_i z_i 1]'; 
% 
% xcoords = linspace(bottomleft_mm(1,1), topright_mm(1), nx * usfactor);
% ycoords = linspace(bottomleft_mm(2,1), topright_mm(2), ny * usfactor);
% zcoords = linspace(bottomleft_mm(3,1), topright_mm(3), nz * usfactor);
% 
% [Xmm, Ymm, Zmm] = meshgrid(ycoords, xcoords, zcoords);

if length(V) > 1
    warning('space-defining structure V should not be a vector! This may apply incorrect transformations if space is different for different images indexex in V')
    V = V(1);
end

% ---------------------------------------------------------------------
% dimensions
% ---------------------------------------------------------------------

dim = V.dim(1:3);

% ---------------------------------------------------------------------
% voxel index grids (MATLAB indexing)
% ---------------------------------------------------------------------
% Xvox, Yvox, Zvox in original voxel coordinates
% Matlab meshgrid order is V(row, column, slice) = V(y,x,z)
% so must exchange dims 1 and 2

[Xo, Yo, Zo] = meshgrid(1:dim(2), 1:dim(1), 1:dim(3));

% ---------------------------------------------------------------------
% convert to voxel-center coordinates expected by NIfTI affine
%
% NIfTI affine assumes voxel centers at integer coordinates starting
% at 0 (0..dim-1). MATLAB arrays start at 1.
% BUT SPM starts at 1, 1, 1 too so do not subtract 1!!
% ---------------------------------------------------------------------

ijk = [Xo(:) Yo(:) Zo(:)]';

% ---------------------------------------------------------------------
% convert voxel indices -> mm coordinates
% ---------------------------------------------------------------------

XYZmm = V.mat * [ijk; ones(1, size(ijk,2))];

Xmm = reshape(XYZmm(1,:), dim);
Ymm = reshape(XYZmm(2,:), dim);
Zmm = reshape(XYZmm(3,:), dim);

xcoords = reshape(XYZmm(1,:), dim);
ycoords = reshape(XYZmm(2,:), dim);
zcoords = reshape(XYZmm(3,:), dim);

% ---------------------------------------------------------------------
% Up-sampled version for display
% ---------------------------------------------------------------------

% usfactor should be 3-element for x, y, z dims
if ~isempty(varargin) && ~isempty(varargin{1})
    usfactor = varargin{1}; % desired voxel res in mm for SPACE.X, SPACE.Y, SPACE.Z
else
    usfactor = 1;  % upsample factor for 1 mm res
end

voxSize = abs(diag(V.mat(1:3, 1:3)));

usfactor = voxSize ./ usfactor; % in mm
new_voxSize = usfactor; % 1 mm

x = linspace(1, dim(1), dim(1) * usfactor(1));
y = linspace(1, dim(2), dim(2) * usfactor(2));
z = linspace(1, dim(3), dim(3) * usfactor(3));

[X, Y, Z] = meshgrid(y, x, z); % reverse y and x 


SPACE = struct('V', V, 'Xo', Xo, 'Yo', Yo, 'Zo', Zo, 'X', X, 'Y', Y, 'Z', Z, ...
    'usfactor', usfactor, 'new_voxSize', new_voxSize, ...
    'xcoords', xcoords, 'ycoords', ycoords, 'zcoords', zcoords, 'Xmm', Xmm, 'Ymm', Ymm, 'Zmm', Zmm);

end
