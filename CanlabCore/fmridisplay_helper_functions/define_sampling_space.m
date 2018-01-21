function SPACE = define_sampling_space(V, varargin)
% Define the sampling space of an image, with an upsampled space to 0.5 mm
% resolution
%
% :Usage:
% ::
%
%     SPACE = define_sampling_space(V, [upsamplefactor])
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
% :Outputs:
%
%   **Xo, Yo:**
%        Meshgrid for original voxel space
%
%   **X, Y:**
%        Meshgrid for upsampled voxel space at 0.5 mm resolution
%
%   **Xmm, Ymm:**
%        Meshgrid for upsampled space in mm 
%
%   **xcoords, ycoords:**
%        mm coordinates for rows and cols for slice locations
%
%   **new_voxSize:**
%        new voxel size in mm for upsampled space
%
%   **usfactor:**
%        Upsampling factor for new sampleing space
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

V = V(1);
voxSize = abs(diag(V.mat(1:3, 1:3))); % in mm

%Zo = vol3(:, :, 32); % x by y mm

% nx = size(vol3, 2);
% ny = size(vol3, 1);

nx = V.dim(1);
ny = V.dim(2);
nz = V.dim(3);

[Xo, Yo, Zo] = meshgrid(1:ny, 1:nx, 1:nz);

if length(varargin) > 0 && ~isempty(varargin{1})
    usfactor = varargin{1};
else
    usfactor = ceil(max(voxSize) ./ 1);  % upsample factor for 1 mm res
end

new_voxSize = voxSize ./ usfactor;

x = linspace(1, nx, nx * usfactor);
y = linspace(1, ny, ny * usfactor);
z = linspace(1, nz, nz * usfactor);

[X, Y, Z] = meshgrid(y, x, z);

bottomleft_mm = V.mat(1:3, 4);
topright_mm = voxel2mm([V.dim(1:3)]', V.mat);

xcoords = linspace(bottomleft_mm(1,1), topright_mm(1), nx * usfactor);
ycoords = linspace(bottomleft_mm(2,1), topright_mm(2), ny * usfactor);
zcoords = linspace(bottomleft_mm(3,1), topright_mm(3), nz * usfactor);

[Xmm, Ymm, Zmm] = meshgrid(ycoords, xcoords, zcoords);

SPACE = struct('V', V, 'Xo', Xo, 'Yo', Yo, 'Zo', Zo, 'X', X, 'Y', Y, 'Z', Z, ...
    'usfactor', usfactor, 'new_voxSize', new_voxSize, ...
    'xcoords', xcoords, 'ycoords', ycoords, 'zcoords', zcoords, 'Xmm', Xmm, 'Ymm', Ymm, 'Zmm', Zmm);

end
