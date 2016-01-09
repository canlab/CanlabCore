function [resampled_dat, SPACEto] = resample_space(dat, V, targetsp)
% :Usage:
% ::
%
%     [resampled_dat, SPACEto] = resample_space(dat, V, [target V or target SPACE])
%
% :Inputs:
%
%   **imdat:**
%        3-D volume data
%
%   **V:**
%        spm-style .mat structure for dat, e.g., from spm_vol
%
%   **V.mat:**
%        4 x 4 matrix of voxel sizes and mm coords for the bottom
%                   back left vox
%
%   **V.dim:**
%        dimensions of image
%
%
%   **targetsp:**
%
%      **target V:**
%          spm-style .mat structure defining space to transform to
%
%      -- OR --
%
%      **target SPACE:**
%          target SPACE, with Xmm, Ymm, Zmm; see map_to_world_space.m
%
% :Outputs:
%
%   **resampled_dat:**
%        data sampled in new space
%
%   **SPACE structure, with fields:**
%
%      **Xmm, Ymm, Zmm:**
%          Meshgrid for voxel volume in mm space
%
%      **xcoords, ycoords, zcoords:**
%          mm coordinates for rows, cols, slices
%
%      **V:**
%          Vol info structure for image in new space

SPACE = map_to_world_space(V);

if isfield(targetsp, 'Xmm') && isfield(targetsp, 'Ymm') && isfield(targetsp, 'Zmm')
    SPACEto = targetsp;

elseif isfield(targetsp, 'mat') && isfield(targetsp, 'dim')
    SPACEto = map_to_world_space(Vto);

else
    error('targetsp is not valid; use spm_vol.m or map_to_world_space.m to define this.');
end

resampled_dat = interp3(SPACE.Xmm, SPACE.Ymm, SPACE.Zmm, dat, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm);

end

