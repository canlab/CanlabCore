function [obj_resamp, space_mapper_struct] = resample_space(obj, target_obj, varargin)
% Resample image data into the voxel space of another image object.
%
% :Usage:
% ::
%
%     [obj_resamp, space_mapper_struct] = resample_space(obj, target_obj, [interp_method, extrapval])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2026  Tor Wager - coding assistance from ChatGPT
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **obj:**
%        Source image object. Expected to be an image_vector-compatible
%        object with fields .dat and .volInfo.
%
%   **target_obj:**
%        Target image object defining the output voxel space. Expected to
%        contain .volInfo and, for masked output, .volInfo.wh_inmask.
%
% :Optional Inputs:
%
%   **interp_method:**
%        Interpolation method passed to interp3, e.g., 'linear', 'nearest',
%        'cubic', or 'spline'. If omitted, MATLAB uses interp3 defaults
%        unless additional varargin entries are supplied.
%
%   **extrapval:**
%        Optional extrapolation value passed through to interp3.
%
% :Outputs:
%
%   **obj_resamp:**
%        Resampled image object in the voxel space of target_obj.
%
%   **space_mapper_struct:**
%        Structure containing intermediate variables used for resampling,
%        including source and target voxel coordinates, target mm
%        coordinates, source-space query coordinates, dimensions, and
%        affine matrices.
%
% :Examples:
% ::
%
%    obj = load_image_set('emotionreg');
%    m = mean(obj);
%    template_obj = fmri_data(which('brainmask.nii'));
%    [obj_resamp, space_mapper_struct] = resample_space2(obj, template_obj, 'linear', 0);
%    [obj_resamp, space_mapper_struct] = resample_space2(obj, template_obj, 'spline');
%    m2 = mean(obj_resamp);
%    orthviews_multiple_objs({m m2});
%
% obj = load_image_set('emotionreg');
% npsvals = apply_mask(obj, nps, 'pattern_expression');
% npsvals2 = apply_mask(resample_space(obj, nps), nps, 'pattern_expression');
% npsvals3 = apply_mask(obj, resample_space(nps, obj), 'pattern_expression');
% npsvals4 = apply_mask(resample_space(obj, nps, 'spline'), nps, 'pattern_expression');
% corr([npsvals npsvals2 npsvals3 npsvals4])
% 
% :References:
%
% :See also:
%   - image_vector/resample_space
%   - reconstruct_image
%   - interp3

% ..
%    Programmers' notes:
%    2026-03-08: Renamed voxel-coordinate variables for clarity:
%    x_vox_orig/y_vox_orig/z_vox_orig are source-grid voxel coordinates;
%    x_vox_target/y_vox_target/z_vox_target are continuous source-space
%    query coordinates corresponding to target voxel centers.
%
%    Added second output space_mapper_struct with relevant intermediate variables.
%
%    Uses SPM-style affine inversion with left matrix division:
%    ijk_from = Vfrom.mat \ XYZmm_homogeneous
%
%    The final interp3 call intentionally uses (y, x, z) ordering because
%    MATLAB arrays are indexed as (row, column, slice), whereas the affine
%    logic here tracks dimension-1, dimension-2, and dimension-3 voxel
%    coordinates explicitly.
% ..

% -------------------------------------------------------------------------
% Define mapping from obj to target
% -------------------------------------------------------------------------

space_mapper_struct = define_space_mapping(obj, target_obj);

% NOTE:
% SPM does
% Mslice  = spm_matrix([0 0 j]);      % Matrix specifying this slice
% Mtrans  = Mto \ Mmap \ Mslice;          % Affine mappping mtx: Mask -> TOvol
% Mtrans = same as inv(Mto) * inv(Mmap) * Mslice
%
% v_source = inv(Mmap) * Mto * Mslice * target_pixels
%
% T = M_{source}^{-1} M_{target}
% T = A^{-1} * B --> maps space of A to space of B
% T = A \ B; --> maps space of B to space of A
% source_voxel = inv(Mfrom) * Mto * target_voxel
% T = Vfrom.mat \ Vto.mat;
% T = Msource \ Mtarget;
% v_source = T * v_target

% T = Vfrom.mat \ Vto.mat;
% ijk_target = [i j k 1]';
% ijk_source = T * ijk_target;

% -------------------------------------------------------------------------
% Interpolate each image
% -------------------------------------------------------------------------

nimgs = size(obj.dat, 2);
resampled_dat = zeros(prod(space_mapper_struct.dim_to), nimgs);

voldata = reconstruct_image(obj);

for i = 1:nimgs

    voldatai = squeeze(voldata(:, :, :, i));

    vals = interp3(space_mapper_struct.y_vox_orig, space_mapper_struct.x_vox_orig, space_mapper_struct.z_vox_orig, ...
        voldatai, ...
        space_mapper_struct.y_vox_query, space_mapper_struct.x_vox_query, space_mapper_struct.z_vox_query, varargin{:});

    resampled_dat(:, i) = vals(:);

end

% Save only in-mask values (in target mask)
resampled_dat = resampled_dat(space_mapper_struct.Vto.wh_inmask, :);

% -------------------------------------------------------------------------
% Create output object
% -------------------------------------------------------------------------

obj_resamp = replace_empty(obj);
obj_resamp.removed_voxels = 0;
obj_resamp.dat = single(resampled_dat);
obj_resamp.volInfo = space_mapper_struct.Vto;

% Replace mask for fmri_data objects
if isa(obj_resamp, 'fmri_data')
    obj_resamp.mask = target_obj;
end


end  % Main function

% Now stand-alone function
% function space_mapper_struct = define_space_mapping(obj, target_obj)
% 
% % -------------------------------------------------------------------------
% % Get volume info
% % -------------------------------------------------------------------------
% 
% Vfrom = obj.volInfo;
% Vto   = target_obj.volInfo;
% 
% dim_from = Vfrom.dim(1:3);
% dim_to   = Vto.dim(1:3);
% 
% % -------------------------------------------------------------------------
% % Create target voxel grid in MATLAB voxel coordinates
% % -------------------------------------------------------------------------
% 
% [Xt, Yt, Zt] = ndgrid(1:dim_to(1), 1:dim_to(2), 1:dim_to(3));
% 
% % Convert target voxel centers to NIfTI/SPM voxel-center coordinates
% % ijk_to = [Xt(:) - 1, Yt(:) - 1, Zt(:) - 1]';
% 
% % NOTE: don't do this
% % No extra origin offset is needed if you stay in SPM voxel 
% % coordinates throughout. SPM’s own documentation and old manual 
% % text describe the first voxel as (1,1,1), and SPM source code 
% % also contains explicit conversions only when exporting to formats 
% % that want a (0,0,0) origin.
% ijk_to = [Xt(:), Yt(:), Zt(:)]';
% 
% % [Xt,Yt,Zt] = ndgrid(1:dim_to(1), 1:dim_to(2), 1:dim_to(3));
% % ijk_target = [Xt(:)'; Yt(:)'; Zt(:)'; ones(1,numel(Xt))];
% 
% % -------------------------------------------------------------------------
% % Convert target voxel centers to mm coordinates
% % -------------------------------------------------------------------------
% 
% XYZmm_homogeneous = Vto.mat * [ijk_to; ones(1, size(ijk_to, 2))];
% XYZmm = XYZmm_homogeneous(1:3, :);
% 
% % -------------------------------------------------------------------------
% % Convert target mm coordinates into source voxel coordinates
% % -------------------------------------------------------------------------
% 
% ijk_from_homogeneous = Vfrom.mat \ XYZmm_homogeneous;
% ijk_from = ijk_from_homogeneous(1:3, :);
% 
% % Convert source-space query points back to MATLAB voxel coordinates
% % Note: Don't add 1 because we are in SPM affine space not nifti. first vox
% % is 1, 1, 1
% x_vox_target = ijk_from(1, :); % + 1;
% y_vox_target = ijk_from(2, :); % + 1;
% z_vox_target = ijk_from(3, :); % + 1;
% 
% % -------------------------------------------------------------------------
% % Define source voxel grid vectors
% % -------------------------------------------------------------------------
% 
% x_vox_orig = 1:dim_from(1);
% y_vox_orig = 1:dim_from(2);
% z_vox_orig = 1:dim_from(3);
% 
% % -------------------------------------------------------------------------
% % Collect diagnostic / intermediate variables
% % -------------------------------------------------------------------------
% 
% space_mapper_struct = struct();
% space_mapper_struct.Vfrom = Vfrom;
% space_mapper_struct.Vto = Vto;
% space_mapper_struct.dim_from = dim_from;
% space_mapper_struct.dim_to = dim_to;
% 
% % space_mapper_struct.Xt = Xt;
% % space_mapper_struct.Yt = Yt;
% % space_mapper_struct.Zt = Zt;
% 
% % space_mapper_struct.ijk_to = ijk_to;
% % space_mapper_struct.XYZmm_homogeneous = XYZmm_homogeneous;
% space_mapper_struct.XYZmm = XYZmm;
% 
% % space_mapper_struct.ijk_from_homogeneous = ijk_from_homogeneous;
% % space_mapper_struct.ijk_from = ijk_from;
% 
% space_mapper_struct.x_vox_orig = x_vox_orig;
% space_mapper_struct.y_vox_orig = y_vox_orig;
% space_mapper_struct.z_vox_orig = z_vox_orig;
% 
% space_mapper_struct.x_vox_target = x_vox_target;
% space_mapper_struct.y_vox_target = y_vox_target;
% space_mapper_struct.z_vox_target = z_vox_target;
% 
% space_mapper_struct.interp3_grid_order = {'y_vox_orig', 'x_vox_orig', 'z_vox_orig'};
% space_mapper_struct.interp3_query_order = {'y_vox_target', 'x_vox_target', 'z_vox_target'};
% 
% end % define_space_mapping

