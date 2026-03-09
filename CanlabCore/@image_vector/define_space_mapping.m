function space_mapper_struct = define_space_mapping(obj, target_obj)
% Define mapping of voxel spaces from obj to target
%
% :Usage:
% ::
%
%     space_mapper_struct = define_space_mapping(obj, target_obj)
%
% There are two basic steps to resampling:
%   1. Determine transform
%    = affine transform source->target, T = inv(Msource) * Mtarget
%   ijk_target = grid of voxels in target space
%   ijk_query = voxel query points to sample the same mm coordinates as in target space
%
% 2. Interpolate to resample in target space
%   ijk_target -> apply T -> ijk_query
%   interp(ijk_orig, voldata_orig, ijk_query)
%
% This function implements Step 1, saving relevant variables in a structure
% called space_mapper_struct that allows for application of the transform
% to resample any volume data in the same source space to the target space.
%
% This function computes the voxel-to-voxel coordinate transform needed
% to resample data from the voxel space of **obj** into the voxel space
% of **target_obj**.
%
% The transform is defined using the affine matrices stored in the
% volume information structures:
%
%     T = Vfrom.mat \ Vto.mat
%
%     This is equivalent to (but numerically more stable than)
%     T = inv(Vfrom.mat) * Vto.mat
%     T = inv(Msource) * Mtarget
%
% where
%
%     Vfrom.mat maps source voxel coordinates → world (mm)
%     Vto.mat   maps target voxel coordinates → world (mm)
%
% The resulting matrix **T** maps target voxel coordinates directly
% into source voxel coordinates:
%
%     v_query = T * v_target
%     v_query is also called v_source or ijk_source in code, as it maps voxels in the source
%     space at the query points needed to sample to the target space.
%     These are source voxel coordinates corresponding to target voxel centers
%
% This allows interpolation of the source image at locations corresponding
% to voxel centers in the target space.
%
% The function computes and stores all coordinate grids required to
% perform interpolation.
%
% The key fields in define_space_mapping that you need to interpolate data are:
% x_vox_orig, y_vox_orig, z_vox_orig -> grid of voxels in the original space
% x_vox_query, y_vox_query, z_vox_query -> grid of query points (relative to original space) for interpolation
% dim_to -> 3-D dimensions of target space for reshaping into 3-D matrix
%
% After resampling (done in resample_space, which calls this function),
% you additionally need to : remove voxels not in the target space mask,
% re-set empty voxels removed.  All this is handled in resample_space.
% Define voxel-space mapping between two images.
%
% To resample:
% 3D vox values = interp3(y_vox_orig, x_vox_orig, z_vox_orig, voldata, y_vox_query, x_vox_query, z_vox_query)
%
% :Inputs:
%
%   **obj**
%        Source image object (image_vector, including fmri_data, atlas, statistic_image) with fields:
%            obj.dat
%            obj.volInfo
%
%   **target_obj**
%        Target image object defining the desired output space.
%
% :Outputs:
%
%   **space_mapper_struct**
%        Structure containing voxel-space mapping variables:
%
%        **T**
%            4×4 affine transform mapping target voxels → source voxels
%
%        **dim_from**
%            Dimensions of source image [Nx Ny Nz]
%
%        **dim_to**
%            Dimensions of target image [Nx Ny Nz]
%
%        **x_vox_orig, y_vox_orig, z_vox_orig**
%            Source voxel grid axes
%
%        **x_vox_query, y_vox_query, z_vox_query**
%            Source voxel coordinates used as interpolation query points
%
% :Examples:
% ::
%
%    obj = load_image_set('emotionreg');
%    template_obj = fmri_data(which('brainmask.nii'));
%    space_mapper_struct = define_space_mapping(obj, template_obj);
%
% :See also:
%    resample_space
%    reconstruct_image
%    interpn, interp3

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

% NOTE:
% SPM does
% Mslice  = spm_matrix([0 0 j]);      % Matrix specifying this slice
% Mtrans  = Mto \ Mmap \ Mslice;          % Affine mappping mtx: Mask -> TOvol
% Mtrans = same as inv(Mto) * inv(Mmap) * Mslice
%
% v_query = inv(Mmap) * Mto * Mslice * target_pixels
%
% T = M_{source}^{-1} M_{target}
% T = A^{-1} * B --> maps space of A to space of B
% T = A \ B; --> maps space of B to space of A
% source_voxel = inv(Mfrom) * Mto * target_voxel
% T = Vfrom.mat \ Vto.mat;
% T = Msource \ Mtarget;
% v_query = T * v_target

% T = Vfrom.mat \ Vto.mat;
% ijk_target = [i j k 1]';
% ijk_query = T * ijk_target;

% -------------------------------------------------------------------------
% Get volume info
% -------------------------------------------------------------------------

Vfrom = obj.volInfo;
Vto   = target_obj.volInfo;

dim_from = Vfrom.dim(1:3);
dim_to   = Vto.dim(1:3);

% -------------------------------------------------------------------------
% Affine transform from target voxel space -> source voxel space
% -------------------------------------------------------------------------

T = Vfrom.mat \ Vto.mat;     % T = inv(Msource) * Mtarget

% -------------------------------------------------------------------------
% Create target voxel grid in MATLAB voxel coordinates
% -------------------------------------------------------------------------

[Xt, Yt, Zt] = ndgrid(1:dim_to(1), 1:dim_to(2), 1:dim_to(3));

% Convert target voxel centers to NIfTI/SPM voxel-center coordinates
% ijk_target = [Xt(:) - 1, Yt(:) - 1, Zt(:) - 1]';

% NOTE: don't do the above
% No extra origin offset is needed if you stay in SPM voxel 
% coordinates throughout. SPM’s own documentation and old manual 
% text describe the first voxel as (1,1,1), and SPM source code 
% also contains explicit conversions only when exporting to formats 
% that want a (0,0,0) origin.

ijk_target = [Xt(:)'; Yt(:)'; Zt(:)'; ones(1,numel(Xt))];

% save for output - don't need to save now
% x_vox_target = ijk_target(1, :); % + 1;
% y_vox_target = ijk_target(2, :); % + 1;
% z_vox_target = ijk_target(3, :); % + 1;

% -------------------------------------------------------------------------
% Transform target voxel coordinates into source voxel query coordinates
% -------------------------------------------------------------------------

ijk_query = T * ijk_target;

% ijk_query(1:3,:) = ijk_query(1:3,:) + 1;

x_vox_query = ijk_query(1,:);
y_vox_query = ijk_query(2,:);
z_vox_query = ijk_query(3,:);

% -------------------------------------------------------------------------
% Define source voxel grid vectors
% -------------------------------------------------------------------------

x_vox_orig = 1:dim_from(1);
y_vox_orig = 1:dim_from(2);
z_vox_orig = 1:dim_from(3);


% The code below defines an alternate, less efficient way of accomplishing
% the same thing, by transforming into XYZmm coordinates. The
% transformation matrix T makes it unnecessary to transform to XYZmm
% directly.

% -------------------------------------------------------------------------
% Convert target voxel centers to mm coordinates
% -------------------------------------------------------------------------
% 
% XYZmm = Vto.mat * [ijk_target; ones(1, size(ijk_target, 2))];
% XYZmm = XYZmm(1:3, :);
%
% -------------------------------------------------------------------------
% Convert target mm coordinates into source voxel coordinates
% (Alternate way, redundant)
% -------------------------------------------------------------------------
% Convert source-space query points back to MATLAB voxel coordinates
% Note: Don't add 1 because we are in SPM affine space not nifti. first vox
% is 1, 1, 1
% 
% ijk = Vfrom.mat \ XYZmm_homogeneous;
% ijk = ijk_query(1:3, :);
% 
% x_vox_target = ijk_from(1, :); % + 1;
% y_vox_target = ijk_from(2, :); % + 1;
% z_vox_target = ijk_from(3, :); % + 1;

% -------------------------------------------------------------------------
% Collect diagnostic / intermediate variables
% -------------------------------------------------------------------------

space_mapper_struct = struct();
space_mapper_struct.Vfrom = Vfrom;
space_mapper_struct.Vto = Vto;
space_mapper_struct.dim_from = dim_from;
space_mapper_struct.dim_to = dim_to;

space_mapper_struct.T = T;

% space_mapper_struct.Xt = Xt;
% space_mapper_struct.Yt = Yt;
% space_mapper_struct.Zt = Zt;

% space_mapper_struct.ijk_target = ijk_target;
% space_mapper_struct.XYZmm_homogeneous = XYZmm_homogeneous;
% space_mapper_struct.XYZmm = XYZmm;
% space_mapper_struct.ijk_from_homogeneous = ijk_from_homogeneous;
% space_mapper_struct.ijk_from = ijk_from;

space_mapper_struct.x_vox_orig = x_vox_orig;
space_mapper_struct.y_vox_orig = y_vox_orig;
space_mapper_struct.z_vox_orig = z_vox_orig;

space_mapper_struct.x_vox_query = x_vox_query;
space_mapper_struct.y_vox_query = y_vox_query;
space_mapper_struct.z_vox_query = z_vox_query;

space_mapper_struct.interp3_grid_order = {'y_vox_orig', 'x_vox_orig', 'z_vox_orig'};
space_mapper_struct.interp3_query_order = {'y_vox_query', 'x_vox_query', 'z_vox_query'};

space_mapper_struct.interp_instructions = 'To resample voldata: 3D vox values = interp3(y_vox_orig, x_vox_orig, z_vox_orig, voldata, y_vox_query, x_vox_query, z_vox_query)';

end % define_space_mapping

