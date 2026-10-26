function [surface_handles, colors, patch_cell] = isosurface(r, varargin)
% isosurface Create a series of surfaces in different colors, one per region.
%
% Render each element of a region object array as a 3-D isosurface using
% imageCluster, with options for a single color across regions and for
% matching colors across left/right hemispheres.
%
% :Usage:
% ::
%
%     [surface_handles, colors, patch_cell] = isosurface(r, [optional arguments])
%
% :Inputs:
%
%   **r:**
%        A region-class object array.
%
% :Optional Inputs:
%
%   **Any optional inputs to imageCluster:**
%        e.g., 'alpha' followed by transparency value, 'sd', or
%        'smoothbox'.
%
%   **'colors':**
%        Followed by a single color in { } or a cell array of multiple
%        colors.
%
%   **'nomatchleftright' or 'nosymmetric':**
%        Do not match colors across hemispheres (left/right). The
%        default matches L/R and may override user-supplied colors.
%
%   **'smooth':**
%        Build very smooth, volume-preserving surfaces instead of the
%        standard imageCluster isosurfaces. Steps:
%          1. Gaussian-smooth the binary region mask (smooth3)
%          2. Upsample the smoothed volume (cubic) and extract the
%             isosurface at the level that preserves the region's volume
%          3. Taubin (lambda/mu) mesh smoothing, which removes the voxel
%             staircase without shrinking the surface
%          4. Rescale the mesh about its centroid so its enclosed volume
%             equals the region's voxel volume
%        Related optional inputs (used only with 'smooth'):
%          'smooth_sd'   Gaussian SD in voxels. Default 1
%          'upsample'    Upsampling factor. Default 3
%          'taubin_iter' Taubin iterations. Default 50
%          'reduce'      Fraction of faces kept (reducepatch). Default .5
%        'alpha' is also applied. 'sd' and 'smoothbox' apply only to the
%        standard (imageCluster) surfaces.
%
% :Outputs:
%
%   **surface_handles:**
%        Vector of patch handles, one per region (legacy form; coerces
%        Patch objects to doubles).
%
%   **colors:**
%        Cell array of [r g b] color triplets used for each region.
%
%   **patch_cell:**
%        Cell array of patch object handles, one per region.
%
% :Examples:
% ::
%
%     atlasfile = which('Morel_thalamus_atlas_object.mat');
%     load(atlasfile)
%
%     surface_handles = isosurface(r);
%     surface_handles = isosurface(r, 'alpha', .5);
%     surface_handles = isosurface(r, 'alpha', .5, 'nomatchleftright');
%
%     view(135, 30);
%     lightRestoreSingle;
%     lightFollowView;
%
%     load(which('CIT168_MNI_subcortical_atlas_object.mat'));
%     r = atlas2region(atlas_obj);
%
%     surface_handles = isosurface(r, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]});
%     surface_handles = isosurface(r, 'alpha', .5, 'color', {[.3 .6 .4] [.5 .4 .2]}, 'nomatchleftright');
%     p = addbrain('hires right');
%     lightFollowView;
%
%     % Very smooth, volume-preserving surfaces
%     atl = load_atlas('canlab2024');
%     r = region(fmri_data(select_atlas_subset(atl, {'Ctx_a24pr_R'}, 'flatten')));
%     figure; isosurface(r, 'color', {[1 .25 .05]}, 'smooth');
%     figure; isosurface(r, 'color', {[1 .25 .05]}, 'smooth', 'smooth_sd', 1.5, 'taubin_iter', 100);
%
% :See also:
%   - atlas/isosurface
%   - imageCluster
%   - match_colors_left_right


k = length(r);

% ..
%    DEFAULTS AND INPUTS
% ..

% input_args = varargin;           % save for later
colors = scn_standard_colors(k); % generate colors
matchcolorsleftright = true;

dosmooth = false;      % 'smooth' option: volume-preserving smooth surfaces
smooth_opts = struct('sd', 1, 'upsample', 3, 'taubin_iter', 50, 'reduce', .5, 'alpha', 1);

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            % do nothing for inputs passed on to imageCluster
            case {'alpha'}
                smooth_opts.alpha = varargin{i+1};

            case 'smooth', dosmooth = true; varargin{i} = [];
            case 'smooth_sd', smooth_opts.sd = varargin{i+1}; varargin{i+1} = []; varargin{i} = [];
            case 'upsample', smooth_opts.upsample = varargin{i+1}; varargin{i+1} = []; varargin{i} = [];
            case 'taubin_iter', smooth_opts.taubin_iter = varargin{i+1}; varargin{i+1} = []; varargin{i} = [];
            case 'reduce', smooth_opts.reduce = varargin{i+1}; varargin{i+1} = []; varargin{i} = [];

            case {'color', 'colors'}
                colors = varargin{i+1}; varargin{i+1} = []; varargin{i} = [];

                if ~iscell(colors), colors = {colors}; end
                
                if length(colors) == 1 % single color; matchcolorsleftright will overwrite this
                    matchcolorsleftright = false;
                end
                
            case {'nomatchleftright', 'nosymmetric'}, matchcolorsleftright = false; varargin{i} = [];
                
            case {'sd' 'smoothbox'}
                % do nothing, but pass these into imageCluster later
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if iscolumn(colors), colors = colors'; end

% Handle single-color input, and other cases where too few colors:
while length(colors) < k, colors = [colors colors]; end

if matchcolorsleftright
    colors = match_colors_left_right(r);
end

colors = colors(1:length(r));

% -------------------------------------------------------------------------
% Make surfaces
% -------------------------------------------------------------------------

cl = region2struct(r);

surface_handles = []; % This appears to be a bug, as it coerces the Patch object into a double. 05/23/2023 MS
patch_cell = {};

for i = 1:k
    
    try
        if dosmooth
            out = smooth_region_patch(r(i), colors{i}, smooth_opts);
        else
            out = imageCluster('cluster', cl(i), 'color', colors{i}, varargin{:});
        end
        
        surface_handles(i) = out; % This appears to be a bug, as it
%         coerces the Patch object into a double. 05/23/2023 MS
        patch_cell{i} = out;
    catch
        disp('Error imaging isosurface; too few voxels?')
        surface_handles(i) = [];
        patch_cell{i}=[];
    end
    
%     % Isocaps, if needed
%     
%     % pp = [];
%     isocap = isocaps(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, V, mythresh);
%     p(end + 1) = patch(isocap, 'FaceColor', 'interp','EdgeColor', 'none', 'FaceAlpha',1);


end


% -------------------------------------------------------------------------
% Lighting, etc.
% -------------------------------------------------------------------------

lightRestoreSingle;
lighting gouraud;
axis vis3d image
material dull

end % function


% -------------------------------------------------------------------------
% Smooth, volume-preserving isosurface for one region
% -------------------------------------------------------------------------

function p = smooth_region_patch(r, color, opts)

xyz = r.XYZ;                             % 3 x n voxel coordinates
nvox = size(xyz, 2);
voxvol = abs(det(r.M(1:3, 1:3)));        % mm^3 per voxel

% Binary volume cropped to the region, padded so the surface closes
pad = ceil(3 * opts.sd) + 2;
mn = min(xyz, [], 2) - pad;
dims = (max(xyz, [], 2) + pad - mn + 1)';
vol = zeros(dims);
vol(sub2ind(dims, xyz(1, :) - mn(1) + 1, xyz(2, :) - mn(2) + 1, xyz(3, :) - mn(3) + 1)) = 1;

% 1. Gaussian smoothing
vol = smooth3(vol, 'gaussian', 2 * ceil(3 * opts.sd) + 1, opts.sd);

% 2. Upsample (meshgrid order: X = dim 2, Y = dim 1)
step = 1 / opts.upsample;
[Xq, Yq, Zq] = meshgrid(1:step:dims(2), 1:step:dims(1), 1:step:dims(3));
volq = interp3(vol, Xq, Yq, Zq, 'cubic');

% Isolevel at which the enclosed volume matches the region's volume
vals = sort(volq(:), 'descend');
level = vals(min(numel(vals), max(1, round(nvox * opts.upsample ^ 3))));
level = min(max(level, .05), .95);

fv = isosurface(Xq, Yq, Zq, volq, level);
if isempty(fv.vertices), error('Empty isosurface'); end

if opts.reduce < 1, fv = reducepatch(fv, opts.reduce); end

% 3. Taubin smoothing
fv.vertices = taubin_smooth(fv.faces, fv.vertices, opts.taubin_iter, .5, -.53);

% Voxel -> mm. isosurface vertices are [dim2 dim1 dim3]
ijk = fv.vertices(:, [2 1 3]) + mn' - 1;
mm = [ijk ones(size(ijk, 1), 1)] * r.M';
v = mm(:, 1:3);

% 4. Rescale about the centroid so the enclosed volume equals the voxel volume
f = fv.faces;
meshvol = abs(sum(dot(v(f(:, 1), :), cross(v(f(:, 2), :), v(f(:, 3), :), 2), 2)) / 6);
c = mean(v, 1);
v = c + (v - c) * (nvox * voxvol / meshvol) ^ (1/3);
fv.vertices = v;

p = patch(fv, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', opts.alpha, 'FaceLighting', 'gouraud', ...
    'AmbientStrength', .35, 'DiffuseStrength', .6, 'SpecularStrength', .25, 'SpecularExponent', 10);

end


function v = taubin_smooth(f, v, niter, lambda, mu)
% Umbrella-operator Taubin smoothing: lambda > 0 shrinks, mu < 0 re-inflates

nv = size(v, 1);
A = sparse([f(:, 1); f(:, 2); f(:, 3)], [f(:, 2); f(:, 3); f(:, 1)], 1, nv, nv);
A = double((A + A') > 0);
deg = max(sum(A, 2), 1);
W = spdiags(1 ./ deg, 0, nv, nv) * A;    % row-normalized adjacency

for k = 1:niter
    v = v + lambda * (W * v - v);
    v = v + mu * (W * v - v);
end

end


