function vol = surf2vol(obj, varargin)
% surf2vol Project an fsaverage surface object back to an MNI152 volume (fmri_data).
%
% :Usage:
% ::
%     vol = surf2vol(surf_obj)
%     vol = surf2vol(surf_obj, 'reference', fmri_data_obj)
%
% Projects cortical-surface data on the fsaverage 164k surface into an MNI152
% volume, returning an fmri_data object (which can then be written to .nii with
% its write method, montaged, etc.). Uses the same vendored CBIG RF-ANTs
% MNI152<->fsaverage per-vertex coordinates as vol2surf, scattering each vertex's
% value into its MNI voxel and averaging co-located vertices (accumarray). Fully
% native -- no FreeSurfer / Connectome Workbench.
%
% This is the inverse of vol2surf and is self-consistent with it. It requires an
% fsaverage_164k object (the CBIG warp is fsaverage-based). For the subcortical
% (volumetric) part of a grayordinate object, use to_fmri_data instead.
%
% :Inputs:
%   **obj:** an fmri_surface_data object in surface_space 'fsaverage_164k'
%            (e.g. produced by vol2surf), with CORTEX_LEFT / CORTEX_RIGHT models.
%
% :Optional Inputs:
%   **'reference':** an image_vector/fmri_data whose volInfo (mat + dim) defines
%                    the target grid. Default: MNI152 2 mm, [91 109 91].
%   **'interp':** 'nearest' (default) -- scatter assignment is nearest-voxel.
%
% :Outputs:
%   **vol:** fmri_data over the cortical voxels hit by the projection
%            [n_hit_voxels x nMaps], in the target MNI space.
%
% :Examples:
% ::
%     s   = vol2surf(fmri_data(which('weights_NSF_grouppred_cvpcr.img')));
%     v   = surf2vol(s);
%     % v.write('fname', '/tmp/back.nii');
%
% :See also: vol2surf, to_fmri_data, fmri_surface_data

% ---- Options ----
refobj = [];
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'reference', refobj = varargin{i+1};
        case 'interp'   % nearest only in v1; accepted for API parity
        otherwise, error('surf2vol:badopt', 'Unknown option: %s', varargin{i});
    end
end

if ~strcmp(obj.surface_space, 'fsaverage_164k')
    error('surf2vol:space', ...
        ['surf2vol requires an fsaverage_164k object (the CBIG warp is fsaverage-based). ' ...
         'This object is "%s". For the subcortical/volumetric part of a grayordinate ' ...
         'object, use to_fmri_data instead.'], obj.surface_space);
end

% ---- Target grid ----
if ~isempty(refobj) && ~isempty(refobj.volInfo) && isfield(refobj.volInfo, 'mat')
    tmat = refobj.volInfo.mat;
    dims = refobj.volInfo.dim;
else
    % Standard MNI152 2 mm grid (1-based SPM .mat), matching the CIFTI subcortical grid
    dims = [91 109 91];
    tmat = [-2 0 0 92; 0 2 0 -128; 0 0 2 -74; 0 0 0 1];
end
nvox = prod(dims);

% ---- Pull cortical hemisphere data ----
[lh_dat, rh_dat] = local_hemi_data(obj);
nMaps = size(obj.dat, 2);

NV = 163842;
L = load(canlab_cbig_warp_path('lh_ras')); lh_ras = L.ras;
R = load(canlab_cbig_warp_path('rh_ras')); rh_ras = R.ras;

acc = zeros(nvox, nMaps);
cnt = zeros(nvox, 1);
for h = 1:2
    if h == 1, ras = lh_ras; hd = lh_dat; else, ras = rh_ras; hd = rh_dat; end
    vox = round(tmat \ [ras; ones(1, NV)]);          % 1-based voxel indices
    in = vox(1,:) >= 1 & vox(1,:) <= dims(1) & ...
         vox(2,:) >= 1 & vox(2,:) <= dims(2) & ...
         vox(3,:) >= 1 & vox(3,:) <= dims(3);
    lin = sub2ind(dims, vox(1,in), vox(2,in), vox(3,in))';
    cnt = cnt + accumarray(lin, 1, [nvox 1]);
    for k = 1:nMaps
        vk = hd(in, k);
        acc(:, k) = acc(:, k) + accumarray(lin, double(vk), [nvox 1]);
    end
end

inmask = cnt > 0;
datm = acc(inmask, :) ./ cnt(inmask);

% ---- Build fmri_data ----
[ix, iy, iz] = ind2sub(dims, find(inmask));
iv = image_vector;
iv.volInfo = struct('mat', tmat, 'dim', dims, 'dt', [16 0], ...
    'xyzlist', [ix iy iz], 'nvox', nvox, ...
    'image_indx', inmask, 'wh_inmask', find(inmask), ...
    'n_inmask', nnz(inmask), 'fname', '');
iv.dat = single(datm);
iv.removed_voxels = false(size(iv.dat, 1), 1);
iv.removed_images = false(size(iv.dat, 2), 1);
iv.image_names = obj.image_names;
iv.history = obj.history;
iv.history{end+1} = sprintf('surf2vol: projected fsaverage_164k -> MNI %dx%dx%d (%d cortical voxels)', ...
    dims(1), dims(2), dims(3), nnz(inmask));

vol = fmri_data(iv);
end


% =========================================================================
function [lh_dat, rh_dat] = local_hemi_data(obj)
NV = 163842;
lh_dat = []; rh_dat = [];
for i = 1:numel(obj.brain_model.models)
    m = obj.brain_model.models{i};
    if ~strcmp(m.type, 'surf'), continue; end
    rows = m.start:(m.start + m.count - 1);
    dense = zeros(m.numvert, size(obj.dat, 2));
    dense(m.vertlist + 1, :) = double(obj.dat(rows, :));
    if strcmpi(m.struct, 'CORTEX_LEFT'),  lh_dat = dense; end
    if strcmpi(m.struct, 'CORTEX_RIGHT'), rh_dat = dense; end
end
if isempty(lh_dat) || isempty(rh_dat) || size(lh_dat,1) ~= NV
    error('surf2vol:models', 'Expected CORTEX_LEFT and CORTEX_RIGHT models with %d vertices.', NV);
end
end
