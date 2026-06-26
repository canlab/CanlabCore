function surf = vol2surf(obj, varargin)
% vol2surf Project a volumetric image to the fsaverage cortical surface.
%
% :Usage:
% ::
%     surf = vol2surf(fmri_data_obj)
%     surf = vol2surf(fmri_data_obj, 'interp', 'nearest')
%
% Projects a volumetric image_vector / fmri_data (in an MNI152 space) onto the
% fsaverage 164k cortical surface, returning an fmri_surface_data object. Uses
% the vendored CBIG Registration-Fusion (RF-ANTs) MNI152->fsaverage mapping
% (per-vertex MNI RAS coordinates) and samples the volume at those coordinates
% with interpn. Fully native -- no FreeSurfer or Connectome Workbench required.
%
% This is the native v1 mapper. It is a fixed group-template MNI152<->fsaverage
% correspondence (correct for group MNI maps; not a per-subject ribbon mapper),
% and it targets fsaverage-164k (fsaverage<->fs_LR deformation is a later step).
%
% :Inputs:
%   **obj:** an image_vector / fmri_data / statistic_image in an MNI152 volume
%            space (uses obj.volInfo.mat for the voxel->mm affine).
%
% :Optional Inputs:
%   **'interp':** 'linear' (default) or 'nearest' (use 'nearest' for label maps).
%
% :Outputs:
%   **surf:** fmri_surface_data, surface_space 'fsaverage_164k', cortex-only
%             ([2*163842 x nMaps], left then right).
%
% :Examples:
% ::
%     dat = fmri_data(which('weights_NSF_grouppred_cvpcr.img'));
%     s = vol2surf(dat);            % fsaverage surface object
%     % surface(s)                  % (rendering: M5)
%
% :See also: surf2vol, fmri_surface_data, canlab_cbig_warp_path

interp = 'linear';
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'interp', interp = varargin{i+1};
        otherwise, error('vol2surf:badopt', 'Unknown option: %s', varargin{i});
    end
end

if isempty(obj.volInfo) || ~isfield(obj.volInfo, 'mat')
    error('vol2surf:novolinfo', 'Input has no volInfo.mat (needs a volumetric space).');
end

% Reconstruct the full 3-D/4-D volume and get the voxel->mm affine
obj = replace_empty(obj);
voldata = reconstruct_image(obj);          % [X Y Z nMaps]
if ndims(voldata) == 3, voldata = voldata(:, :, :, 1); end
mat = obj.volInfo.mat;
nMaps = size(voldata, 4);

NV = 163842;                                % fsaverage vertices per hemisphere
L = load(canlab_cbig_warp_path('lh_ras')); lh_ras = L.ras;   % 3 x NV (MNI mm)
R = load(canlab_cbig_warp_path('rh_ras')); rh_ras = R.ras;

datmat = zeros(2 * NV, nMaps, 'single');
hemis = {lh_ras, rh_ras};
offsets = [0, NV];
for h = 1:2
    ras = [hemis{h}; ones(1, NV)];
    vox = mat \ ras;                        % 4 x NV, 1-based voxel coords
    rows = offsets(h) + (1:NV);
    for k = 1:nMaps
        V = voldata(:, :, :, k);
        vals = interpn(V, vox(1, :), vox(2, :), vox(3, :), interp, 0);
        datmat(rows, k) = single(vals(:));
    end
end

% Build the fsaverage cortex-only brain_model
mL = struct('struct', 'CORTEX_LEFT', 'type', 'surf', 'start', 1, 'count', NV, ...
    'numvert', NV, 'vertlist', 0:NV-1, 'voxlist', []);
mR = struct('struct', 'CORTEX_RIGHT', 'type', 'surf', 'start', NV+1, 'count', NV, ...
    'numvert', NV, 'vertlist', 0:NV-1, 'voxlist', []);
bm = struct('type', 'dense', 'length', 2*NV, 'models', {{mL, mR}}, 'vol', []);
bm.grayordinate_type = 'cortex_only';
bm.cluster = [];

surf = fmri_surface_data('dat', datmat, 'brain_model', bm, ...
    'surface_space', 'fsaverage_164k', 'intent', 'dscalar');
if ~isempty(obj.image_names), surf.image_names = obj.image_names; end
surf.history = obj.history;
surf.history{end+1} = sprintf(['vol2surf: projected %d maps to fsaverage_164k via CBIG ' ...
    'RF-ANTs MNI152 mapping (interp=%s)'], nMaps, interp);
end
