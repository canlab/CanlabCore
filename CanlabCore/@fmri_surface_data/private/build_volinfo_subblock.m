function [volInfo, rows] = build_volinfo_subblock(bm)
% Build an SPM-style volInfo describing ONLY the subcortical/volumetric voxel
% models of a CIFTI brain_model, plus the grayordinate row indices they occupy.
%
% Used by fmri_surface_data to (a) populate the inherited .volInfo slot for the
% volume sub-block (so volInfo-aware code and to_fmri_data work) and (b) export
% the subcortex to an fmri_data object. Returns empty if there is no volume part.
%
% :Inputs:
%   **bm:** a brain_model struct (fmri_surface_data.brain_model / CIFTI diminfo{1}),
%           with .models{i} (.type 'surf'|'vox', .start, .count, .voxlist 3xN 0-based)
%           and .vol (.dims, .sform).
%
% :Outputs:
%   **volInfo:** struct (.mat 1-based SPM affine, .dim, .dt, .xyzlist Nx3 1-based,
%               .nvox, .image_indx, .wh_inmask, .n_inmask, .fname). [] if no volume.
%   **rows:**    column vector of row indices into the full grayordinate .dat that
%               correspond, in order, to volInfo.xyzlist / volInfo.wh_inmask.
%
% :See also: fmri_surface_data, to_fmri_data, reconstruct_image, extract_vol_from_cifti

volInfo = [];
rows = [];

if isempty(bm) || ~isfield(bm, 'vol') || isempty(bm.vol)
    return
end
voxmodels = find(cellfun(@(m) strcmp(m.type, 'vox'), bm.models));
if isempty(voxmodels)
    return
end

dims = bm.vol.dims(:)';

% CIFTI sform maps 0-based IJK -> mm. SPM .mat maps 1-based voxel -> mm.
% Convert: mat(:,4) = sform(:,4) - sform(:,1:3)*[1;1;1].
affine = bm.vol.sform;
affine(:, 4) = affine(:, 4) - sum(affine(:, 1:3), 2);

allvox = zeros(0, 3);    % Nx3, 1-based IJK, in grayordinate (voxlist) order
rows   = zeros(0, 1);
for k = voxmodels(:)'
    m = bm.models{k};
    vx = m.voxlist;                 % 3xN, 0-based
    allvox = [allvox; vx' + 1];     %#ok<AGROW> Nx3, 1-based
    rows   = [rows; (m.start:(m.start + m.count - 1))']; %#ok<AGROW>
end

ind = sub2ind(dims, allvox(:, 1), allvox(:, 2), allvox(:, 3));   % linear, voxlist order

volInfo = struct();
volInfo.mat = affine;
volInfo.dim = dims;
volInfo.dt = [16 0];
volInfo.xyzlist = allvox;
volInfo.nvox = prod(dims);
volInfo.image_indx = false(prod(dims), 1);
volInfo.image_indx(ind) = true;
volInfo.wh_inmask = ind;            % aligned with xyzlist / rows (NOT sorted)
volInfo.n_inmask = size(allvox, 1);
volInfo.fname = '';
end
