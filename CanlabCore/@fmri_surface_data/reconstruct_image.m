function out = reconstruct_image(obj)
% reconstruct_image Reconstruct surface + volume data from a grayordinate object.
%
% :Usage:
% ::
%     out = reconstruct_image(surf_obj)
%
% Reconstructs the flat [grayordinates x maps] .dat back into its natural spatial
% representations: dense per-hemisphere vertex arrays for the cortical surface
% models (with the medial wall filled with NaN), and a 3-D/4-D volume for the
% subcortical voxel models. Unlike image_vector.reconstruct_image (which assumes
% a single volumetric space), this returns a struct, because grayordinate data
% spans surface and volume spaces.
%
% No replace_empty pre-step is needed: .dat is always the full grayordinate set
% (decision D5b).
%
% :Inputs:
%   **obj:** an fmri_surface_data object.
%
% :Outputs:
%   **out:** struct with fields (present only when that model exists):
%       .<struct_name>   dense [numvert x nMaps] per surface model, lowercased
%                        (e.g. .cortex_left, .cortex_right); medial wall = NaN.
%       .volume          [X x Y x Z x nMaps] subcortical volume (NaN outside).
%       .volume_volInfo  the SPM-style volInfo for .volume (1-based affine).
%       .models          the brain_model.models list (for reference).
%
% :Examples:
% ::
%     s = fmri_surface_data(which('transcriptomic_gradients.dscalar.nii'));
%     r = reconstruct_image(s);
%     size(r.cortex_left)      % [32492 x nMaps], medial wall NaN
%     size(r.volume)           % [91 109 91 x nMaps]
%
% :See also: fmri_surface_data, to_fmri_data, surface, render_on_surface

bm = obj.brain_model;
if isempty(bm)
    error('fmri_surface_data:reconstruct_image:nobm', 'Object has no brain_model.');
end

nMaps = size(obj.dat, 2);
out = struct();
out.models = bm.models;

% ---- Surface models -> dense per-hemisphere vertex arrays (medial wall NaN) ----
for i = 1:numel(bm.models)
    m = bm.models{i};
    if ~strcmp(m.type, 'surf'), continue; end
    dense = nan(m.numvert, nMaps);
    rows = m.start:(m.start + m.count - 1);
    dense(m.vertlist + 1, :) = obj.dat(rows, :);    % vertlist is 0-based
    fld = matlab.lang.makeValidName(lower(m.struct));
    out.(fld) = dense;
end

% ---- Volume models -> 3-D/4-D volume ----
[vi, rows] = build_volinfo_subblock(bm);
if ~isempty(vi) && ~isempty(rows)
    dims = vi.dim;
    vol = nan([dims nMaps]);
    nv = prod(dims);
    for k = 1:nMaps
        v = nan(nv, 1);
        v(vi.wh_inmask) = obj.dat(rows, k);
        vol(:, :, :, k) = reshape(v, dims);
    end
    out.volume = vol;
    out.volume_volInfo = vi;
end
end
