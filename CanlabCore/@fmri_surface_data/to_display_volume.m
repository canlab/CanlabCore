function vol = to_display_volume(obj)
% to_display_volume Project a surface/grayordinate object to an MNI volume for
% rendering on an ARBITRARY (non-standard) MNI isosurface.
%
% :Usage:
% ::
%     vol = to_display_volume(surf_obj)
%
% When surface data must be shown on a mesh that is not one of its own standard
% surfaces (e.g. an addbrain 'hires left' pial surface, a 'cutaway', or any MNI
% isosurface), the data is first turned into a volume so a patch can be coloured
% by sampling that volume at each vertex's MNI coordinate (see
% image_vector.render_on_surface). This is the projection step:
%
%   - Cortex: resampled to fsaverage-164k (resample_surface) if not already, then
%     surf2vol (CBIG registration fusion) to MNI 2 mm.
%   - Subcortex: to_fmri_data (its grayordinate voxels are already volumetric).
%   - Mixed grayordinate objects: the two are merged (subcortex overlaid on the
%     cortical volume; the regions are spatially disjoint).
%
% Label/binary data resamples by nearest neighbour automatically (resample_surface
% detects .dlabel / binary), preserving discrete values.
%
% :Inputs:
%   **obj:** an fmri_surface_data object (any surface_space, with cortex and/or
%            subcortical voxel models).
%
% :Outputs:
%   **vol:** an fmri_data in MNI space (2 mm) covering the projected cortex and/or
%            subcortex, ready for image_vector.render_on_surface.
%
% :See also: render_on_surface, surf2vol, resample_surface, to_fmri_data

has_cortex = ~isempty(obj.brain_model) && ...
    any(cellfun(@(m) strcmp(m.type, 'surf'), obj.brain_model.models));
has_vox = ~isempty(obj.brain_model) && ...
    any(cellfun(@(m) strcmp(m.type, 'vox'), obj.brain_model.models));

volc = []; vols = [];

if has_cortex
    fs = obj;
    if ~strcmp(obj.surface_space, 'fsaverage_164k')
        % surf2vol needs fsaverage (the CBIG warp is fsaverage-based); get there
        % natively. resample_surface auto-selects nearest for label/binary data.
        fs = resample_surface(obj, 'fsaverage_164k');
    end
    volc = surf2vol(fs);
end

if has_vox
    vols = to_fmri_data(obj);
end

if isempty(volc) && isempty(vols)
    error('fmri_surface_data:to_display_volume:empty', ...
        'Object has no cortex or subcortex to project to a volume.');
elseif isempty(vols)
    vol = volc;
elseif isempty(volc)
    vol = vols;
else
    vol = local_merge(volc, vols);
end
end


% -------------------------------------------------------------------------
function vol = local_merge(volc, vols)
% Overlay the (spatially disjoint) subcortex onto the cortical volume. Falls back
% to the cortical volume if the subcortex cannot be placed on the cortical grid.
try
    vols = resample_space(vols, volc);
    volc = replace_empty(volc);
    vols = replace_empty(vols);
    d = volc.dat; s = vols.dat;
    m = any(s ~= 0 & isfinite(s), 2);      % voxels the subcortex actually fills
    d(m, :) = s(m, :);
    volc.dat = d;
    vol = remove_empty(volc);
catch
    vol = volc;
end
end
