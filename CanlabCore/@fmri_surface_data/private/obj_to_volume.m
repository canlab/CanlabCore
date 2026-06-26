function vol = obj_to_volume(obj)
% obj_to_volume Convert a surface/grayordinate object to a volumetric fmri_data.
%
% Used by rendering when the target is an arbitrary MNI surface (e.g. an addbrain
% pial surface) rather than the object's own mesh: the data is first projected to
% a volume so the standard image_vector.render_on_surface (which samples a volume
% at patch vertices) can be reused.
%
% - fsaverage_164k cortex  -> surf2vol (CBIG RF) to MNI 2 mm.
% - volume-only grayordinate -> to_fmri_data (subcortex).
% - other spaces (e.g. fsLR_32k cortex) are not yet supported here.
%
% :See also: surf2vol, to_fmri_data, surface, render_on_surface

if strcmp(obj.surface_space, 'fsaverage_164k')
    vol = surf2vol(obj);
elseif ~isempty(obj.brain_model) && any(cellfun(@(m) strcmp(m.type,'vox'), obj.brain_model.models)) ...
        && ~any(cellfun(@(m) strcmp(m.type,'surf'), obj.brain_model.models))
    vol = to_fmri_data(obj);
else
    error('fmri_surface_data:obj_to_volume:space', ...
        ['Cannot project surface_space ''%s'' to a volume for rendering on an arbitrary ' ...
         'MNI surface. Supported: fsaverage_164k (via surf2vol) or volume-only objects ' ...
         '(via to_fmri_data). Render fsLR_32k data on its native mesh instead, or pass ' ...
         'an fsaverage object.'], obj.surface_space);
end
end
