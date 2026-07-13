function vol = obj_to_volume(obj)
% obj_to_volume Convert a surface/grayordinate object to a volumetric fmri_data.
%
% Used by rendering when the target is an arbitrary MNI surface (e.g. an addbrain
% pial surface) rather than the object's own mesh: the data is first projected to
% a volume so the standard image_vector.render_on_surface (which samples a volume
% at patch vertices) can be reused.
%
% Delegates to the public to_display_volume, which handles any surface space
% (fs_LR is resampled to fsaverage first) and merges cortex + subcortex.
%
% :See also: to_display_volume, surf2vol, to_fmri_data, surface, render_on_surface

vol = to_display_volume(obj);
end
