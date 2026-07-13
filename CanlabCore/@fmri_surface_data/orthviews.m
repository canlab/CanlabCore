function varargout = orthviews(obj, varargin)
% orthviews Show the SUBCORTICAL volume of a grayordinate object in SPM orthviews.
%
% :Usage:
% ::
%     orthviews(surf_obj)
%
% Cortical-surface data has no voxel representation, so a 3-plane volume viewer
% only applies to the subcortical/volumetric grayordinates. This method routes
% those through to_fmri_data and calls the standard fmri_data/orthviews, so you
% see the subcortical structures on slices. To view the cortical surface, use
% surface(obj) instead.
%
% :Inputs:
%   **obj:** an fmri_surface_data object with subcortical (volume) grayordinates.
%
% :See also: surface, to_fmri_data, montage, image_vector.orthviews

if ~local_has_volume(obj)
    error('fmri_surface_data:orthviews:cortexonly', ...
        ['orthviews shows volumetric data, but this object is cortex-only ' ...
         '(no subcortical grayordinates). Use surface(obj) to render the cortical surface.']);
end

vol = to_fmri_data(obj);
[varargout{1:nargout}] = orthviews(vol, varargin{:});
end


function tf = local_has_volume(obj)
tf = ~isempty(obj.brain_model) && ...
     any(cellfun(@(m) strcmp(m.type, 'vox'), obj.brain_model.models));
end
