function varargout = slices(obj, varargin)
% slices Plot slices of the SUBCORTICAL volume of a grayordinate object.
%
% :Usage:
% ::
%     slices(surf_obj)
%
% A slice plot only applies to the volumetric (subcortical) grayordinates;
% cortical-surface data is rendered with surface(obj). This routes the
% subcortical grayordinates through to_fmri_data and calls image_vector/slices.
%
% :Inputs:
%   **obj:** an fmri_surface_data object with subcortical (volume) grayordinates.
%
% :See also: surface, to_fmri_data, orthviews, montage, image_vector.slices

if ~local_has_volume(obj)
    error('fmri_surface_data:slices:cortexonly', ...
        ['slices shows volumetric slices, but this object is cortex-only ' ...
         '(no subcortical grayordinates). Use surface(obj) to render the cortical surface.']);
end

vol = to_fmri_data(obj);
[varargout{1:nargout}] = slices(vol, varargin{:});
end


function tf = local_has_volume(obj)
tf = ~isempty(obj.brain_model) && ...
     any(cellfun(@(m) strcmp(m.type, 'vox'), obj.brain_model.models));
end
