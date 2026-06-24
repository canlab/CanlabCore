function obj = remove_legend(obj)
% Remove surface colorbar legends from an fmridisplay object.
%
% Deletes the colorbar/legend axes created for surface blobs (stored per blob
% layer in activation_maps{}.legendhandle) and any stray ColorBar objects in
% the object's registered surface figures. The blobs themselves are kept; only
% the colorbar legends are removed. Useful when a surface colorbar overlaps the
% surface or you simply don't want it.
%
% :Usage:
% ::
%
%     o2 = remove_legend(o2)
%
% :Inputs:
%
%   **obj:** an fmridisplay object (handle).
%
% :Outputs:
%
%   **obj:** the same handle, with surface colorbar legends deleted.
%
% :See also:
%   - addblobs, removeblobs, render_layer_surfaces, render_on_surface
%
% ..
%    2026 visualization overhaul
% ..

% Delete per-layer legend handles
for k = 1:numel(obj.activation_maps)
    if isfield(obj.activation_maps{k}, 'legendhandle') && ~isempty(obj.activation_maps{k}.legendhandle)
        lh = obj.activation_maps{k}.legendhandle;
        delete(lh(ishandle(lh)));
        obj.activation_maps{k}.legendhandle = [];
    end
end

% Sweep any remaining ColorBar objects from the registered surface figures
for i = 1:numel(obj.surface)
    h = obj.surface{i}.object_handle;
    h = h(ishandle(h));
    if isempty(h), continue, end
    fig = ancestor(h(1), 'figure');
    if isempty(fig) || ~isvalid(fig), continue, end
    cbs = findobj(fig, 'Type', 'colorbar');
    delete(cbs);
end

end
