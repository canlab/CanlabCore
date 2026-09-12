function obj = remove_legend(obj)
% Remove colorbar legends (montage and surface) from an fmridisplay object.
%
% Deletes the legend axes tracked per blob layer in
% activation_maps{}.legendhandle (drawn by legend() on the montage figure --
% continuous ramps and the discrete indexmap legends alike -- and the surface
% colorbars from render_on_surface), any legend axes tagged
% 'fmridisp_fig_legend' on the object's montage figures, and any stray ColorBar
% objects in the object's registered surface figures. The blobs themselves are
% kept; only the legends are removed. Useful when a surface colorbar overlaps
% the surface or you simply don't want it. This is what the display controller's
% 'Toggle legend' button calls to turn figure legends off.
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
    % per-view tracking used by legend() / render_layer_surfaces to replace only
    % their own legend on redraw
    for f = {'montage_legendhandle', 'surface_legendhandle'}
        if isfield(obj.activation_maps{k}, f{1}) && ~isempty(obj.activation_maps{k}.(f{1}))
            lh = obj.activation_maps{k}.(f{1});
            delete(lh(ishandle(lh)));
            obj.activation_maps{k}.(f{1}) = [];
        end
    end
end

% Sweep legend axes tagged by the controller's Toggle legend on the montage figures
for i = 1:numel(obj.montage)
    ah = obj.montage{i}.axis_handles;
    ah = ah(ishandle(ah));
    if isempty(ah), continue, end
    fig = ancestor(ah(1), 'figure');
    if isempty(fig) || ~isvalid(fig), continue, end
    delete(findobj(fig, 'Tag', 'fmridisp_fig_legend'));
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
