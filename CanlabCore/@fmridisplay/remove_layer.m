function obj = remove_layer(obj, k)
% Remove a single blob layer from an fmridisplay object.
%
% Deletes layer k's montage blob graphics and colorbar legend, drops the layer
% from activation_maps, re-renders the remaining layers onto the surfaces, and
% refreshes an open controller. (removeblobs removes ALL layers; this removes
% just one — e.g. the controller's per-layer "Remove layer" button.)
%
% :Usage:
% ::
%
%     o2 = remove_layer(o2, 2)   % remove the 2nd blob layer
%
% :Inputs:
%
%   **obj:** an fmridisplay object (handle).
%   **k:**   index of the layer to remove (into obj.activation_maps).
%
% :Outputs:
%
%   **obj:** the same handle, with layer k removed.
%
% :See also:
%   - addblobs, removeblobs, controller
%
% ..
%    2026 visualization overhaul
% ..

if k < 1 || k > numel(obj.activation_maps), return, end

L = obj.activation_maps{k};

% Delete this layer's montage blob graphics and any colorbar legend
if isfield(L, 'blobhandles') && ~isempty(L.blobhandles)
    bh = L.blobhandles; delete(bh(ishandle(bh)));
end
if isfield(L, 'legendhandle') && ~isempty(L.legendhandle)
    lh = L.legendhandle; delete(lh(ishandle(lh)));
end

obj.activation_maps(k) = [];

% Surfaces show an RGB composite of the layers; recomposite from the remaining
% layers (resets to gray, repaints each), so the removed layer's colour is cleared.
if ~isempty(obj.surface)
    obj = composite_surfaces(obj);
end

% Layer count changed -> rebuild an open controller
obj = update_controller(obj);

end
