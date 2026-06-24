function obj = refresh(obj, varargin)
% Re-render stored blob layers in place on their montages, from retained source.
%
% This is the engine behind the stateful visualization overhaul: because
% fmridisplay is now a handle class and each blob layer retains its source
% data and render options (see addblobs / VISUALIZATION_OVERHAUL_NOTES.md),
% a layer can be redrawn without the caller re-supplying anything. set_colormap,
% set_opacity, and rethreshold all mutate stored layer state and then call
% refresh to propagate the change to the figure(s).
%
% :Usage:
% ::
%
%     obj = refresh(obj)              % re-render all blob layers
%     obj = refresh(obj, wh_layers)   % re-render only the listed layer indices
%
% :Inputs:
%
%   **obj:**
%        An fmridisplay object (handle) with one or more activation_maps.
%
%   **wh_layers:** (optional)
%        Numeric vector of layer indices (into obj.activation_maps) to
%        refresh. Default: all layers.
%
% :Outputs:
%
%   **obj:**
%        The same handle, with each refreshed layer's blob graphics deleted
%        and redrawn on its registered montages. cmaprange is updated.
%
% :Note:
%   refresh operates on MONTAGE views. Surface re-rendering is a documented
%   follow-up; surface blobs are left untouched by refresh.
%
% :See also:
%   - addblobs, removeblobs, rethreshold, set_colormap, set_opacity
%
% ..
%    2026 visualization overhaul
% ..

wh_layers = 1:numel(obj.activation_maps);
if ~isempty(varargin) && isnumeric(varargin{1}) && ~isempty(varargin{1})
    wh_layers = varargin{1};
end

for k = wh_layers

    if k < 1 || k > numel(obj.activation_maps), continue, end

    layer = obj.activation_maps{k};

    if ~isfield(layer, 'render_args') || isempty(layer.render_args)
        % Layer predates source retention (or has no stored options); it
        % cannot be re-rendered from source. Leave it as-is.
        warning('fmridisplay:refresh', ...
            'Layer %d has no stored render options; cannot refresh. Re-add via addblobs.', k);
        continue
    end

    % Delete this layer's existing montage blob handles
    if isfield(layer, 'blobhandles') && ~isempty(layer.blobhandles)
        wh = ishandle(layer.blobhandles);
        if any(wh), delete(layer.blobhandles(wh)); end
    end
    obj.activation_maps{k}.blobhandles = [];

    % Target montages: those the layer was originally drawn on (fall back to all)
    wh_montage = layer.wh_montage;
    if isempty(wh_montage), wh_montage = 1:numel(obj.montage); end

    currentmap = obj.activation_maps{k};

    for i = wh_montage

        if numel(obj.montage) < i, continue, end

        [blobhan, cmaprange] = render_blobs(currentmap, obj.montage{i}, obj.SPACE, layer.render_args{:});

        if ~isempty(blobhan)
            obj.activation_maps{k}.blobhandles = [obj.activation_maps{k}.blobhandles; blobhan];
            obj.activation_maps{k}.cmaprange = cmaprange;
        end
    end

end

if ~iscell(obj.history), obj.history = {}; end
obj.history{end + 1} = 'refresh: re-rendered blob layers from retained source';

end
