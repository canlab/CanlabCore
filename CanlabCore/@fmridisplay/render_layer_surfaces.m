function obj = render_layer_surfaces(obj, k, wh_surface)
% Render blob layer k onto registered surface view(s), in place.
%
% Shared surface-rendering engine for the stateful visualization overhaul.
% addblobs (new layer -> all surfaces), surface (pull existing layers onto a
% newly-added surface), and refresh (re-render after rethreshold / recolor)
% all route through here so montages and surfaces stay in sync. Colors are
% derived from the layer's stored render_args (the source of truth that
% set_colormap updates), so re-coloring propagates to surfaces too. See
% VISUALIZATION_OVERHAUL_NOTES.md.
%
% :Usage:
% ::
%
%     obj = render_layer_surfaces(obj, k)              % onto all surfaces
%     obj = render_layer_surfaces(obj, k, wh_surface)  % onto listed surfaces
%
% :Inputs:
%
%   **obj:**  an fmridisplay object (handle).
%   **k:**    index into obj.activation_maps (the blob layer to draw).
%
% :Optional Inputs:
%
%   **wh_surface:** numeric vector of surface indices. Default: all surfaces.
%
% :Note:
%   render_on_surface paints the surface object's vertex colors directly, so
%   when several layers target the same surface the last one rendered wins
%   (surfaces do not composite multiple layers); montages do show every layer.
%
% :See also:
%   - addblobs, surface, refresh, removeblobs, render_on_surface
%
% ..
%    2026 visualization overhaul
% ..

if isempty(obj.surface) || k < 1 || k > numel(obj.activation_maps)
    return
end

layer = obj.activation_maps{k};

if ~isfield(layer, 'source_region') || isempty(layer.source_region)
    return   % legacy layer with no retained source; nothing to render from
end

if nargin < 3 || isempty(wh_surface)
    wh_surface = 1:numel(obj.surface);
end

args = {};
if isfield(layer, 'render_args') && ~isempty(layer.render_args), args = layer.render_args; end

cmaprange = [];
if isfield(layer, 'cmaprange'), cmaprange = layer.cmaprange; end

img = region2imagevec(layer.source_region);

% Color mode, derived from render_args (kept current by set_colormap)
wh_split = find(strcmp(args, 'splitcolor'), 1);
wh_index = find(strcmp(args, 'indexmap'),   1);

% Constant opacity for the surface patch (set_opacity / controller slider).
% render_on_surface drives blob ALPHA via 'transvalue' on montages, but on a
% surface the overlay IS the vertex coloring, so we dim the patch via FaceAlpha.
wh_tv = find(strcmp(args, 'transvalue'), 1);
surf_alpha = [];
if ~isempty(wh_tv), surf_alpha = args{wh_tv + 1}; end

for i = wh_surface

    if i < 1 || i > numel(obj.surface)
        warning('Requested surface does not exist! Check input surface indices.');
        continue
    end

    surfh = obj.surface{i}.object_handle;
    surfh = surfh(ishandle(surfh));        % skip handles whose figure was closed
    if isempty(surfh), continue, end

    % If this surface already carries blobs from a prior render, restore its
    % saved gray FIRST. Otherwise render_on_surface blends onto — and re-saves
    % as its restore data — the previous blob colors, which would defeat a later
    % removeblobs (the surface could never return to anatomy gray). All patches
    % of a surface get their gray stashed together on first paint, so a non-empty
    % UserData on the first patch reliably means "already painted".
    if ~isempty(surfh) && ishandle(surfh(1)) && ~isempty(get(surfh(1), 'UserData'))
        addbrain('eraseblobs', surfh);
    end

    if ~isempty(wh_split)
        % Split +/- colormap: build pos/neg maps from the stored 4-color spec
        sc = args{wh_split + 1};                 % {minneg maxneg minpos maxpos}
        pos_colormap = colormap_tor(sc{3}, sc{4});
        neg_colormap = colormap_tor(sc{1}, sc{2});
        if ~isempty(cmaprange)
            [~, bar1axis, bar2axis] = render_on_surface(img, surfh, 'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap, 'cmaprange', cmaprange, args{:});
        else
            [~, bar1axis, bar2axis] = render_on_surface(img, surfh, 'pos_colormap', pos_colormap, 'neg_colormap', neg_colormap, args{:});
        end

    elseif ~isempty(wh_index)
        idxmap = args{wh_index + 1};
        [~, bar1axis, bar2axis] = render_on_surface(img, surfh, args{:}, 'colormap', idxmap, 'indexmap');

    else
        if ~isempty(cmaprange)
            [~, bar1axis, bar2axis] = render_on_surface(img, surfh, 'cmaprange', cmaprange, args{:});
        else
            [~, bar1axis, bar2axis] = render_on_surface(img, surfh, args{:});
        end
    end

    % Constant-opacity overlay: dim the whole surface patch to the requested
    % value so the controller's opacity slider has an effect on surfaces.
    if ~isempty(surf_alpha)
        valid = surfh(ishandle(surfh));
        if ~isempty(valid), set(valid, 'FaceAlpha', surf_alpha); end
    end

    % Track the colorbar axes as the layer's legend; drop a prior one so
    % legends from multiple surfaces don't stack up.
    if isfield(obj.activation_maps{k}, 'legendhandle') && ~isempty(obj.activation_maps{k}.legendhandle)
        old = obj.activation_maps{k}.legendhandle;
        delete(old(ishandle(old)));
    end
    obj.activation_maps{k}.legendhandle = [bar1axis, bar2axis];

end

end
