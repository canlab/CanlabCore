function obj = render_layer_surfaces(obj, k, wh_surface, show_legend)
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

% A layer toggled invisible is skipped on surfaces (lower layers show through).
if isfield(layer, 'visible') && ~isempty(layer.visible) && ~layer.visible
    return
end

if nargin < 3 || isempty(wh_surface)
    wh_surface = 1:numel(obj.surface);
end

% Surface colorbar legends are OFF by default (the controller shows a per-layer
% legend in its own panel). The 'Toggle legend' button re-renders with
% show_legend=true to put colorbars on the surface figures for export.
if nargin < 4 || isempty(show_legend), show_legend = false; end
leg_args = {}; if ~show_legend, leg_args = {'nolegend'}; end

args = {};
if isfield(layer, 'render_args') && ~isempty(layer.render_args), args = layer.render_args; end

cmaprange = [];
if isfield(layer, 'cmaprange'), cmaprange = layer.cmaprange; end

img = region2imagevec(layer.source_region);

% Translate the layer's colour spec (in render_args) into colormaps that
% render_on_surface understands. It knows pos_colormap / neg_colormap / colormap /
% color / splitcolor but NOT maxcolor / mincolor, so we build pos/neg ramps here.
% Without this, warm/cool/winter (maxcolor/mincolor) and solid colours were
% silently ignored on surfaces and render_on_surface fell back to its default
% split hot/cool (e.g. blue negatives under a "warm" map). See render_on_surface.m:222.
wh_index = find(strcmp(args, 'indexmap'), 1);
[surf_pos_cm, surf_neg_cm] = surface_colormaps_from_args(args);
clean_args = strip_color_tokens(args);
color_args = {};
if ~isempty(surf_pos_cm), color_args = [color_args, {'pos_colormap', surf_pos_cm}]; end
if ~isempty(surf_neg_cm), color_args = [color_args, {'neg_colormap', surf_neg_cm}]; end

% Central value->RGB map (single source of truth): render_on_surface colours the
% vertices in true-colour RGB via this, so surfaces match the montage exactly and
% multi-layer compositing becomes possible. The pos/neg colormaps above are still
% passed so render_on_surface can build the colorbar legend. See canlab_colormap.
tc_map = canlab_colormap.from_render_args(args, cmaprange);

% Layer opacity (set_opacity / controller slider). On surfaces the overlay IS the
% vertex colouring, so opacity BLENDS this layer's colours with what's underneath
% (lower layers / gray) rather than making the whole patch translucent — matching
% the per-layer blending montages already do. Passed to render_on_surface as
% 'truecolor_alpha'.
wh_tv = find(strcmp(args, 'transvalue'), 1);
layer_alpha = 1;
if ~isempty(wh_tv) && isnumeric(args{wh_tv + 1}), layer_alpha = args{wh_tv + 1}; end

for i = wh_surface

    if i < 1 || i > numel(obj.surface)
        warning('Requested surface does not exist! Check input surface indices.');
        continue
    end

    surfh = obj.surface{i}.object_handle;
    surfh = surfh(ishandle(surfh));        % skip handles whose figure was closed
    if isempty(surfh), continue, end

    % NOTE: this paints layer k onto the CURRENT surface colours (no erase), so a
    % new layer composites on top of lower layers (true-colour, top wins per
    % vertex). The anatomy gray is saved once (render_on_surface) so removeblobs/
    % erase still restore gray. A full recomposite (reset to gray, repaint all
    % layers in order) is done by composite_surfaces, used by refresh / remove_layer.

    if ~isempty(wh_index)
        % Index/atlas colormap: pass straight through (render_on_surface handles it).
        idxmap = args{wh_index + 1};
        [~, bar1axis, bar2axis] = render_on_surface(img, surfh, args{:}, 'colormap', idxmap, 'indexmap', leg_args{:});
    else
        call_args = [clean_args, color_args, {'truecolor', tc_map, 'truecolor_alpha', layer_alpha}, leg_args];
        if ~isempty(cmaprange), call_args = [call_args, {'cmaprange', cmaprange}]; end
        [~, bar1axis, bar2axis] = render_on_surface(img, surfh, call_args{:});
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


function [pos_cm, neg_cm] = surface_colormaps_from_args(args)
% Build pos/neg colormaps for render_on_surface from the layer's colour spec, so
% surfaces match montages. Returns [] for both when there is no explicit spec
% (render_on_surface then uses its default split hot/cool, matching the
% render_blobs default splitcolor).
pos_cm = []; neg_cm = [];
hask = @(key) any(strcmp(args, key));
valk = @(key) args{find(strcmp(args, key), 1) + 1};

if hask('splitcolor')
    sc = valk('splitcolor');                      % {minneg maxneg minpos maxpos}
    pos_cm = colormap_tor(sc{3}, sc{4});
    neg_cm = colormap_tor(sc{1}, sc{2});

elseif hask('maxcolor') || hask('mincolor')
    maxc = [1 1 0]; minc = [1 0 0];               % render_blobs-style defaults
    if hask('maxcolor'), maxc = valk('maxcolor'); end
    if hask('mincolor'), minc = valk('mincolor'); end
    ramp = colormap_tor(minc, maxc);
    pos_cm = ramp; neg_cm = ramp;                 % single ramp for both signs (matches montage)

elseif hask('color')
    c = valk('color');
    solid = repmat(c(:)', 256, 1);
    pos_cm = solid; neg_cm = solid;               % solid on both signs (no default cool negatives)

end
% 'colormap' (a raw n x 3 matrix) and the no-spec default are left to
% render_on_surface via clean_args.
end


function out = strip_color_tokens(args)
% Remove the colour-spec tokens we translate ourselves (plus cmaprange, re-added
% from the layer), so they are not double-handled by render_on_surface. Keeps
% all other options (sourcespace, targetsurface, nooutline, smooth, colormap, ...).
keys = {'splitcolor', 'color', 'maxcolor', 'mincolor', 'pos_colormap', 'neg_colormap', 'cmaprange'};
out = {};
i = 1;
while i <= numel(args)
    if ischar(args{i}) && any(strcmp(args{i}, keys))
        i = i + 2;                                % all of these are value-bearing
        continue
    end
    out{end + 1} = args{i}; %#ok<AGROW>
    i = i + 1;
end
end
