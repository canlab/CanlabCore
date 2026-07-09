function obj = add_surface_blobs(obj, surf_obj, varargin)
% add_surface_blobs Add an fmri_surface_data as a SURFACE-NATIVE blob layer.
%
% :Usage:
% ::
%     o2 = add_surface_blobs(o2, surf_obj, [color options])
%     o2 = addblobs(o2, surf_obj, ...)      % addblobs dispatches here
%
% Registers a grayordinate/surface object as a managed layer on an fmridisplay,
% painted DIRECTLY on any registered cortical mesh whose vertex count matches the
% object's space (e.g. the fs_LR-32k 'hcp inflated' surfaces for a native CIFTI
% object, or the fsaverage-164k surfaces for a vol2surf result). No volume
% resampling -- the per-vertex values color the vertices at full fidelity, using
% the same central canlab_colormap value->color map as montages, so colors match.
%
% The layer participates in the stateful display: set_colormap, set_opacity,
% refresh, removeblobs, and the controller act on it like a volume layer. It has
% no volume representation, so it does not appear on slice montages (use
% to_fmri_data / surf2vol for the volumetric part), and it is skipped on any
% registered surface whose mesh does not match the object's space.
%
% :Inputs:
%   **obj:**      an fmridisplay object with one or more surface views (add them
%                 with surface(obj, 'hcp inflated') etc.).
%   **surf_obj:** an fmri_surface_data object.
%
% :Optional Inputs:
%   **'which_image':** map (column) to display. Default 1.
%   **'cmaprange':**   [lo hi] (or 4-element for splitcolor) color range. Default:
%                      robust per-arm range from the data (canlab_default_cmaprange).
%   Any color option understood by canlab_colormap.from_render_args:
%   'colormap', 'pos_colormap'/'neg_colormap', 'splitcolor', 'maxcolor'/'mincolor',
%   'color'. Opacity via 'transvalue' (set_opacity / controller).
%
% :Examples:
% ::
%     o2 = fmridisplay;
%     o2 = surface(o2, 'hcp inflated left'); o2 = surface(o2, 'hcp inflated right');
%     t  = ttest(cat(sub1, sub2, sub3));          % a surface statistic map
%     o2 = addblobs(o2, t, 'colormap', 'hot');    % paint it on the fs_LR surfaces
%     o2 = set_colormap(o2, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);   % recolor
%
% :See also: addblobs, render_layer_surfaces, set_colormap, fmri_surface_data.surface

% ---- which_image (a single layer shows one map) ----
which_image = 1;
wh = find(strcmp(varargin, 'which_image'), 1);
if ~isempty(wh), which_image = varargin{wh + 1}; varargin(wh:wh + 1) = []; end

% ---- color range: explicit, else robust default from the data ----
cmaprange = [];
wh = find(strcmp(varargin, 'cmaprange'), 1);
if ~isempty(wh), cmaprange = varargin{wh + 1}; end
if isempty(cmaprange)
    v = double(surf_obj.dat(:, which_image));
    v = v(v ~= 0 & isfinite(v));
    if ~isempty(v)
        split_flag = {};
        if any(strcmp(varargin, 'splitcolor')), split_flag = {'splitcolor'}; end
        cmaprange = canlab_default_cmaprange(v, split_flag{:});
    end
end

% ---- which surfaces (default all) ----
wh_surface = 1:numel(obj.surface);
whs = find(strcmp(varargin, 'wh_surface') | strcmp(varargin, 'wh_surfaces'), 1);
if ~isempty(whs), wh_surface = varargin{whs + 1}; end

% ---- build a layer struct compatible with the volume layers ----
layer = struct('mapdata', [], 'V', [], 'SPACE', [], 'blobhandles', [], ...
    'cmaprange', cmaprange, 'mincolor', [0 0 1], 'maxcolor', [1 0 0], 'color', [], ...
    'source_region', [], 'wh_montage', [], 'wh_surface', wh_surface, ...
    'applied_threshold', [], 'which_image', which_image, 'visible', true, ...
    'legendhandle', []);
layer.source_object  = surf_obj;
layer.source_surface = surf_obj;      % marks this as a surface-native layer
layer.render_args    = varargin;

obj.activation_maps{end + 1} = layer;
k = numel(obj.activation_maps);

if isempty(obj.surface)
    warning('fmridisplay:add_surface_blobs:nosurface', ...
        ['No surface views on this fmridisplay. Add one first, e.g. ' ...
         'o2 = surface(o2, ''hcp inflated'');']);
    return
end

obj = render_layer_surfaces(obj, k, wh_surface);
end
