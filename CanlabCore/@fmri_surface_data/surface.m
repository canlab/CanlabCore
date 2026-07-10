function han = surface(obj, varargin)
% surface Render grayordinate/surface data on cortical surfaces.
%
% :Usage:
% ::
%     han = surface(obj)                              % native inflated, 4 views
%     han = surface(obj, 'surftype', 'midthickness')  % native, other mesh
%     han = surface(obj, 'which_image', 2, 'clim', [-3 3])
%     han = surface(obj, 'existingsurface', patch_handles)   % render on given patches
%     han = surface(obj, 'mni_surface', 'left')       % render on an addbrain MNI surface
%
% Renders an fmri_surface_data object on cortical surfaces, mirroring
% image_vector.surface. Three modes:
%
%   1. NATIVE (default): builds a managed, stateful fmridisplay whose surface
%      views are the mesh set that matches the object's surface_space (fs_LR-32k
%      -> 'foursurfaces_hcp', fsaverage-164k -> 'foursurfaces_freesurfer'), then
%      paints the data as a MANAGED surface-native layer (colored DIRECTLY from
%      the per-vertex data -- no resampling; medial wall and zeros render gray).
%      Because the returned object is a stateful fmridisplay under a controller,
%      set_colormap / set_opacity / rethreshold / removeblobs / refresh act on
%      the surfaces (this is why the colormap is now changeable after the fact).
%
%   2. EXISTING SURFACE ('existingsurface', han): colors patch handles you
%      already have (e.g. from addbrain or a prior surface call). Matching-space
%      meshes are colored directly; arbitrary MNI surfaces are handled by
%      projecting to a volume (see render_on_surface).
%
%   3. MNI SURFACE ('mni_surface', name): creates an addbrain surface (e.g.
%      'left', 'right', 'hcp inflated') and renders onto it, projecting to a
%      volume when the surface is not the object's native mesh.
%
% :Optional Inputs:
%   **'surftype':**     'inflated' (default), 'midthickness', 'sphere' (fs_LR).
%   **'which_image':**  map (column) to render. Default 1.
%   **'existingsurface':** vector of patch handles to color.
%   **'mni_surface':**  an addbrain keyword (string).
%
%   Color options (harmonized with the volume pipeline; forwarded to
%   render_on_surface -- see there for details):
%   **'clim' / 'cmaprange':** [lo hi] color limits (default symmetric from data).
%   **'colormap' / 'colormapname':** named MATLAB colormap or [n x 3] matrix
%        (a single sequential map over clim).
%   **'pos_colormap' / 'neg_colormap':** [n x 3] split colormaps (default hot/cool).
%   **'splitcolor':**   {neg_low, neg_high, pos_low, pos_high} colors.
%   **'maxcolor' / 'mincolor':** endpoint colors -> a single gradient over clim.
%   **'color':**        a single solid color for all in-data vertices.
%
% :Outputs:
%   **han:** In native mode (default) a stateful **fmridisplay** object with the
%        surfaces registered as managed views and the data added as a layer; use
%        set_colormap / removeblobs / refresh on it. In the 'existingsurface' and
%        'mni_surface' modes a struct with fields .figure, .axes, .surfaces.
%
% :Examples:
% ::
%     s = fmri_surface_data(which('S1200.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii'));
%     surface(s, 'which_image', 1);
%
%     v = ttest(load_image_set('emotionreg'));   % a volumetric statistic_image
%     surface(vol2surf(v));                      % native fsaverage render
%     surface(vol2surf(v), 'mni_surface', 'left'); % on an addbrain MNI surface
%
% :See also: render_on_surface, load_surface_geom, addbrain, fmri_surface_data

surftype = 'inflated';
which_image = 1;
existing = [];
mni_surface = '';
coloropts = {};             % forwarded to render_on_surface (harmonized colors)

i = 1;
while i <= numel(varargin)
    switch lower(char(varargin{i}))
        case 'surftype',        surftype = varargin{i+1};    i = i + 2;
        case 'which_image',     which_image = varargin{i+1}; i = i + 2;
        case {'existingsurface','surface_handles'}, existing = varargin{i+1}; i = i + 2;
        case 'mni_surface',     mni_surface = varargin{i+1}; i = i + 2;
        otherwise
            % Any other option (clim, colormap, cmaprange, pos_colormap /
            % neg_colormap, splitcolor, maxcolor / mincolor, color, ...) is
            % forwarded to render_on_surface, which harmonizes them with the
            % volume visualization color pipeline.
            coloropts = [coloropts, varargin(i:i+1)]; %#ok<AGROW>
            i = i + 2;
    end
end

ropts = [{'which_image', which_image}, coloropts];

% ---- Mode 2: existing handles ----
if ~isempty(existing)
    render_on_surface(obj, existing, ropts{:});
    han = struct('figure', ancestor(existing(1), 'figure'), 'axes', [], 'surfaces', existing);
    return
end

% ---- Mode 3: addbrain MNI surface ----
if ~isempty(mni_surface)
    hp = addbrain(mni_surface);
    render_on_surface(obj, hp, ropts{:});
    han = struct('figure', gcf, 'axes', gca, 'surfaces', hp);
    return
end

% ---- Mode 1: native managed render (returns a stateful fmridisplay) ----
% Build a managed display whose surface views are the mesh set that matches this
% object's space, then add the data as a managed surface-native layer. This is
% what makes the colormap changeable afterwards (set_colormap / removeblobs act
% on the returned object). The heavy lifting -- pick the matching surfaces, add
% them, paint at full fidelity -- is done by @fmridisplay/surface's
% fmri_surface_data path, so there is a single code path for surface(obj) and
% surface(o2, obj).
o2 = fmridisplay;
o2 = surface(o2, obj, ropts{:});

% Honor a non-default surftype (midthickness / sphere) by swapping the managed
% patches' geometry to the requested mesh. The foursurfaces keyword draws
% inflated meshes; vertex ordering is shared across fs_LR surftypes, so swapping
% both vertices AND faces (self-consistent, from the same file) keeps the
% painted per-vertex data aligned. Best-effort: leave inflated on any failure.
if ~strcmpi(surftype, 'inflated')
    try
        o2 = swap_managed_geom(o2, obj.surface_space, surftype);
    catch err
        warning('fmri_surface_data:surface:surftype', ...
            'Could not apply surftype ''%s'' (%s); showing inflated.', surftype, err.message);
    end
end

han = o2;
end


function o2 = swap_managed_geom(o2, surface_space, surftype)
% Replace the Vertices/Faces of the managed cortical patches with the requested
% surftype's mesh (same space, same vertex count -> painted data stays aligned).
geom = load_surface_geom(surface_space, surftype);
for i = 1:numel(o2.surface)
    h = o2.surface{i}.object_handle;
    for hh = h(ishandle(h))'
        if ~strcmp(get(hh, 'Type'), 'patch'), continue; end
        nv  = size(get(hh, 'Vertices'), 1);
        tag = lower(get(hh, 'Tag'));
        if nv == size(geom.vertices_rh, 1) && contains(tag, 'right')
            set(hh, 'Vertices', geom.vertices_rh, 'Faces', geom.faces_rh);
        elseif nv == size(geom.vertices_lh, 1)
            set(hh, 'Vertices', geom.vertices_lh, 'Faces', geom.faces_lh);
        end
    end
end
end
