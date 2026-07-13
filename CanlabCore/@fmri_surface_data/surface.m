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
%      views are the four-surface set that MATCHES the object's surface_space
%      (see "Surface spaces" below), then paints the data as a MANAGED
%      surface-native layer (colored DIRECTLY from the per-vertex data -- no
%      resampling when the mesh is the object's own space; medial wall and zeros
%      render gray). Because the returned object is a stateful fmridisplay under
%      a controller, set_colormap / set_opacity / rethreshold / removeblobs /
%      refresh act on the surfaces (this is why the colormap is changeable after
%      the fact). With NO explicit surface argument, the default view is chosen
%      to match the object's space (below); passing a surface keyword overrides it.
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
% :Surface spaces and their default display meshes:
%
%   The object's .surface_space determines the default four-surface view. Two
%   spaces have NATIVE display meshes bundled with CanlabCore and paint directly;
%   the rest render on the four-surface view of their parent / aligned space via
%   a fast, cached nearest-neighbour resample (a one-line message prints when
%   that happens). All mesh files live under
%   canlab_canonical_brains/Canonical_brains_surfaces/.
%
%   ======================================================================================
%   surface_space   verts/hemi  default view              display mesh files
%   ======================================================================================
%   fsLR_32k        32492       foursurfaces_hcp          S12000.L/R.inflated_MSMAll.32k_fsl_LR.mat   (native)
%   fsaverage_164k  163842      foursurfaces_freesurfer   surf_freesurf_inflated_Left/Right.mat        (native)
%   fsaverage6      40962       foursurfaces_freesurfer   (resampled up to fsaverage-164k, then above)
%   fsaverage5      10242       foursurfaces_freesurfer   (resampled up to fsaverage-164k, then above)
%   fsaverage4      2562        foursurfaces_freesurfer   (resampled up to fsaverage-164k, then above)
%   onavg_41k       40962       foursurfaces_hcp          (onavg is fs_LR-aligned -> resampled to fs_LR)
%   onavg_10k       10242       foursurfaces_hcp          (onavg is fs_LR-aligned -> resampled to fs_LR)
%   ======================================================================================
%
%   Other fs_LR surftypes (fsLR_32k only): 'midthickness' uses
%   S1200.L/R.midthickness_MSMAll.32k_fs_LR.surf.gii and 'sphere' uses
%   S1200.L/R.sphere.32k_fs_LR.mat (see 'surftype' below). fsaverage-164k has an
%   inflated mesh only. To render in a space with no bundled mesh at true native
%   fidelity, resample first with resample_surface (see :See also).
%
% :Optional Inputs:
%   **'surftype':**     'inflated' (default), 'midthickness', 'sphere' (fs_LR only).
%   **'which_image':**  map (column) to render. Default 1.
%   **'existingsurface':** vector of patch handles to color.
%   **'mni_surface':**  an addbrain keyword (string).
%   Any addbrain surface keyword given as a bare token (e.g. 'foursurfaces_hcp',
%   'hcp inflated left') overrides the space-matched default view.
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
%     % ---- Native four-surface render, one example per surface space ----
%     % Each object's space picks its matching four-surface view automatically.
%
%     % fs_LR-32k (HCP CIFTI grayordinates)  -> foursurfaces_hcp
%     s = fmri_surface_data(which('S1200.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii'));
%     surface(s);
%
%     % fsaverage-164k (from a volume via CBIG registration fusion) -> foursurfaces_freesurfer
%     v  = ttest(load_image_set('emotionreg'));    % a volumetric statistic_image
%     sf = vol2surf(v);                            % fsaverage-164k object
%     surface(sf, 'clim', [-4 4]);
%
%     % Nested fsaverage6 (resamples up onto the fsaverage-164k four-surface view)
%     s6 = resample_surface(s, 'fsaverage6');
%     surface(s6);
%
%     % onavg (equal-area; fs_LR-aligned, so it renders on the fs_LR view)
%     son = resample_surface(s, 'onavg_41k');
%     surface(son);
%
%     % Other fs_LR surftypes (fsLR_32k meshes only)
%     surface(s, 'surftype', 'midthickness');
%     surface(s, 'surftype', 'sphere');
%
%     % ---- Rendering onto a specific existing / MNI surface ----
%     surface(s, 'foursurfaces_hcp');              % explicit view (overrides default)
%     surface(sf, 'mni_surface', 'left');          % an addbrain MNI pial surface (via volume)
%     hp = addbrain('hcp inflated left');          % a native fs_LR patch
%     surface(s, 'existingsurface', hp);           % color it directly (no resampling)
%
% :See also: render_on_surface, resample_surface, load_surface_geom, addbrain, fmri_surface_data

surftype = 'inflated';
which_image = 1;
existing = [];
mni_surface = '';
coloropts = {};             % colour opts -> render_on_surface / the managed layer
surf_args = {};             % surface directives (which view) -> @fmridisplay/surface

% Colour options that take a VALUE (consume a pair). Everything else that is not
% a recognized mode keyword is a surface directive forwarded to the native
% managed path -- either a bare addbrain/foursurfaces keyword (single token) or a
% 'direction'/'orientation'/'axes' pair -- so surface(obj, 'foursurfaces_hcp'),
% surface(obj, 'hcp inflated left'), etc. override the space-matched default view.
color_value_keys = {'clim', 'cmaprange', 'colormap', 'colormapname', ...
    'pos_colormap', 'neg_colormap', 'splitcolor', 'maxcolor', 'mincolor', ...
    'color', 'transvalue'};
surf_pair_keys = {'direction', 'orientation', 'axes'};

i = 1;
while i <= numel(varargin)
    tok = varargin{i};
    key = ''; if ischar(tok) || isstring(tok), key = lower(char(tok)); end
    switch key
        case 'surftype',        surftype = varargin{i+1};    i = i + 2;
        case 'which_image',     which_image = varargin{i+1}; i = i + 2;
        case {'existingsurface','surface_handles'}, existing = varargin{i+1}; i = i + 2;
        case 'mni_surface',     mni_surface = varargin{i+1}; i = i + 2;
        case {'unique','solid'}                 % value-less colour-mode flags
            coloropts = [coloropts, varargin(i)]; %#ok<AGROW>
            i = i + 1;
        otherwise
            if any(strcmp(key, color_value_keys))
                coloropts = [coloropts, varargin(i:i+1)]; %#ok<AGROW>
                i = i + 2;
            elseif any(strcmp(key, surf_pair_keys))
                surf_args = [surf_args, varargin(i:i+1)]; %#ok<AGROW>
                i = i + 2;
            else
                surf_args = [surf_args, varargin(i)]; %#ok<AGROW>   % bare surface directive
                i = i + 1;
            end
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
% surface(o2, obj). A surf_args entry (e.g. 'foursurfaces_hcp', a bare addbrain
% keyword, or a direction/orientation/axes pair) overrides the space-matched
% default view; with none, @fmridisplay/surface picks the matching four-surface
% set (see surface_default_keyword).
o2 = fmridisplay;
o2 = surface(o2, obj, ropts{:}, surf_args{:});

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
