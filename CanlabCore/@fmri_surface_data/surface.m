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
%   1. NATIVE (default): loads the bundled mesh that matches the object's
%      surface_space (fs_LR-32k or fsaverage-164k) and colors vertices DIRECTLY
%      from the data -- no resampling. Produces a 4-panel figure (L/R x
%      lateral/medial). The medial wall and zero values render gray.
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
%   **'clim':**         [lo hi] color limits (default symmetric from data).
%   **'pos_colormap' / 'neg_colormap':** [n x 3] colormaps (default hot / cool).
%   **'existingsurface':** vector of patch handles to color.
%   **'mni_surface':**  an addbrain keyword (string).
%
% :Outputs:
%   **han:** struct with fields .figure, .axes, .surfaces (graphics handles).
%
% :Examples:
% ::
%     s = fmri_surface_data(which('transcriptomic_gradients.dscalar.nii'));
%     surface(s, 'which_image', 1);
%
%     v = fmri_data(which('weights_NSF_grouppred_cvpcr.img'));
%     surface(vol2surf(v));                      % native fsaverage render
%     surface(vol2surf(v), 'mni_surface', 'left'); % on an addbrain MNI surface
%
% :See also: render_on_surface, load_surface_geom, addbrain, fmri_surface_data

surftype = 'inflated';
which_image = 1;
clim = [];
poscm = [];
negcm = [];
existing = [];
mni_surface = '';

i = 1;
while i <= numel(varargin)
    switch lower(varargin{i})
        case 'surftype',        surftype = varargin{i+1};   i = i + 2;
        case 'which_image',     which_image = varargin{i+1}; i = i + 2;
        case 'clim',            clim = varargin{i+1};       i = i + 2;
        case 'pos_colormap',    poscm = varargin{i+1};      i = i + 2;
        case 'neg_colormap',    negcm = varargin{i+1};      i = i + 2;
        case {'existingsurface','surface_handles'}, existing = varargin{i+1}; i = i + 2;
        case 'mni_surface',     mni_surface = varargin{i+1}; i = i + 2;
        otherwise
            error('fmri_surface_data:surface:badopt', 'Unknown option: %s', num2str(varargin{i}));
    end
end

ropts = {'which_image', which_image, 'clim', clim, 'pos_colormap', poscm, 'neg_colormap', negcm};

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

% ---- Mode 1: native 4-panel render ----
geom = load_surface_geom(obj.surface_space, surftype);

panels = { 'lh', [270 0], 1; ...   % left lateral
           'lh', [90 0],  2; ...   % left medial
           'rh', [90 0],  3; ...   % right lateral
           'rh', [270 0], 4 };     % right medial

fig = figure('Color', 'w', 'Name', sprintf('fmri_surface_data: %s (%s)', obj.surface_space, surftype));
surfaces = gobjects(1, 4);
axes_h = gobjects(1, 4);
for p = 1:4
    ax = subplot(2, 2, panels{p, 3});
    hold(ax, 'on');
    if strcmp(panels{p, 1}, 'lh')
        V = geom.vertices_lh; F = geom.faces_lh; tag = 'left';
    else
        V = geom.vertices_rh; F = geom.faces_rh; tag = 'right';
    end
    hp = patch('Parent', ax, 'Faces', F, 'Vertices', V, 'EdgeColor', 'none', ...
        'FaceColor', [.5 .5 .5], 'Tag', tag, 'SpecularStrength', .2, 'SpecularExponent', 200);
    axis(ax, 'off', 'image', 'vis3d');
    view(ax, panels{p, 2}(1), panels{p, 2}(2));
    try, lightRestoreSingle(ax); catch, camlight(ax); end %#ok<NOCOM>
    material(ax, 'dull');
    surfaces(p) = hp;
    axes_h(p) = ax;
end

render_on_surface(obj, surfaces, ropts{:});

han = struct('figure', fig, 'axes', axes_h, 'surfaces', surfaces);
end
