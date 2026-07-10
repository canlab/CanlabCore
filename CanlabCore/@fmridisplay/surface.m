function obj = surface(obj, varargin)
% Add one or more brain surfaces to an fmridisplay object (managed display).
%
% Each surface is registered as a view on the object, so blobs added with
% addblobs (and refresh / removeblobs / rethreshold) act on all surfaces and
% montages together. The actual surface geometry and default camera view come
% from addbrain, so ANY surface keyword addbrain understands works here.
%
% :Usage:
% ::
%
%     obj = surface(obj, 'direction', DIR, 'orientation', ORN, 'axes', POS)
%     obj = surface(obj, DIR)          % bare direction token (e.g. 'thalamus')
%
% :Inputs:
%
%   **obj:**
%        an fmridisplay object (handle). Initialize first, e.g. with montage.
%
% :Optional Inputs:
%
%   **'direction', DIR:**
%        Which surface to draw. DIR is passed through to addbrain, so the full
%        set of addbrain surface/region keywords is available (run `help
%        addbrain` for the complete, current list). DIR may also be given as a
%        bare token: surface(obj, 'thalamus') == surface(obj, 'direction',
%        'thalamus'). Common cortical choices:
%          'hires left' / 'hires right'              (default: 'hires right')
%          'inflated left' / 'inflated right'        (freesurfer/fsaverage)
%          'hcp inflated left' / 'hcp inflated right'
%          'freesurfer inflated left' / '... right'
%          'flat left' / 'flat right'
%        Composite / subcortical examples (all from addbrain):
%          'brainstem left' / 'brainstem right'      (brainstem + midbrain group)
%          'caudate left' / 'caudate right'          (basal ganglia group)
%          'thalamus', 'amygdala', 'cerebellum', 'limbic', 'bg', ...
%        Multi-surface keywords add a SET of views to this same object:
%          'foursurfaces', 'foursurfaces_hcp', 'foursurfaces_freesurfer'
%
%   **'orientation', ORN:**
%        'lateral' (default) or 'medial'. For an L/R cortical surface, 'medial'
%        mirrors the camera azimuth 180 degrees. Ignored for surfaces without a
%        meaningful medial view.
%
%   **'axes', POS:**
%        An existing axes handle, or a [left bottom width height] position
%        vector for a new axes. Default: current axes (gca).
%
%   Any remaining inputs are passed straight through to addbrain.
%
% :Outputs:
%
%   **obj:**
%        the fmridisplay object, with the new surface view(s) registered in
%        obj.surface and any existing blob layers rendered onto them.
%
% :Examples:
% ::
%
%     o2 = fmridisplay; o2 = montage(o2, 'axial');
%     o2 = surface(o2, 'direction', 'hires right', 'orientation', 'lateral');
%     o2 = surface(o2, 'thalamus');                 % bare addbrain keyword
%     o2 = surface(o2, 'foursurfaces_hcp');         % 2x2 set in its own figure
%
% :See also:
%   - addbrain (full surface/region keyword list), addblobs, montage, refresh,
%     removeblobs, render_on_surface, canlab_results_fmridisplay
%
% See help fmridisplay

if nargin == 0 || isempty(obj) || ~isa(obj, 'fmridisplay')
    obj = fmridisplay;
end
% initialize, if nothing passed in; but you would have to call overloaded
% method, fmridisplay.montage, to invoke this.


if isempty(obj.SPACE) || ~isstruct(obj.SPACE.V)
    error('fmridisplay is not initialized correctly. run obj = fmridisplay; first and then pass in to this method.')

end

% Surfaces cannot be drawn into a uifigure (e.g. the auto-opened controller
% window): addbrain's subplot/gca then errors ("Adding subplots to a container
% with 'AutoResizeChildren' set to 'on' is not supported"). If the active figure
% is a uifigure and the caller did not target an explicit axes, open a dedicated
% traditional figure so gca/subplot land there. (The multi-surface path below
% manages its own figure; this covers the single-surface path and pre-empties the
% current figure so the multi path opens a clean one.)
if ~any(strcmp(varargin, 'axes'))
    cf = get(groot, 'CurrentFigure');
    if ~isempty(cf) && is_uifigure(cf), figure; end
end

% An fmri_surface_data argument: add its native surface(s) AND paint it as a
% managed surface-native layer. surface(o2, surf_obj, ...) then behaves like
% surface(surf_obj) but on the stateful display -- the surfaces are registered
% views under the controller, and the data is a real layer (set_colormap /
% set_opacity / removeblobs / refresh act on it). With no surface keyword, the
% set that MATCHES the object's space is added automatically (foursurfaces_hcp
% for fs_LR, foursurfaces_freesurfer for fsaverage), so the data always renders
% at full fidelity. An explicit non-matching surface (e.g. fsaverage meshes for
% fs_LR data) is added but cannot be painted, and add_surface_blobs warns.
% See add_surface_blobs, fmri_surface_data.surface.
is_surf = cellfun(@(a) isa(a, 'fmri_surface_data'), varargin);
if any(is_surf)
    surf_data = varargin{find(is_surf, 1)};
    rest = varargin(~is_surf);

    % Split the remaining args into surface directives (WHICH surfaces to draw)
    % and colour options (HOW to paint the layer). Colour options are the
    % value-bearing keys add_surface_blobs / canlab_colormap understand; every
    % other token (foursurfaces*, direction/orientation/axes pairs, a bare
    % addbrain surface keyword) is a surface directive.
    color_keys = {'which_image', 'clim', 'cmaprange', 'colormap', 'colormapname', ...
        'pos_colormap', 'neg_colormap', 'splitcolor', 'maxcolor', 'mincolor', ...
        'color', 'transvalue', 'wh_surface', 'wh_surfaces'};
    surf_args = {}; color_args = {};
    j = 1;
    while j <= numel(rest)
        a = rest{j};
        if ischar(a) && any(strcmpi(a, color_keys))
            color_args = [color_args, rest(j:j+1)]; j = j + 2; %#ok<AGROW>
        else
            surf_args{end + 1} = a; j = j + 1; %#ok<AGROW>
        end
    end

    if isempty(surf_args)
        surf_args = {surface_default_keyword(surf_data)};   % matching space
    end

    % Shared colour range across cortex and subcortex, so both layers use ONE
    % scale (unless the caller passed clim/cmaprange). Computed from all
    % grayordinates via the standard robust policy.
    which_img = 1;
    whi = find(strcmpi(color_args, 'which_image'), 1);
    if ~isempty(whi), which_img = color_args{whi + 1}; end
    if ~any(strcmpi(color_args, 'cmaprange')) && ~any(strcmpi(color_args, 'clim'))
        vv = double(surf_data.dat(:, min(which_img, size(surf_data.dat, 2))));
        vv = vv(vv ~= 0 & isfinite(vv));
        if ~isempty(vv)
            sf = {}; if any(strcmpi(color_args, 'splitcolor')), sf = {'splitcolor'}; end
            color_args = [color_args, {'cmaprange', canlab_default_cmaprange(vv, sf{:})}];
        end
    end

    n0 = numel(obj.surface);
    obj = surface(obj, surf_args{:});                       % add the view(s)
    new_idx = (n0 + 1):numel(obj.surface);
    if isempty(new_idx), new_idx = 1:numel(obj.surface); end

    % Cortex: surface-native layer painted directly on matching meshes.
    obj = add_surface_blobs(obj, surf_data, color_args{:}, 'wh_surface', new_idx);

    % Subcortex (if present): its grayordinate voxels have no surface, so render
    % them as a standard volume layer. On these surface views that paints the
    % subcortical anatomical meshes (thalamus, brainstem, cerebellum, ...) that
    % the foursurfaces set draws, via the volume->surface projection; on a montage
    % it shows the slices. It is a normal volume layer, so rethreshold / opacity /
    % set_colormap / the controller all drive it too.
    if surface_has_volume(surf_data)
        subvol = get_wh_image(to_fmri_data(surf_data), which_img);
        vol_color = volume_color_args(color_args);
        obj = addblobs(obj, subvol, vol_color{:}, 'wh_surface', new_idx);
        % Confine this layer to the subcortical meshes: it must not repaint the
        % cortical surfaces (owned by the surface-native layer above). Recomposite
        % so the cortex is restored to the surface-native colouring / medial-wall
        % gray after addblobs' incremental paint.
        obj.activation_maps{end}.skip_cortex_nv = cortex_vertex_counts(surf_data);
        obj = composite_surfaces(obj, new_idx);
    end
    return
end

% Multi-surface keywords (e.g. 'foursurfaces', 'foursurfaces_hcp') add a SET
% of surface views to THIS SAME object, laid out in the current figure. Each
% becomes its own registered view, so blobs/refresh act on all of them
% together. (fmridisplay is a handle class, so han is not replaced.)
multi = '';
for kw = {'foursurfaces_hcp', 'foursurfaces_freesurfer', 'foursurfaces'}
    if any(strcmp(varargin, kw{1})), multi = kw{1}; break; end
end
if ~isempty(multi)
    obj = add_multiple_surfaces(obj, multi, varargin{:});
    return
end

% Parse inputs. Three reserved keyword pairs ('direction', 'orientation',
% 'axes') are consumed here; the FIRST unrecognized char token is taken as a
% BARE direction (so surface(han, 'thalamus') works, not just surface(han,
% 'direction', 'thalamus')), and any remaining args pass straight through to
% addbrain. This is what lets ANY eligible addbrain surface keyword be used in
% the managed display (see `help addbrain` for the full surface/region list).
ax = gca;
dir = '';
orn = 'lateral';
extra = {};                                            % pass-through to addbrain
i = 1;
while i <= numel(varargin)
    a = varargin{i};
    if ischar(a)
        switch a
            case 'direction',   dir = varargin{i + 1}; i = i + 2; continue
            case 'orientation', orn = varargin{i + 1}; i = i + 2; continue
            case 'axes',        ax  = varargin{i + 1}; i = i + 2; continue
            otherwise
                if isempty(dir)
                    dir = a;                           % bare direction token
                else
                    extra{end + 1} = a; %#ok<AGROW>    % addbrain pass-through
                end
                i = i + 1; continue
        end
    else
        extra{end + 1} = a; %#ok<AGROW>                % non-char (e.g. addbrain flag value)
        i = i + 1;
    end
end
if isempty(dir), dir = 'hires right'; end

if isa(ax,'matlab.graphics.axis.Axes')
    axh = ax;
    axes(axh);
else
    axh = axes('Position', ax);
end

% Build the surface(s). addbrain owns the surface geometry AND the default
% (lateral) camera view + lighting for every keyword, including its many
% composites (brainstem/caudate groups, limbic, basal ganglia, ...), so this
% method is essentially a pass-through: ANY eligible addbrain surface works
% here automatically. We only (1) special-case 'bigbrain left/right', which map
% to a DIFFERENT addbrain surface than 'bigbrain left' would, and (2) mirror the
% azimuth 180 degrees for a medial view (medial = lateral azimuth + 180).
switch dir
    case 'bigbrain left'
        h = addbrain('bigbrain', extra{:});
        view(-137, 18); lightRestoreSingle;

    case 'bigbrain right'
        h = addbrain('bigbrain', extra{:});
        view(137, 18); lightRestoreSingle;

    otherwise
        try
            h = addbrain(dir, extra{:});
        catch err
            % addbrain rejects unknown keywords with a terse 'Unknown method.';
            % translate that into an informative message that names the offending
            % token and points at the full list (the rest of addbrain's errors,
            % e.g. a missing surface file, are passed through unchanged).
            if strcmp(err.message, 'Unknown method.')
                error('fmridisplay:surface:unknownDirection', ...
                    ['surface: ''%s'' is not a recognized addbrain surface/region keyword.\n' ...
                     'Run `help addbrain` for the full list. Common choices: ''hires left''/''hires right'', ' ...
                     '''inflated left''/''inflated right'', ''flat left''/''flat right''; cutaways are ' ...
                     '''left_cutaway''/''right_cutaway'' (note the order, not ''cutaway_left''); composites ' ...
                     'include ''insula surfaces'', ''brainstem left'', ''caudate right'', ''limbic''.'], dir);
            else
                rethrow(err);
            end
        end
        if strcmp(orn, 'medial')
            [az, el] = view(axh);
            view(axh, mod(az + 180, 360), el);
        end
end





% Normalize face colouring per object. addbrain may return modern graphics
% objects OR old-style NUMERIC handles (doubles) — dot indexing (h(1).FaceColor)
% fails on the latter — and a single request can mix patch types with different
% colouring: e.g. 'coronal_slabs_4' interleaves isosurface slabs (a flat
% FaceColor) with isocaps that use 'interp' to show the anatomy cross-section.
% So use get/set and decide per object: keep any 'interp'/textured colouring,
% set everything else to flat gray. Make all faces opaque.
for ii = 1:numel(h)
    hh = h(ii);
    if ~isprop(hh, 'FaceColor'), continue, end
    fc = get(hh, 'FaceColor');
    if ~(ischar(fc) && strcmp(fc, 'interp'))
        set(hh, 'FaceColor', [.5 .5 .5]);
    end
    set(hh, 'FaceAlpha', 1);
end
%[
lightRestoreSingle;
axis image;
axis off;
lighting gouraud;
material dull;

% Hide the "..." interaction toolbar. Sweep the whole figure, not just axh, since
% some keywords (cutaways, coronal_slabs_*) draw their patches into extra axes
% they create themselves rather than into axh.
canlab_hide_axes_toolbar(ancestor(axh, 'figure'));


% register info in fmridisplay object
% 'axis_handles', axh,
obj.surface{end + 1} = struct('axis_handles', axh, 'direction', dir, 'orientation', orn, 'object_handle', h);

% Pull-in: if blob layers already exist on this object, render them onto the
% surface we just added, so a surface added AFTER addblobs still shows the
% blobs and add/remove/refresh act on all views together (design notes 4.2).
new_surf_idx = numel(obj.surface);
for k = 1:numel(obj.activation_maps)
    obj = render_layer_surfaces(obj, k, new_surf_idx);
end

end % main function


function obj = add_multiple_surfaces(obj, multi, varargin)
% Add a canonical set of four surface views (L/R lateral + L/R medial) to the
% same fmridisplay object, laid out 2x2 in a dedicated figure.
%
% For the keywords addbrain builds directly ('foursurfaces' and
% 'foursurfaces_hcp'), the geometry is routed THROUGH addbrain so the surfaces
% are identical to what fmri_data.surface / region.surface produce (same
% meshes, subplot layout, camera views, lighting, and the added thalamus).
% Each addbrain subplot is then registered as its own fmridisplay view, so
% blobs added with addblobs (and refresh / removeblobs) act on all four panels.
% 'foursurfaces_freesurfer' is NOT an addbrain keyword, so it is still built
% here from individual inflated surfaces via single-surface surface() calls.

% Ensure a dedicated figure for the 2x2 layout. If the current figure already
% has content (a montage's axes, OR a uifigure like the controller with uipanel
% children), the four surfaces would draw over it or error, so open a fresh one.
% Reuse only a genuinely empty traditional figure (e.g. the user typed `figure;`).
curfig = get(groot, 'CurrentFigure');
if isempty(curfig) || is_uifigure(curfig) || ~isempty(allchild(curfig))
    figure;
end
layout_fig = gcf;
set(layout_fig, 'Color', 'w');                       % white background behind the surfaces

rest = varargin;
rest(strcmp(rest, multi)) = [];                       % drop the multi keyword
for kw = {'direction', 'orientation', 'axes'}         % drop any single-view controls
    wh = find(strcmp(rest, kw{1}));
    for j = sort(wh, 'descend'), rest(j:j+1) = []; end
end

if any(strcmp(multi, {'foursurfaces', 'foursurfaces_hcp'}))

    % Build the canonical layout through addbrain (drawn into layout_fig via
    % subplots), so the surfaces match fmri_data.surface exactly.
    h = addbrain(multi, rest{:});
    h = h(ishandle(h));

    % Register one fmridisplay view per axes addbrain populated. Group the
    % returned patch handles by their parent axes, preserving creation order,
    % so each panel (including its thalamus) becomes a view that addblobs and
    % refresh render onto.
    ax_list = gobjects(0);
    grp = zeros(numel(h), 1);
    for ii = 1:numel(h)
        a = ancestor(h(ii), 'axes');
        idx = find(ax_list == a, 1);
        if isempty(idx), ax_list(end + 1) = a; idx = numel(ax_list); end %#ok<AGROW>
        grp(ii) = idx;
    end

    for a = 1:numel(ax_list)
        these = h(grp == a);
        obj.surface{end + 1} = struct('axis_handles', ax_list(a), ...
            'direction', multi, 'orientation', '', 'object_handle', these(:)');

        % Pull any existing blob layers onto the newly registered surface, so
        % a surface added after addblobs still shows the blobs (design notes 4.2).
        new_surf_idx = numel(obj.surface);
        for k = 1:numel(obj.activation_maps)
            obj = render_layer_surfaces(obj, k, new_surf_idx);
        end
    end

else

    % 'foursurfaces_freesurfer' is not an addbrain keyword, so build the 2x2
    % from individual inflated surfaces. {direction, orientation, axes-position}:
    % top row = lateral views, bottom row = medial views; left col = left hemi.
    Lsurf = 'freesurfer inflated left'; Rsurf = 'freesurfer inflated right';
    configs = { ...
        {Lsurf, 'lateral', [0.03 0.50 0.45 0.45]}, ...
        {Rsurf, 'lateral', [0.52 0.50 0.45 0.45]}, ...
        {Lsurf, 'medial',  [0.03 0.03 0.45 0.45]}, ...
        {Rsurf, 'medial',  [0.52 0.03 0.45 0.45]} };

    for c = 1:numel(configs)
        cfg = configs{c};
        obj = surface(obj, 'direction', cfg{1}, 'orientation', cfg{2}, 'axes', cfg{3}, rest{:});
    end

end

% Turn off the axis box/background for EVERY axes in the layout figure — the
% surface axes, the hidden colorbar axes, and any stray default axes created
% during rendering (which otherwise shows a frame behind the surfaces). The
% surface patches are children of the axes and stay visible; colorbars are
% separate objects and are unaffected.
set(findobj(layout_fig, 'Type', 'axes'), 'Visible', 'off');
canlab_hide_axes_toolbar(layout_fig);   % hide the "..." toolbar on all surface panels

end


function tf = is_uifigure(f)
% Thin wrapper on the shared canlab_is_uifigure (keeps figure-creation behaviour
% consistent with image_vector.montage / canlab_results_fmridisplay).
tf = canlab_is_uifigure(f);
end


function tf = surface_has_volume(surf_data)
% True when the grayordinate object has a subcortical (voxel) sub-block.
tf = ~isempty(surf_data.brain_model) && ...
     any(cellfun(@(m) strcmp(m.type, 'vox'), surf_data.brain_model.models));
end


function nv = cortex_vertex_counts(surf_data)
% Per-hemisphere vertex counts of the object's cortical models (e.g. 32492 for
% fs_LR-32k). Used to keep a subcortical-volume layer off the cortical meshes.
nv = [];
for i = 1:numel(surf_data.brain_model.models)
    m = surf_data.brain_model.models{i};
    if strcmp(m.type, 'surf') && ~isempty(m.numvert) && ~isnan(m.numvert)
        nv(end + 1) = m.numvert; %#ok<AGROW>
    end
end
nv = unique(nv);
end


function out = volume_color_args(color_args)
% Translate the surface color options into the subset the volume blob path
% (render_blobs) understands, so cortex and subcortex share a colour spec.
% which_image is already applied; clim -> cmaprange; drop surface-only tokens.
out = {};
i = 1;
while i <= numel(color_args)
    key = color_args{i};
    if ~ischar(key), i = i + 1; continue; end
    switch lower(key)
        case 'clim'
            out = [out, {'cmaprange', color_args{i + 1}}]; %#ok<AGROW>
        case {'cmaprange', 'maxcolor', 'mincolor', 'splitcolor', 'color', ...
              'pos_colormap', 'neg_colormap'}
            out = [out, color_args(i:i + 1)]; %#ok<AGROW>
        case 'colormap'
            if i + 1 <= numel(color_args) && isnumeric(color_args{i + 1})
                out = [out, color_args(i:i + 1)]; %#ok<AGROW>
            end
        % which_image / colormapname / transvalue / wh_surface: not forwarded.
    end
    i = i + 2;
end
end


function kw = surface_default_keyword(surf_data)
% Multi-surface keyword whose meshes MATCH the object's surface space, so the
% data paints directly (no cross-space resampling). Used when surface(o2,
% surf_obj) is called with no explicit surface directive.
switch surf_data.surface_space
    case 'fsLR_32k',       kw = 'foursurfaces_hcp';          % 32492 verts/hemi
    case 'fsaverage_164k', kw = 'foursurfaces_freesurfer';   % 163842 verts/hemi
    otherwise,             kw = 'foursurfaces_hcp';          % best-effort default
end
end

