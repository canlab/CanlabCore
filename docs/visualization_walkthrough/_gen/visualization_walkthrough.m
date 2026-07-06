function visualization_walkthrough
% VISUALIZATION_WALKTHROUGH  Recreate every figure in the CANlab visualization walkthrough.
%
% One section per walkthrough page (1.1, 1.2, ... 6.5). Running this script with
% CanlabCore (and Neuroimaging_Pattern_Masks, for atlases) on the path regenerates
% all PNGs under ../figures/. It is also meant to be read: each section is the
% copy-pasteable code shown on the corresponding walkthrough page, plus the one
% wt_save() line that captures the figure.
%
% Pages:
%   1. Getting started        docs/visualization_walkthrough/01_getting_started.md
%   2. Montages and slices    .../02_montages.md
%   3. Surfaces and 3-D       .../03_surfaces.md
%   4. Colors and colormaps   .../04_colormaps.md
%   5. The display controller .../05_controller.md
%   6. Atlases and regions    .../06_atlases_and_regions.md
%
% This supersedes the per-section gen_*.m scripts; those remain as small,
% independently runnable per-page generators.

here   = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(fileparts(here), 'figures');
if ~exist(figdir, 'dir'), mkdir(figdir); end

% Reload edited classdefs cleanly if instances linger from a prior run.
close all force

% Shared handles: capture the RIGHT figure for multi-figure display methods.
montfig = @(o) ancestor(o.montage{1}.axis_handles(1), 'figure');
surffig = @(o) ancestor(local_surfh(o), 'figure');

% Running dataset used throughout (bundled; no download needed).
obj     = load_image_set('emotionreg', 'noverbose');   % 30 contrast images (fmri_data)
t       = threshold(ttest(obj), .05, 'unc');           % statistic_image, p < .05 unc
tstrict = threshold(ttest(obj), .001, 'unc');          % stricter map for overlays
r       = region(t);                                   % contiguous blobs


%% ===== 1. Getting started =================================================

%% 1.1 orthviews (SPM-based three-plane viewer)
orthviews(t);
gwin = spm_figure('FindWin', 'Graphics');
if ~isempty(gwin), wt_save(gwin, fullfile(figdir, '01_orthviews.png')); end
close all force

%% 1.2 canlab_orthviews (SPM-free three-plane viewer)
canlab_orthviews(t);
wt_save(gcf, fullfile(figdir, '01_canlab_orthviews.png'));
close all force

% (canlab_niivue(t) opens an interactive browser viewer; nothing to capture.)


%% ===== 2. Montages and slices =============================================

%% 2.2 Prepackaged layouts: compact2 and full
o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '02_compact2.png'));
close all force

o = canlab_results_fmridisplay(t, 'full', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '02_full.png'));
close all force

%% 2.3 regioncenters: each blob in its own axis, colored by value
o = montage(r, 'regioncenters', 'colormap', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '02_regioncenters.png'));
close all force

%% 2.4 Customizing blobs: red outline + transparency
o = canlab_results_fmridisplay(t, 'compact2', 'color', [1 0 0], 'outline', ...
    'linewidth', 2, 'trans', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '02_outline_trans.png'));
close all force

%% 2.5 Overlaying two maps in two colors
o = canlab_results_fmridisplay(t, 'compact2', 'color', [1 .3 0], 'noverbose');
o = addblobs(o, region(tstrict), 'color', [0 .4 1], 'no_surface');
wt_save(montfig(o), fullfile(figdir, '02_two_maps.png'));
close all force

%% 2.6 Custom montage: fmridisplay.montage with orientation / range / spacing
o = fmridisplay;
o = montage(o, 'axial', 'slice_range', [-30 60], 'spacing', 10, 'onerow');
o = addblobs(o, region(t), 'nooutline');
wt_save(montfig(o), fullfile(figdir, '02_custom_montage.png'));
close all force


%% ===== 3. Surfaces and 3-D rendering ======================================

%% 3.1 The surface() method: foursurfaces, and the default cutaway
create_figure('four');
sh = surface(t, 'foursurfaces', 'noverbose');          % returns surface handles
wt_save(gcf, fullfile(figdir, '03_foursurfaces.png'), 'nobars');

%% 3.2 Recoloring existing handles with render_on_surface (summer LUT)
render_on_surface(t, sh, 'colormap', 'summer', 'nolegend');
wt_save(gcf, fullfile(figdir, '03_foursurfaces_summer.png'), 'nobars');
close all force

create_figure('def'); axis off
surface(t, 'noverbose');                               % no layout -> cutaway
wt_save(gcf, fullfile(figdir, '03_surface_default.png'), 'nobars');
close all force

%% 3.3 Cutaways and slabs: coronal_slabs_4
create_figure('cs'); axis off
surface(t, 'coronal_slabs_4', 'noverbose');
wt_save(gcf, fullfile(figdir, '03_coronal_slabs_4.png'), 'nobars');
close all force

%% 3.4 Subcortical close-up composed with addbrain
create_figure('sub'); axis off
sh = addbrain('thalamus_group');
sh = [sh addbrain('brainstem')];
sh = [sh addbrain('amygdala')];
render_on_surface(t, sh, 'clim', [-3 3], 'nolegend');
set(sh, 'FaceAlpha', .9); view(135, 12); lightRestoreSingle; camzoom(1.3);
wt_save(gcf, fullfile(figdir, '03_subcortical_closeup.png'), 'nobars');
close all force

%% 3.5 Isosurfaces from an atlas: bare nuclei, then t-map rendered onto them
create_figure('iso'); axis off
atl = load_atlas('thalamus');
sh  = isosurface(atl);
axis image vis3d off; material dull; view(210, 20); lightRestoreSingle
wt_save(gcf, fullfile(figdir, '03_isosurface_thalamus.png'), 'nobars');
render_on_surface(t, sh, 'colormap', 'hot', 'nolegend');
wt_save(gcf, fullfile(figdir, '03_isosurface_thalamus_rendered.png'), 'nobars');
close all force


%% ===== 4. Colors and colormaps ============================================
% One fmridisplay object, restyled in place: split -> mango -> warm ramp ->
% continuous LUT -> solid. The value->color source drives montage and legend.

o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '04_hotcool.png'));           % default split

set_colormap(o, 'splitcolor', {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]});
wt_save(montfig(o), fullfile(figdir, '04_mango.png'));             % mango split

set_colormap(o, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);
wt_save(montfig(o), fullfile(figdir, '04_warm.png'));             % single warm ramp

set_colormap(o, 'colormap', hot(256));
wt_save(montfig(o), fullfile(figdir, '04_perceptual.png'));        % continuous LUT

set_colormap(o, 'color', [1 .4 .9]);
wt_save(montfig(o), fullfile(figdir, '04_solid.png'));            % solid color
close all force

%% 4.x Color limits: same data, two cmaprange settings
o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '04_cmaprange_default.png'));
set_colormap(o, 'cmaprange', [-6 6]);
wt_save(montfig(o), fullfile(figdir, '04_cmaprange_wide.png'));
close all force

%% 4.x Uniformity: one color source drives montage AND surface identically
o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
o = surface(o, 'foursurfaces');
wt_save(montfig(o), fullfile(figdir, '04_uniform_montage.png'));
sh = o.surface{1}.object_handle; sh = sh(ishandle(sh));
wt_save(ancestor(sh(1), 'figure'), fullfile(figdir, '04_uniform_surface.png'));
close all force


%% ===== 5. The display controller ==========================================

%% 5.1 Build a composite (montage + surface + 2nd layer) and open the controller
o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
o = surface(o, 'foursurfaces');
o = addblobs(o, region(tstrict), 'color', [0 .4 1]);   % 2nd layer, solid blue
cf = controller(o);
wt_save(cf,          fullfile(figdir, '05_controller.png'), 'app');
wt_save(montfig(o),  fullfile(figdir, '05_montage.png'));
wt_save(surffig(o),  fullfile(figdir, '05_surface.png'));

%% 5.3 Command line == GUI: recolor layer 1 (mango), fade layer 2; views sync
set_colormap(o, 'splitcolor', {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]}, 'layers', 1);
set_opacity(o, 0.4, 'layers', 2);
cf = controller(o);
wt_save(cf,         fullfile(figdir, '05_controller_after.png'), 'app');
wt_save(montfig(o), fullfile(figdir, '05_montage_after.png'));
close all force

%% 5.4 multiview re-composes the object's layers onto a canonical layout
% o = multiview(o, 'full');   % (interactive; not captured here)


%% ===== 6. Atlases and regions =============================================

%% 6.1 region(t) in unique per-blob colors
o = montage(region(t), 'compact2', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '06_region_unique.png'));
close all force

%% 6.2 Atlas on slices (thalamic nuclei, unique colors)
atl = load_atlas('thalamus');
o = montage(atl, 'compact2', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '06_atlas_montage.png'));
close all force

%% 6.4 Atlas as 3-D isosurfaces (subcortical close-up)
create_figure('iso'); axis off
isosurface(load_atlas('thalamus'));
axis image vis3d off; material dull; view(210, 20); lightRestoreSingle; camzoom(1.4);
wt_save(gcf, fullfile(figdir, '06_atlas_isosurface.png'));
close all force

%% 6.3 Cortical atlas on surfaces in unique (indexed) colors
ctx  = select_atlas_subset(load_atlas('canlab2018_2mm'), {'Ctx'});
cols = scn_standard_colors(num_regions(ctx));
cmap = cat(1, cols{:});                                % n-by-3 color table
create_figure('atlsurf');
sh = addbrain('foursurfaces');
render_on_surface(ctx, sh, 'colormap', cmap, 'indexmap', 'interp', 'nearest', 'nolegend');
wt_save(gcf, fullfile(figdir, '06_atlas_surface.png'), 'nobars');
close all force

disp('visualization_walkthrough: ALL figures regenerated.');
end

% -------------------------------------------------------------------------
function h = local_surfh(o)
% First valid surface graphics handle on an fmridisplay object.
sh = o.surface{1}.object_handle; sh = sh(ishandle(sh)); h = sh(1);
end
