function gen_revise
% Targeted regeneration: surface figures WITHOUT colorbar legends, plus the
% canlab_orthviews and custom-montage figures. (Temporary; superseded by
% visualization_walkthrough.m.)
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(fileparts(here), 'figures');
close all force
montfig = @(o) ancestor(o.montage{1}.axis_handles(1), 'figure');

obj = load_image_set('emotionreg', 'noverbose');
t   = threshold(ttest(obj), .05, 'unc');

% --- canlab_orthviews (SPM-free 3-panel) ------------------------------------
try
    canlab_orthviews(t);
    wt_save(gcf, fullfile(figdir, '01_canlab_orthviews.png'));
    close all force
catch e, fprintf('canlab_orthviews FAILED: %s\n', e.message); end

% --- custom montage: fmridisplay.montage with slice range / spacing ---------
try
    o = fmridisplay;
    o = montage(o, 'axial', 'slice_range', [-30 60], 'spacing', 10, 'onerow');
    o = addblobs(o, region(t), 'nooutline');
    wt_save(montfig(o), fullfile(figdir, '02_custom_montage.png'));
    close all force
catch e, fprintf('custom_montage FAILED: %s\n', e.message); end

% --- surfaces WITHOUT colorbars ---------------------------------------------
create_figure('four');
sh = surface(t, 'foursurfaces', 'noverbose');
wt_save(gcf, fullfile(figdir, '03_foursurfaces.png'), 'nobars');
render_on_surface(t, sh, 'colormap', 'summer', 'nolegend');
wt_save(gcf, fullfile(figdir, '03_foursurfaces_summer.png'), 'nobars');
close all force

create_figure('def'); axis off
surface(t, 'noverbose');
wt_save(gcf, fullfile(figdir, '03_surface_default.png'), 'nobars');
close all force

create_figure('cs'); axis off
surface(t, 'coronal_slabs_4', 'noverbose');
wt_save(gcf, fullfile(figdir, '03_coronal_slabs_4.png'), 'nobars');
close all force

create_figure('sub'); axis off
sh = addbrain('thalamus_group');
sh = [sh addbrain('brainstem')];
sh = [sh addbrain('amygdala')];
render_on_surface(t, sh, 'clim', [-3 3], 'nolegend');
set(sh, 'FaceAlpha', .9); view(135, 12); lightRestoreSingle; camzoom(1.3);
wt_save(gcf, fullfile(figdir, '03_subcortical_closeup.png'), 'nobars');
close all force

create_figure('iso'); axis off
atl = load_atlas('thalamus');
sh = isosurface(atl);
axis image vis3d off; material dull; view(210, 20); lightRestoreSingle
wt_save(gcf, fullfile(figdir, '03_isosurface_thalamus.png'), 'nobars');
render_on_surface(t, sh, 'colormap', 'hot', 'nolegend');
wt_save(gcf, fullfile(figdir, '03_isosurface_thalamus_rendered.png'), 'nobars');
close all force

% --- atlas on cortical surfaces, no colorbar --------------------------------
ctx  = select_atlas_subset(load_atlas('canlab2018_2mm'), {'Ctx'});
cmap = cat(1, scn_standard_colors(num_regions(ctx)){:});
create_figure('atlsurf');
sh = addbrain('foursurfaces');
render_on_surface(ctx, sh, 'colormap', cmap, 'indexmap', 'interp', 'nearest', 'nolegend');
wt_save(gcf, fullfile(figdir, '06_atlas_surface.png'), 'nobars');
close all force

disp('gen_revise DONE');
end
