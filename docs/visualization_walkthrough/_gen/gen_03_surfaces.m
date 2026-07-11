function gen_03_surfaces
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(fileparts(here), 'figures');
close all force

obj = load_image_set('emotionreg', 'noverbose');
t   = threshold(ttest(obj), .05, 'unc');

% --- default surface() : a preset series incl. subcortical cutaway ----------
create_figure('surf'); axis off
surface(t, 'noverbose');
wt_save(gcf, fullfile(figdir, '03_surface_default.png'));
close all force

% --- foursurfaces : lateral + medial cortical set ---------------------------
create_figure('four');
sh = surface(t, 'foursurfaces', 'noverbose');
wt_save(gcf, fullfile(figdir, '03_foursurfaces.png'));
% recolor the SAME handles with a perceptual colormap
render_on_surface(t, sh, 'colormap', 'summer');
wt_save(gcf, fullfile(figdir, '03_foursurfaces_summer.png'));
close all force

% --- cutaways ---------------------------------------------------------------
create_figure('lc'); axis off
surface(t, 'left_cutaway', 'noverbose');
wt_save(gcf, fullfile(figdir, '03_left_cutaway.png'));
close all force

create_figure('cs'); axis off
surface(t, 'coronal_slabs_4', 'noverbose');
wt_save(gcf, fullfile(figdir, '03_coronal_slabs_4.png'));
close all force

% --- subcortical close-up via addbrain + render_on_surface ------------------
create_figure('sub'); axis off
sh = addbrain('thalamus_group');
sh = [sh addbrain('brainstem')];
sh = [sh addbrain('amygdala')];
render_on_surface(t, sh, 'clim', [-3 3]);
set(sh, 'FaceAlpha', .9);
view(135, 12); lightRestoreSingle; camzoom(1.3);
wt_save(gcf, fullfile(figdir, '03_subcortical_closeup.png'));
close all force

% --- isosurface of an atlas (thalamic nuclei), coloured by the t-map --------
create_figure('iso'); axis off
atl = load_atlas('thalamus');
sh = isosurface(atl);
axis image vis3d off; material dull; view(210, 20); lightRestoreSingle
wt_save(gcf, fullfile(figdir, '03_isosurface_thalamus.png'));
render_on_surface(t, sh, 'colormap', 'hot');
wt_save(gcf, fullfile(figdir, '03_isosurface_thalamus_rendered.png'));
close all force

disp('gen_03 DONE');
end
