function gen_04_colormaps
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(fileparts(here), 'figures');
close all force

obj = load_image_set('emotionreg', 'noverbose');
t   = threshold(ttest(obj), .05, 'unc');

montfig = @(o) ancestor(o.montage{1}.axis_handles(1), 'figure');

% Default split hot/cool
o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '04_hotcool.png'));

% Mango split
set_colormap(o, 'splitcolor', {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]});
wt_save(montfig(o), fullfile(figdir, '04_mango.png'));

% Single warm ramp (mincolor -> maxcolor)
set_colormap(o, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);
wt_save(montfig(o), fullfile(figdir, '04_warm.png'));

% Continuous perceptual LUT
set_colormap(o, 'colormap', hot(256));
wt_save(montfig(o), fullfile(figdir, '04_perceptual.png'));

% Solid colour
set_colormap(o, 'color', [1 .4 .9]);
wt_save(montfig(o), fullfile(figdir, '04_solid.png'));
close all force

% cmaprange: same data, two colour limits (default vs tight)
o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '04_cmaprange_default.png'));
set_colormap(o, 'cmaprange', [-6 6]);
wt_save(montfig(o), fullfile(figdir, '04_cmaprange_wide.png'));
close all force

% Uniformity: SAME map, montage AND surface, one colour source
o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
o = surface(o, 'foursurfaces');
wt_save(montfig(o), fullfile(figdir, '04_uniform_montage.png'));
sh = o.surface{1}.object_handle; sh = sh(ishandle(sh));
wt_save(ancestor(sh(1), 'figure'), fullfile(figdir, '04_uniform_surface.png'));
close all force

disp('gen_04 DONE');
end
