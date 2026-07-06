function gen_05_controller
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(fileparts(here), 'figures');
close all force

obj = load_image_set('emotionreg', 'noverbose');
t   = threshold(ttest(obj), .05, 'unc');
tstrict = threshold(ttest(obj), .001, 'unc');

montfig = @(o) ancestor(o.montage{1}.axis_handles(1), 'figure');
surffig = @(o) ancestor(local_surfh(o), 'figure');

% Build a composite: montage + surface, two blob layers
o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
o = surface(o, 'foursurfaces');
o = addblobs(o, region(tstrict), 'color', [0 .4 1]);   % 2nd layer, solid blue
cf = controller(o);

wt_save(cf, fullfile(figdir, '05_controller.png'), 'app');
wt_save(montfig(o), fullfile(figdir, '05_montage.png'));
wt_save(surffig(o), fullfile(figdir, '05_surface.png'));

% Command-line action that mirrors a GUI control: recolor layer 1, drop
% opacity of layer 2. Montage AND surface update in place.
set_colormap(o, 'splitcolor', {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]}, 'layers', 1);
set_opacity(o, 0.4, 'layers', 2);
cf = controller(o);
wt_save(cf, fullfile(figdir, '05_controller_after.png'), 'app');
wt_save(montfig(o), fullfile(figdir, '05_montage_after.png'));
close all force

disp('gen_05 DONE');
end

function h = local_surfh(o)
sh = o.surface{1}.object_handle; sh = sh(ishandle(sh)); h = sh(1);
end
