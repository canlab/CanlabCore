function gen_02_montages
% Figures for "2. Montages and slices".
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(fileparts(here), 'figures');
close all force
montfig = @(o) ancestor(o.montage{1}.axis_handles(1), 'figure');

obj = load_image_set('emotionreg', 'noverbose');
t   = threshold(ttest(obj), .05, 'unc');
r   = region(t);

% compact2 and full layouts
o = canlab_results_fmridisplay(t, 'compact2', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '02_compact2.png'));
close all force

o = canlab_results_fmridisplay(t, 'full', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '02_full.png'));
close all force

% regioncenters: each blob in its own axis, coloured by value
o = montage(r, 'regioncenters', 'colormap', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '02_regioncenters.png'));
close all force

% customization: red outline + transparency
o = canlab_results_fmridisplay(t, 'compact2', 'color', [1 0 0], 'outline', ...
    'linewidth', 2, 'trans', 'noverbose');
wt_save(montfig(o), fullfile(figdir, '02_outline_trans.png'));
close all force

% overlay two maps in two colours
o = canlab_results_fmridisplay(t, 'compact2', 'color', [1 .3 0], 'noverbose');
tstrict = threshold(ttest(obj), .001, 'unc');
o = addblobs(o, region(tstrict), 'color', [0 .4 1], 'no_surface');
wt_save(montfig(o), fullfile(figdir, '02_two_maps.png'));
close all force

disp('gen_02 DONE');
end
