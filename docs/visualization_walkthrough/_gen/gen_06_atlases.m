function gen_06_atlases
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(fileparts(here), 'figures');
close all force
montfig = @(o) ancestor(o.montage{1}.axis_handles(1), 'figure');

% region(t) in unique colors -------------------------------------------------
try
    obj = load_image_set('emotionreg', 'noverbose');
    t   = threshold(ttest(obj), .05, 'unc');
    o = montage(region(t), 'compact2', 'noverbose');   % region -> unique per-blob colors
    wt_save(montfig(o), fullfile(figdir, '06_region_unique.png'));
    close all force
catch e, fprintf('region_unique FAILED: %s\n', e.message); end

% Atlas montage in unique colors (subcortical) -------------------------------
try
    atl = load_atlas('thalamus');
    o = montage(atl, 'compact2', 'noverbose');
    wt_save(montfig(o), fullfile(figdir, '06_atlas_montage.png'));
    close all force
catch e, fprintf('atlas_montage FAILED: %s\n', e.message); end

% Atlas isosurface in unique colors (subcortical 3-D close-up) ----------------
try
    atl = load_atlas('thalamus');
    create_figure('iso'); axis off
    isosurface(atl);
    axis image vis3d off; material dull; view(210, 20); lightRestoreSingle; camzoom(1.4);
    wt_save(gcf, fullfile(figdir, '06_atlas_isosurface.png'));
    close all force
catch e, fprintf('atlas_isosurface FAILED: %s\n', e.message); end

% Cortical atlas on surfaces in unique colors --------------------------------
try
    ctx  = select_atlas_subset(load_atlas('canlab2018_2mm'), {'Ctx'});
    cmap = cat(1, scn_standard_colors(num_regions(ctx)){:});
    create_figure('atlsurf');
    sh = addbrain('foursurfaces');
    render_on_surface(ctx, sh, 'colormap', cmap, 'indexmap', 'interp', 'nearest');
    wt_save(gcf, fullfile(figdir, '06_atlas_surface.png'));
    close all force
catch e, fprintf('atlas_surface FAILED: %s\n', e.message); end

disp('gen_06 DONE');
end
