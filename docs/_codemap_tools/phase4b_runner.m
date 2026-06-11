% Phase 4b runner: produce sample PNGs for ~26 documentation pages.

parent = '/Users/f003vz1/Documents/GitHub';
addpath(genpath(fullfile(parent, 'CanlabCore')));
addpath(genpath(fullfile(parent, 'Neuroimaging_Pattern_Masks')));

png_dir = fullfile(parent, 'CanlabCore', 'docs', 'class_method_pngs');
if ~exist(png_dir, 'dir'), mkdir(png_dir); end

function save_first_fig(out_path, varargin)
    p = inputParser;
    p.addParameter('min_height', 0);
    p.addParameter('min_width', 0);
    p.parse(varargin{:});
    figs = findobj('type', 'figure');
    if isempty(figs)
        warning('No figures to save for %s', out_path);
        return
    end
    nums = arrayfun(@(f) f.Number, figs);
    [~, ix] = min(nums);
    first_fig = figs(ix);
    figure(first_fig);
    drawnow;
    pos = first_fig.Position;
    if p.Results.min_height > 0 && pos(4) < p.Results.min_height
        pos(4) = p.Results.min_height;
    end
    if p.Results.min_width > 0 && pos(3) < p.Results.min_width
        pos(3) = p.Results.min_width;
    end
    first_fig.Position = pos;
    drawnow;
    exportgraphics(first_fig, out_path, 'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('Saved %s\n', out_path);
end

% Render a multi-line text string as a PNG (used for text-output 'tables').
function save_text_as_png(txt, out_path, varargin)
    p = inputParser;
    p.addParameter('font_size', 9);
    p.addParameter('width_chars', 100);
    p.parse(varargin{:});
    fs = p.Results.font_size;
    % Strip ANSI / hyperlinks from MATLAB's pretty-print
    txt = regexprep(txt, '<a[^>]*>', '');
    txt = regexprep(txt, '</a>', '');
    txt = regexprep(txt, '\x{0008}', '');  % backspace
    lines = splitlines(txt);
    if numel(lines) > 60, lines = lines(1:60); end  % truncate very long
    n_lines = numel(lines);
    fh = figure('Color', 'white', 'Units', 'pixels', ...
                'Position', [50 50 1200 max(280, 22 * (n_lines + 2))], ...
                'MenuBar', 'none', 'ToolBar', 'none');
    ax = axes('Parent', fh, 'Position', [0 0 1 1], 'Visible', 'off', ...
              'XLim', [0 1], 'YLim', [0 1]);
    text(ax, 0.01, 0.99, strjoin(lines, newline), ...
         'FontName', 'Menlo', 'FontSize', fs, ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'Interpreter', 'none');
    drawnow;
    exportgraphics(fh, out_path, 'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('Saved %s\n', out_path);
end

% --------------- shared cached objects ---------------
imgs = load_image_set('emotionreg');
fprintf('Loaded emotionreg (%d images)\n', size(imgs.dat, 2));

% =================================================================
% IMAGE_VECTOR / FMRI_DATA METHODS
% =================================================================

% 1. image_similarity_plot — NPS-plus polar plot
try
    fprintf('\n=== image_similarity_plot ===\n');
    close all
    t_grp = ttest(imgs);  % statistic_image
    stats = image_similarity_plot(t_grp, 'mapset', 'npsplus', 'average', 'noplot'); %#ok<NASGU>
    image_similarity_plot(t_grp, 'mapset', 'npsplus');
    save_first_fig(fullfile(png_dir, 'fmri_data_image_similarity_plot_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 2. image_similarity_plot_bucknermaps
try
    fprintf('\n=== image_similarity_plot_bucknermaps ===\n');
    close all
    t_grp = ttest(imgs);
    stats = image_similarity_plot_bucknermaps(t_grp); %#ok<NASGU>
    save_first_fig(fullfile(png_dir, 'fmri_data_image_similarity_plot_bucknermaps_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 3. jackknife_similarity
try
    fprintf('\n=== jackknife_similarity ===\n');
    close all
    [sim_values, d, low_agreement, Nvox] = jackknife_similarity(imgs, ...
        'similarity_metric', 'correlation'); %#ok<ASGLU>
    save_first_fig(fullfile(png_dir, 'fmri_data_jackknife_similarity_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 4. montage method (fmri_data → slice montage of group t-map)
try
    fprintf('\n=== montage (fmri_data) ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc', 'k', 10);
    create_figure('m'); axis off; montage(t_grp);
    save_first_fig(fullfile(png_dir, 'fmri_data_montage_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 5. orthviews — capture frame after orthviews opens SPM viewer
try
    fprintf('\n=== orthviews ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc', 'k', 10);
    orthviews(t_grp);
    drawnow;
    pause(1);
    % SPM orthviews uses figure named 'Graphics' — find by name
    spm_fig = findobj('Type', 'figure', 'Name', 'Graphics');
    if isempty(spm_fig)
        save_first_fig(fullfile(png_dir, 'fmri_data_orthviews_sample.png'));
    else
        exportgraphics(spm_fig, fullfile(png_dir, 'fmri_data_orthviews_sample.png'), ...
                       'Resolution', 150, 'BackgroundColor', 'white');
        fprintf('Saved orthviews via Graphics window\n');
    end
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 6. mahal — Mahalanobis distance with diagnostic plot
try
    fprintf('\n=== mahal ===\n');
    close all
    [ds, expectedds, p_vals, wh_outlier_uncorr, wh_outlier_corr] = mahal(imgs); %#ok<ASGLU>
    save_first_fig(fullfile(png_dir, 'fmri_data_mahal_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 7. pca — diagnostic figure
try
    fprintf('\n=== pca ===\n');
    close all
    [scores, eig_obj, explained] = pca(imgs, 'k', 5); %#ok<ASGLU>
    save_first_fig(fullfile(png_dir, 'fmri_data_pca_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 8. rmssd_movie — first frame
try
    fprintf('\n=== rmssd_movie ===\n');
    close all
    rmssd_movie(imgs);
    save_first_fig(fullfile(png_dir, 'fmri_data_rmssd_movie_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 9. surface — render group t on cortical surfaces
try
    fprintf('\n=== surface ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc', 'k', 10);
    create_figure('s'); surface(t_grp);
    save_first_fig(fullfile(png_dir, 'fmri_data_surface_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 10. wedge_plot_by_atlas
try
    fprintf('\n=== wedge_plot_by_atlas ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc');
    [hh, vals] = wedge_plot_by_atlas(t_grp, 'atlases', {'buckner_networks'}); %#ok<ASGLU>
    save_first_fig(fullfile(png_dir, 'fmri_data_wedge_plot_by_atlas_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% =================================================================
% STATISTIC_IMAGE METHODS
% =================================================================

% 11. riverplot
try
    fprintf('\n=== riverplot (statistic_image) ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc');
    riverplot(t_grp);
    save_first_fig(fullfile(png_dir, 'statistic_image_riverplot_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 12. multi_threshold
try
    fprintf('\n=== multi_threshold (statistic_image) ===\n');
    close all
    t_grp = ttest(imgs);
    [o2, sig, poscl, negcl] = multi_threshold(t_grp); %#ok<ASGLU>
    save_first_fig(fullfile(png_dir, 'statistic_image_multi_threshold_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 13. statistic_image.table — text screenshot
try
    fprintf('\n=== statistic_image.table (text) ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc', 'k', 10);
    txt = evalc('table(t_grp);');
    save_text_as_png(txt, fullfile(png_dir, 'statistic_image_table_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 14. statistic_image.threshold — show before/after montage
try
    fprintf('\n=== statistic_image.threshold ===\n');
    close all
    t_grp = ttest(imgs);
    t_grp = threshold(t_grp, .005, 'unc', 'k', 10);
    create_figure('thr'); axis off; montage(t_grp);
    save_first_fig(fullfile(png_dir, 'statistic_image_threshold_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% =================================================================
% REGION METHODS (using a region from the group t-test)
% =================================================================

% Build a region for region examples
t_grp_thr = ttest(imgs);
t_grp_thr = threshold(t_grp_thr, .005, 'unc', 'k', 10);
r = region(t_grp_thr);
fprintf('\n%d regions for region examples\n', numel(r));

% 15. region.table — text screenshot
try
    fprintf('\n=== region.table (text) ===\n');
    close all
    txt = evalc('table(r);');
    save_text_as_png(txt, fullfile(png_dir, 'region_table_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 16. region.surface
try
    fprintf('\n=== region.surface ===\n');
    close all
    create_figure('rs'); surface(r);
    save_first_fig(fullfile(png_dir, 'region_surface_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 17. region.labelled_surface
try
    fprintf('\n=== region.labelled_surface ===\n');
    close all
    create_figure('rls'); labelled_surface(r);
    save_first_fig(fullfile(png_dir, 'region_labelled_surface_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 18. region.isosurface
try
    fprintf('\n=== region.isosurface ===\n');
    close all
    create_figure('ri'); set(gcf, 'Position', [100 100 900 700]);
    isosurface(r);
    save_first_fig(fullfile(png_dir, 'region_isosurface_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 19. region.montage
try
    fprintf('\n=== region.montage ===\n');
    close all
    montage(r, 'regioncenters', 'colormap');
    save_first_fig(fullfile(png_dir, 'region_montage_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 20. region.table_of_atlas_regions_covered — text screenshot
try
    fprintf('\n=== region.table_of_atlas_regions_covered (text) ===\n');
    close all
    txt = evalc('[~] = table_of_atlas_regions_covered(r);');
    save_text_as_png(txt, fullfile(png_dir, 'region_table_of_atlas_regions_covered_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% =================================================================
% OTHER FUNCTIONS
% =================================================================

% 21. addbrain
try
    fprintf('\n=== addbrain ===\n');
    close all
    create_figure('ab'); set(gcf, 'Position', [100 100 900 700]);
    p = addbrain('hires left'); set(p, 'FaceAlpha', .8); %#ok<NASGU>
    addbrain('thalamus'); addbrain('bg');
    view(135, 10); lightRestoreSingle;
    save_first_fig(fullfile(png_dir, 'addbrain_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 22. barplot_columns — small group of columns with errors
try
    fprintf('\n=== barplot_columns ===\n');
    close all
    rng(7);
    Y = [randn(20,1)+1, randn(20,1)+0.3, randn(20,1)-0.5];
    barplot_columns(Y, 'nofig', 'names', {'Cond A', 'Cond B', 'Cond C'});
    save_first_fig(fullfile(png_dir, 'barplot_columns_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 23. cluster_surf — Yeo network on surface
try
    fprintf('\n=== cluster_surf ===\n');
    close all
    atl = load_atlas('yeo17networks');
    cl = atlas2region(atl);
    create_figure('cs'); set(gcf, 'Position', [100 100 1000 700]);
    p = addbrain('hires left'); set(p, 'FaceAlpha', .8); %#ok<NASGU>
    cluster_surf(cl, 5, gca);
    view(-90, 0); lightRestoreSingle;
    save_first_fig(fullfile(png_dir, 'cluster_surf_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 24. image_scatterplot
try
    fprintf('\n=== image_scatterplot ===\n');
    close all
    t_grp = ttest(imgs);
    image_scatterplot(t_grp, t_grp, 'colorpoints');
    save_first_fig(fullfile(png_dir, 'image_scatterplot_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 25. canlab_results_fmridisplay
try
    fprintf('\n=== canlab_results_fmridisplay ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc', 'k', 10);
    o2 = canlab_results_fmridisplay(t_grp); %#ok<NASGU>
    save_first_fig(fullfile(png_dir, 'canlab_results_fmridisplay_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% 26. plot_correlation_matrix
try
    fprintf('\n=== plot_correlation_matrix ===\n');
    close all
    rng(7);
    S = toeplitz([1 .6 .3 .1 0 0]);
    X = mvnrnd([0 0 0 0 0 0], S, 50);
    var_names = {'A' 'B' 'C' 'D' 'E' 'F'};
    OUT = plot_correlation_matrix(X, 'var_names', var_names); %#ok<NASGU>
    save_first_fig(fullfile(png_dir, 'plot_correlation_matrix_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

close all
fprintf('\nPhase 4b runner complete.\n');
