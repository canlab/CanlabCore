% Fixes for Phase 4b failures
parent = '/Users/f003vz1/Documents/GitHub';
addpath(genpath(fullfile(parent, 'CanlabCore')));
addpath(genpath(fullfile(parent, 'Neuroimaging_Pattern_Masks')));
png_dir = fullfile(parent, 'CanlabCore', 'docs', 'class_method_pngs');

function save_first_fig(out_path)
    figs = findobj('type', 'figure');
    if isempty(figs), warning('no fig'); return; end
    nums = arrayfun(@(f) f.Number, figs);
    [~, ix] = min(nums);
    drawnow;
    exportgraphics(figs(ix), out_path, 'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('Saved %s\n', out_path);
end

imgs = load_image_set('emotionreg');

% --- 1. image_similarity_plot — pass image set directly with 'average' ---
try
    fprintf('\n=== FIX: image_similarity_plot ===\n');
    close all
    image_similarity_plot(imgs, 'mapset', 'npsplus', 'average');
    save_first_fig(fullfile(png_dir, 'fmri_data_image_similarity_plot_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% --- 2. jackknife_similarity — needs 'doplot' flag ---
try
    fprintf('\n=== FIX: jackknife_similarity ===\n');
    close all
    [sim_values, d, low_agreement, Nvox] = jackknife_similarity(imgs, ...
        'similarity_metric', 'correlation', 'doplot'); %#ok<ASGLU>
    save_first_fig(fullfile(png_dir, 'fmri_data_jackknife_similarity_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% --- 3. wedge_plot_by_atlas — use 'yeo17networks' instead of 'buckner_networks' ---
try
    fprintf('\n=== FIX: wedge_plot_by_atlas ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc');
    [hh, vals] = wedge_plot_by_atlas(t_grp, 'atlases', {'yeo17networks'}); %#ok<ASGLU>
    save_first_fig(fullfile(png_dir, 'fmri_data_wedge_plot_by_atlas_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% --- 4. riverplot — pass two fmri_data layers ---
try
    fprintf('\n=== FIX: riverplot (statistic_image) ===\n');
    close all
    % Build two fmri_data layers from npsplus signature set & emotionreg t-map
    layer1 = load_image_set('npsplus');
    if size(layer1.dat, 2) > 4
        layer1 = get_wh_image(layer1, 1:4);
    end
    layer2 = ttest(imgs); layer2 = threshold(layer2, .005, 'unc');
    layer2 = fmri_data(layer2);
    layer1.image_names = char({'NPS','SIIPS','GenS','VPS'});
    layer2.image_names = char({'EmoReg group t'});
    riverplot(layer1, 'layer2', layer2);
    save_first_fig(fullfile(png_dir, 'statistic_image_riverplot_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% --- 5. region.table_of_atlas_regions_covered — pass an explicit atlas ---
function save_text_as_png(txt, out_path)
    txt = regexprep(txt, '<a[^>]*>', '');
    txt = regexprep(txt, '</a>', '');
    txt = regexprep(txt, '\x{0008}', '');
    lines = splitlines(txt);
    if numel(lines) > 60, lines = lines(1:60); end
    n_lines = numel(lines);
    fh = figure('Color', 'white', 'Units', 'pixels', ...
                'Position', [50 50 1200 max(280, 22 * (n_lines + 2))], ...
                'MenuBar', 'none', 'ToolBar', 'none');
    ax = axes('Parent', fh, 'Position', [0 0 1 1], 'Visible', 'off', ...
              'XLim', [0 1], 'YLim', [0 1]);
    text(ax, 0.01, 0.99, strjoin(lines, newline), ...
         'FontName', 'Menlo', 'FontSize', 9, ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'Interpreter', 'none');
    drawnow;
    exportgraphics(fh, out_path, 'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('Saved %s\n', out_path);
end

try
    fprintf('\n=== FIX: region.table_of_atlas_regions_covered ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc', 'k', 10);
    r = region(t_grp);
    atl = load_atlas('canlab2018_2mm');
    txt = evalc('[~] = table_of_atlas_regions_covered(r, atl);');
    save_text_as_png(txt, fullfile(png_dir, 'region_table_of_atlas_regions_covered_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% --- 6. cluster_surf — deprecated, use simplest invocation ---
try
    fprintf('\n=== FIX: cluster_surf ===\n');
    close all
    atl = load_atlas('yeo17networks');
    cl = atlas2region(atl);
    create_figure('cs'); set(gcf, 'Position', [100 100 1000 700]);
    cluster_surf(cl, 5, 'left');
    save_first_fig(fullfile(png_dir, 'cluster_surf_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

close all
fprintf('\nDone.\n');
