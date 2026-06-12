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

function save_text_as_png(txt, out_path)
    txt = regexprep(txt, '<a[^>]*>', '');
    txt = regexprep(txt, '</a>', '');
    txt = regexprep(txt, '\x{0008}', '');
    lines = splitlines(txt);
    if numel(lines) > 60, lines = lines(1:60); end
    n_lines = numel(lines);
    fh = figure('Color', 'white', 'Units', 'pixels', ...
                'Position', [50 50 1300 max(280, 22 * (n_lines + 2))], ...
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

imgs = load_image_set('emotionreg');

% jackknife_similarity — pass doplot with explicit value
try
    fprintf('\n=== FIX 2: jackknife_similarity ===\n');
    close all
    [sim_values, d, low_agreement, Nvox] = jackknife_similarity(imgs, ...
        'similarity_metric', 'correlation', 'doplot', true); %#ok<ASGLU>
    save_first_fig(fullfile(png_dir, 'fmri_data_jackknife_similarity_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% region.table_of_atlas_regions_covered — function self-disclaims as broken,
% so we capture the @image_vector version's output text but show the call
% pattern most users will use after region(t).
try
    fprintf('\n=== FIX 2: region.table_of_atlas_regions_covered ===\n');
    close all
    t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc', 'k', 10);
    txt = evalc('[~] = table_of_atlas_regions_covered(t_grp);');
    save_text_as_png(txt, fullfile(png_dir, 'region_table_of_atlas_regions_covered_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

close all
fprintf('\nDone.\n');
