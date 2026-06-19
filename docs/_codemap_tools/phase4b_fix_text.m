parent = '/Users/f003vz1/Documents/GitHub';
addpath(genpath(fullfile(parent, 'CanlabCore')));
addpath(genpath(fullfile(parent, 'Neuroimaging_Pattern_Masks')));
png_dir = fullfile(parent, 'CanlabCore', 'docs', 'class_method_pngs');

function save_text_as_png(txt, out_path)
    % Strip ALL HTML/markup tags MATLAB inserts in formatted output.
    txt = regexprep(txt, '<[^>]*>', '');
    txt = regexprep(txt, '\x{0008}', '');
    % Truncate / wrap if super long
    lines = splitlines(txt);
    % Drop blank lines at end and at start
    while ~isempty(lines) && isempty(strtrim(lines{end})), lines(end) = []; end
    while ~isempty(lines) && isempty(strtrim(lines{1})), lines(1) = []; end
    if numel(lines) > 65, lines = lines(1:65); lines{end+1} = '... (truncated)'; end
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
t_grp = ttest(imgs); t_grp = threshold(t_grp, .005, 'unc', 'k', 10);
r = region(t_grp);

% statistic_image.table
try
    fprintf('\n=== statistic_image.table ===\n');
    close all
    txt = evalc('table(t_grp);');
    save_text_as_png(txt, fullfile(png_dir, 'statistic_image_table_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% region.table
try
    fprintf('\n=== region.table ===\n');
    close all
    txt = evalc('table(r);');
    save_text_as_png(txt, fullfile(png_dir, 'region_table_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% region.table_of_atlas_regions_covered (use @image_vector path; @region is broken)
try
    fprintf('\n=== region.table_of_atlas_regions_covered ===\n');
    close all
    txt = evalc('[~] = table_of_atlas_regions_covered(t_grp);');
    save_text_as_png(txt, fullfile(png_dir, 'region_table_of_atlas_regions_covered_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

close all
fprintf('\nDone.\n');
