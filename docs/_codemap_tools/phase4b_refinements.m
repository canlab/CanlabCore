% Phase 4b refinements: re-render 5 PNGs with user-requested adjustments.
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

% Render multi-line text as a tight PNG; height auto-fits the line count
% so there is no empty whitespace below the final line.
function save_text_as_png_tight(txt, out_path)
    txt = regexprep(txt, '<[^>]*>', '');     % strip HTML markup
    txt = regexprep(txt, '\x{0008}', '');     % backspace
    lines = splitlines(txt);
    while ~isempty(lines) && isempty(strtrim(lines{end})), lines(end) = []; end
    while ~isempty(lines) && isempty(strtrim(lines{1})), lines(1) = []; end
    n = numel(lines);
    if n == 0
        warning('empty text for %s', out_path); return
    end
    line_h = 16;                          % pixels per line at FontSize 9
    top_pad = 8;
    bot_pad = 8;
    fig_h = top_pad + bot_pad + line_h * n;
    fh = figure('Color', 'white', 'Units', 'pixels', ...
                'Position', [50 50 1300 fig_h], 'MenuBar', 'none', 'ToolBar', 'none');
    ax = axes('Parent', fh, 'Position', [0 0 1 1], 'Visible', 'off', ...
              'XLim', [0 1], 'YLim', [0 1]);
    % Place text in normalized coords with a small left/top inset.
    text(ax, 0.006, 1 - (top_pad/fig_h), strjoin(lines, newline), ...
         'FontName', 'Menlo', 'FontSize', 9, ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'Interpreter', 'none', 'Units', 'normalized');
    drawnow;
    exportgraphics(fh, out_path, 'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('Saved %s\n', out_path);
end

imgs = load_image_set('emotionreg');

% --- 1. barplot_columns with seaborn palette ---
try
    fprintf('\n=== REFINE: barplot_columns (with seaborn colors) ===\n');
    close all
    rng(7);
    Y = [randn(20,1)+1, randn(20,1)+0.3, randn(20,1)-0.5];
    [colors, ~] = seaborn_colors(8);
    barplot_columns(Y, 'nofig', 'names', {'Cond A', 'Cond B', 'Cond C'}, ...
                    'colors', colors(1:3));
    save_first_fig(fullfile(png_dir, 'barplot_columns_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% --- 2. cluster_surf with a thresholded t-map ---
try
    fprintf('\n=== REFINE: cluster_surf (thresholded t-map) ===\n');
    close all
    t = ttest(imgs);
    t = threshold(t, .005, 'unc', 'k', 10);
    r = region(t);
    create_figure('cs'); set(gcf, 'Position', [100 100 1000 700]);
    cluster_surf(r, 5, 'left');
    save_first_fig(fullfile(png_dir, 'cluster_surf_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% --- 3. image_scatterplot: OLS t vs robust t ---
try
    fprintf('\n=== REFINE: image_scatterplot (OLS vs robust t) ===\n');
    close all
    obj = imgs;
    % Mean-centered Reappraisal_Success regressor + intercept
    if ismember('Reappraisal_Success', obj.metadata_table.Properties.VariableNames)
        x = obj.metadata_table.Reappraisal_Success;
    else
        x = randn(size(obj.dat, 2), 1);   % fallback if metadata missing
    end
    obj.X = x - mean(x);
    obj.X(:, end+1) = 1;                  % intercept
    out_ols  = regress(obj, 'noverbose', 'nodisplay');
    out_rob  = regress(obj, 'robust', 'noverbose', 'nodisplay');
    t_ols    = get_wh_image(out_ols.t, 1);   % first regressor t
    t_rob    = get_wh_image(out_rob.t, 1);
    image_scatterplot(t_ols, t_rob, 'colorpoints');
    xlabel('OLS t-value'); ylabel('Robust regression t-value');
    save_first_fig(fullfile(png_dir, 'image_scatterplot_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% --- 4. region.table — full content, no truncation, tight crop ---
try
    fprintf('\n=== REFINE: region.table (no truncation) ===\n');
    close all
    t = ttest(imgs);
    t = threshold(t, .005, 'unc', 'k', 10);
    r = region(t);
    txt = evalc('table(r);');
    save_text_as_png_tight(txt, fullfile(png_dir, 'region_table_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

% --- 5. statistic_image.table — full content, no truncation, tight crop ---
try
    fprintf('\n=== REFINE: statistic_image.table (no truncation) ===\n');
    close all
    t = ttest(imgs);
    t = threshold(t, .005, 'unc', 'k', 10);
    txt = evalc('table(t);');
    save_text_as_png_tight(txt, fullfile(png_dir, 'statistic_image_table_sample.png'));
catch ME, fprintf(2, 'FAIL: %s\n', ME.message); end

close all
fprintf('\nRefinements complete.\n');
