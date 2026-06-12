% Phase 4a runner: produce sample PNGs for 4 documentation pages.

parent = '/Users/f003vz1/Documents/GitHub';
addpath(genpath(fullfile(parent, 'CanlabCore')));
addpath(genpath(fullfile(parent, 'Neuroimaging_Pattern_Masks')));

png_dir = fullfile(parent, 'CanlabCore', 'docs', 'class_method_pngs');
if ~exist(png_dir, 'dir'), mkdir(png_dir); end

function save_first_fig(out_path, varargin)
    % Optional name-value 'min_height', H to ensure figure is at least H pixels
    % tall before exporting (avoids clipping of xlabels in shallow figures).
    p = inputParser;
    p.addParameter('min_height', 0);
    p.addParameter('min_width', 0);
    p.parse(varargin{:});
    min_h = p.Results.min_height;
    min_w = p.Results.min_width;

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
    if min_h > 0 && pos(4) < min_h
        pos(4) = min_h;
    end
    if min_w > 0 && pos(3) < min_w
        pos(3) = min_w;
    end
    first_fig.Position = pos;
    drawnow;
    exportgraphics(first_fig, out_path, 'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('Saved %s\n', out_path);
end

% =============================================================
% 1. atlas_methods.md  → atlas_methods_sample.png
% =============================================================
try
    fprintf('\n=== Example 1: atlas_methods (canlab2024 isosurface + montage) ===\n');
    close all
    obj = load_atlas('canlab2024');
    create_figure('fig');
    set(gcf, 'Position', [100 100 1100 850]);
    isosurface(obj);
    view(135, 30); lightFollowView;
    create_figure('fig2'); axis off; montage(obj);
    save_first_fig(fullfile(png_dir, 'atlas_methods_sample.png'));
catch ME
    fprintf(2, 'Example 1 FAILED: %s\n', ME.message);
end

% =============================================================
% 2. atlas_select_atlas_subset.md
% =============================================================
try
    fprintf('\n=== Example 2: atlas_select_atlas_subset ===\n');
    close all
    obj = load_atlas('canlab2024');
    thal = select_atlas_subset(obj, {'Thal'});
    create_figure('fig');
    set(gcf, 'Position', [100 100 900 700]);
    isosurface(thal);
    save_first_fig(fullfile(png_dir, 'atlas_select_atlas_subset_sample.png'));
catch ME
    fprintf(2, 'Example 2 FAILED: %s\n', ME.message);
end

% =============================================================
% 3. fmri_data_annotate_binary_results_map.md
% =============================================================
try
    fprintf('\n=== Example 3: annotate_binary_results_map ===\n');
    close all
    obj = load_image_set('emotionreg');
    t = ttest(obj, .005, 'uncorrected');
    t.dat = single(t.dat > 3);
    t = fmri_data(t);
    RESULTS = annotate_binary_results_map(t); %#ok<NASGU>
    % Figure 1 from this function is intentionally short (940x206). Both the
    % xlabel and the right-side "Transmodal ->" annotation overflow the
    % figure border. Pad width AND height, and shift / shrink the axes so
    % everything stays inside.
    figs = findobj('type', 'figure');
    nums = arrayfun(@(f) f.Number, figs);
    [~, ix] = min(nums);
    fig1 = figs(ix);
    fig1.Position = [fig1.Position(1), fig1.Position(2), 1200, 420];
    ax = findall(fig1, 'type', 'axes');
    for a = ax'
        p = a.Position;
        % Shift up to expose xlabel below; shrink width to expose right label.
        a.Position = [p(1), p(2) + 0.13, p(3) * 0.86, p(4) * 0.65];
    end
    drawnow;
    save_first_fig(fullfile(png_dir, 'fmri_data_annotate_binary_results_map_sample.png'));
catch ME
    fprintf(2, 'Example 3 FAILED: %s\n', ME.message);
end

% =============================================================
% 4. fmri_data_outliers.md
% =============================================================
try
    fprintf('\n=== Example 4: outliers (notimeseries) ===\n');
    close all
    obj = load_image_set('emotionreg');
    [est_outliers_uncorr, est_outliers_corr, outlier_tables] = ...
        outliers(obj, 'notimeseries'); %#ok<ASGLU>
    save_first_fig(fullfile(png_dir, 'fmri_data_outliers_sample.png'));
catch ME
    fprintf(2, 'Example 4 FAILED: %s\n', ME.message);
end

close all
fprintf('\nPhase 4a runner complete.\n');
