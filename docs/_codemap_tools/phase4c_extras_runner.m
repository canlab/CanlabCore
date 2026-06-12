% Phase 4c-extras runner: produce sample PNGs for the new stand-alone
% function pages (roc_plot, clusterdata_permtest). The PNG for
% canlab_force_directed_graph was rendered manually by the user.
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
    if isempty(figs), warning('No figures to save for %s', out_path); return; end
    nums = arrayfun(@(f) f.Number, figs);
    [~, ix] = min(nums);
    first_fig = figs(ix);
    figure(first_fig);
    drawnow;
    pos = first_fig.Position;
    if p.Results.min_height > 0 && pos(4) < p.Results.min_height, pos(4) = p.Results.min_height; end
    if p.Results.min_width > 0 && pos(3) < p.Results.min_width, pos(3) = p.Results.min_width; end
    first_fig.Position = pos;
    drawnow;
    exportgraphics(first_fig, out_path, 'Resolution', 150, 'BackgroundColor', 'white');
    fprintf('Saved %s\n', out_path);
end

% --- 1. roc_plot ----------------------------------------------------------
% Two overlapping Gaussian distributions, treated as continuous "pattern
% expression" values for two classes.
try
    fprintf('\n=== roc_plot ===\n');
    close all
    rng(42);
    n = 60;
    pattern_exp_neg = randn(n, 1) - 0.5;     % "negatives" (class 0)
    pattern_exp_pos = randn(n, 1) + 0.8;     % "positives" (class 1)
    input_values   = [pattern_exp_neg; pattern_exp_pos];
    binary_outcome = [false(n, 1); true(n, 1)];
    create_figure('roc'); set(gcf, 'Position', [100 100 700 600]);
    ROC = roc_plot(input_values, binary_outcome, ...
                   'color', [.2 .4 .8], 'plotmethod', 'observed', ...
                   'noboot', 'nooutput'); %#ok<NASGU>
    save_first_fig(fullfile(png_dir, 'roc_plot_sample.png'));
catch ME, fprintf(2, 'FAIL roc_plot: %s\n', ME.message); end

% --- 2. clusterdata_permtest ---------------------------------------------
% Three well-separated Gaussian blobs in 10-D + a few noise dims; the
% function should pick k = 3 as the best cluster count.
try
    fprintf('\n=== clusterdata_permtest ===\n');
    close all
    rng(7);
    n_per = 30; n_vars = 10;
    X = [randn(n_per, n_vars);
         randn(n_per, n_vars) + 2.0;
         randn(n_per, n_vars) - 2.0];
    stats = clusterdata_permtest(X, 'k', 2:6, 'doplot', true, 'verbose', false); %#ok<NASGU>
    save_first_fig(fullfile(png_dir, 'clusterdata_permtest_sample.png'));
catch ME, fprintf(2, 'FAIL clusterdata_permtest: %s\n', ME.message); end

close all
fprintf('\nPhase 4c-extras runner complete.\n');
