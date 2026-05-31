function fig = plot_hrf_atlas_curves(input_table, varargin)
%PLOT_HRF_ATLAS_CURVES Group HRF curves for the top activating atlas regions.
%
% Reads atlas region-mean columns from the per-(subject, run, model) score
% CSVs referenced by an input_table, pools per-condition per-lag values
% across subjects, ranks regions by activation magnitude, and plots the
% top N as a grid of per-region subplots. Each subplot shows the
% across-subject mean HRF curve per condition with a SEM band.
%
% Usage
% -----
%   fig = plot_hrf_atlas_curves(input_table_lf);                % defaults
%   fig = plot_hrf_atlas_curves(input_table_lf, ...
%             'Model', 'sfir', ...
%             'Conditions', {'nback-stimblock_ttl_1', 'rest_stim_ttl_1'}, ...
%             'AtlasName', 'canlab2024', ...
%             'TopN', 16, ...
%             'RankBy', 'peak_abs');
%
% Name-value parameters
% ---------------------
%   'Model'       - which row's model column to use. Default 'sfir'.
%   'Object'      - 'beta' (default) or 't'.
%   'Conditions'  - cellstr; subset of conditions to plot. Default [] (all).
%   'AtlasObj'    - the atlas object (recommended). When provided, its
%                   .labels drive column matching and region naming so
%                   multi-token atlas names like 'CANLab2024_...' work
%                   without further configuration. Default [].
%   'AtlasName'   - case-insensitive substring of the atlas-name token
%                   embedded in 'atlas_<name>_<region>_<suffix>'. Used
%                   only when AtlasObj is not provided. Default '' (any
%                   atlas; takes every atlas_*_<suffix> column).
%   'Normalize'   - which suffix to read: 'mean' (default), 'l1', 'sum'.
%                   Must match what hrf_score_wholebrain_input_table wrote.
%   'TopN'        - number of regions to plot. Default 16.
%   'RankBy'      - 'peak_abs' (default) -> max |mean curve value| over
%                                            all (condition, lag) pairs
%                   'peak'       -> max signed value (positive activations)
%                   'auc_abs'    -> max sum(|mean curve|) across lags
%                   'sd'         -> max across-lag SD (any deflection)
%   'Regions'     - cellstr; explicit region list (overrides TopN/RankBy).
%   'Layout'      - [nrows ncols]; default auto-grid from TopN.
%   'FigureSize'  - [w h] in pixels. Default scales with grid.
%   'ErrorBand'   - 'sem' (default), 'sd', or 'none'.
%   'YZero'       - true (default) to draw horizontal 0 line per subplot.
%   'Title'       - char/string figure title. Default auto-generated.
%
% Returns
% -------
%   fig - figure handle.
%
% See also: hrf_curve_summaries, hrf_curve_summary_groupstats,
%           hrf_score_wholebrain_input_table.

p = inputParser;
p.addRequired('input_table', @istable);
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('Object', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('Conditions', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('AtlasObj', [], @(x) isempty(x) || isa(x, 'atlas') || isa(x, 'image_vector'));
p.addParameter('AtlasName', '', @(x) ischar(x) || isstring(x));
p.addParameter('Normalize', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('TopN', 16, @(x) isscalar(x) && x >= 1);
p.addParameter('RankBy', 'peak_abs', @(x) ischar(x) || isstring(x));
p.addParameter('Regions', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Layout', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
p.addParameter('FigureSize', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
p.addParameter('ErrorBand', 'sem', @(x) ischar(x) || isstring(x));
p.addParameter('YZero', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Title', '', @(x) ischar(x) || isstring(x));
p.addParameter('Verbose', false, @(x) islogical(x) || isnumeric(x));
p.parse(input_table, varargin{:});
opts = p.Results;

model = lower(char(opts.Model));
object = lower(char(opts.Object));
suffix = local_normalize_suffix(opts.Normalize);

% 1. Walk the input_table, load the matching score CSV per row, extract the
%    atlas columns. Output is a long table per-row.
long = local_collect_atlas_long(input_table, model, object, opts.AtlasObj, char(opts.AtlasName), suffix, logical(opts.Verbose));
if isempty(long) || height(long) == 0
    error('plot_hrf_atlas_curves:NoData', ...
        ['No atlas columns ending in _%s found in any score CSV for ' ...
         'model=%s, object=%s. Pass ''AtlasObj'' (preferred) or ' ...
         '''AtlasName'' substring to disambiguate, or check that ' ...
         'scoring actually wrote atlas columns for this model/object.'], ...
        suffix, model, object);
end

if ~isempty(opts.Conditions)
    requested = cellstr(string(opts.Conditions));
    long = long(ismember(long.condition, string(requested)), :);
end
if height(long) == 0
    error('plot_hrf_atlas_curves:NoMatchingConditions', ...
        'After Conditions filter, no rows remain.');
end

% 2. Pool across (subject, run) per (condition, region, lag).
pooled = local_pool_subjects(long);

% 3. Pick regions.
if ~isempty(opts.Regions)
    top_regions = cellstr(string(opts.Regions));
    top_regions = intersect(top_regions, cellstr(unique(pooled.region)), 'stable');
else
    top_regions = local_rank_regions(pooled, char(opts.RankBy), opts.TopN);
end
if isempty(top_regions)
    error('plot_hrf_atlas_curves:NoRegions', 'No regions to plot.');
end

% 4. Plot.
fig = local_plot_grid(pooled, top_regions, opts);
end


% =========================================================================
% Data collection
% =========================================================================
function long = local_collect_atlas_long(input_table, model, object, atlas_obj, atlas_name, suffix, verbose)
% Returns one row per (subject, run_label, condition, region, lag_seconds, value).
% Column matching is 3-tier:
%   (a) AtlasObj.labels  -> exact end-match per label (region names from labels)
%   (b) AtlasName        -> case-insensitive substring of column name (region names
%                           derived by stripping the matched substring from the middle)
%   (c) neither          -> any atlas_*_<suffix> column (region name = entire middle)
long = table();
v = input_table.Properties.VariableNames;
file_col = sprintf('%s_scores_file', object);
if ~any(strcmp(file_col, v))
    error('plot_hrf_atlas_curves:NoFileColumn', ...
        'input_table is missing %s.', file_col);
end

suffix_str = ['_' suffix];
labels = local_labels_from_atlas(atlas_obj);

if verbose
    fprintf('plot_hrf_atlas_curves: input_table has %d rows; filtering on model=''%s'', object=''%s''\n', ...
        height(input_table), model, object);
    if ~isempty(labels)
        fprintf('  matcher: AtlasObj.labels (%d labels)\n', numel(labels));
    elseif ~isempty(atlas_name)
        fprintf('  matcher: AtlasName substring ''%s''\n', atlas_name);
    else
        fprintf('  matcher: auto (any atlas_*_%s)\n', suffix);
    end
end

chunks = {};
n_loaded = 0;
n_skipped_model = 0;
n_skipped_path = 0;
n_skipped_meta = 0;
n_skipped_no_atlas = 0;
seen_paths = strings(0, 1);
for i = 1:height(input_table)
    row_model = local_get_string(input_table, i, 'model');
    if ~isempty(model) && ~strcmpi(row_model, model)
        n_skipped_model = n_skipped_model + 1;
        continue
    end
    path = char(string(input_table.(file_col)(i)));
    if isempty(path) || exist(path, 'file') ~= 2
        n_skipped_path = n_skipped_path + 1;
        continue
    end
    try
        T = readtable(path, 'TextType', 'string');
    catch
        n_skipped_path = n_skipped_path + 1;
        continue
    end
    cols = T.Properties.VariableNames;
    if ~any(strcmp('condition', cols)) || ~any(strcmp('lag_seconds', cols))
        n_skipped_meta = n_skipped_meta + 1;
        continue
    end

    [region_cols, region_labels] = local_match_atlas_columns(cols, labels, atlas_name, suffix_str);
    if isempty(region_cols)
        n_skipped_no_atlas = n_skipped_no_atlas + 1;
        if verbose && n_skipped_no_atlas <= 2
            atlas_like = cols(startsWith(cols, 'atlas_'));
            fprintf('  row %d (model=%s): no matching atlas columns. atlas_* present (first 3): %s\n', ...
                i, row_model, strjoin(atlas_like(1:min(3, numel(atlas_like))), ', '));
        end
        continue
    end

    n_loaded = n_loaded + 1;
    seen_paths(end + 1, 1) = string(path); %#ok<AGROW>
    if verbose && n_loaded <= 3
        fprintf('  loaded row %d  model=%s  subj=%s  cols matched: %d  path: %s\n', ...
            i, row_model, local_get_string(input_table, i, 'subject'), numel(region_cols), path);
    end

    n = height(T);
    subj = string(local_get_string(input_table, i, 'subject'));
    run  = string(local_get_string(input_table, i, 'run_label'));

    for c = 1:numel(region_cols)
        chunk = table();
        chunk.subject = repmat(subj, n, 1);
        chunk.run_label = repmat(run, n, 1);
        chunk.condition = string(T.condition);
        chunk.lag_seconds = double(T.lag_seconds);
        chunk.region = repmat(string(region_labels{c}), n, 1);
        chunk.value = double(T.(region_cols{c}));
        chunks{end + 1} = chunk; %#ok<AGROW>
    end
end

if ~isempty(chunks)
    long = vertcat(chunks{:});
end

if verbose
    fprintf('  loaded %d row(s); skipped model=%d, path=%d, meta=%d, no_atlas=%d\n', ...
        n_loaded, n_skipped_model, n_skipped_path, n_skipped_meta, n_skipped_no_atlas);
    n_unique_paths = numel(unique(seen_paths));
    if n_loaded > 0 && n_unique_paths < n_loaded
        fprintf(['  WARNING: %d loaded rows reference only %d unique file paths -- ' ...
                 'multiple rows point at the same score CSV. The scoring step probably ' ...
                 'wrote one file per (subject, run) regardless of model, so every model ' ...
                 'filter ends up reading the same data.\n'], ...
                n_loaded, n_unique_paths);
    end
    if ~isempty(long)
        fprintf('  long table: %d rows; %d unique (subject, run, condition, region, lag)\n', ...
            height(long), height(unique(long(:, {'subject','run_label','condition','region','lag_seconds'}))));
        fprintf('  unique regions: %d; unique conditions: %d; lag count per (subj, run, region, cond): %d\n', ...
            numel(unique(long.region)), numel(unique(long.condition)), ...
            mode(splitapply(@numel, long.lag_seconds, findgroups(long.subject, long.run_label, long.region, long.condition))));
    end
end
end


function labels = local_labels_from_atlas(atlas_obj)
labels = {};
if isempty(atlas_obj), return; end
if isprop(atlas_obj, 'labels') && ~isempty(atlas_obj.labels)
    labels = cellstr(string(atlas_obj.labels));
    labels = labels(:);
end
end


function [matched_cols, matched_labels] = local_match_atlas_columns(cols, labels, atlas_name, suffix_str)
matched_cols = {};
matched_labels = {};

% Tier A: AtlasObj.labels drive matching. For each label, look for a column
% that ends with _<makeValidName(label)>_<suffix>. Region names come from
% the original labels (preserving any non-validname chars like '/').
if ~isempty(labels)
    for L = 1:numel(labels)
        tok = matlab.lang.makeValidName(char(labels{L}));
        end_pat = ['_' tok suffix_str];
        is_match = startsWith(cols, 'atlas_') & endsWith(cols, end_pat);
        idx = find(is_match);
        for k = 1:numel(idx)
            matched_cols{end + 1} = cols{idx(k)}; %#ok<AGROW>
            matched_labels{end + 1} = labels{L}; %#ok<AGROW>
        end
    end
    return
end

% Tier B: AtlasName substring match (case-insensitive). Region name is
% the column middle with the matched substring stripped from the front.
if ~isempty(atlas_name)
    is_atlas = startsWith(cols, 'atlas_') & endsWith(cols, suffix_str);
    is_match = is_atlas & contains(lower(cols), lower(atlas_name));
    hit_cols = cols(is_match);
    for k = 1:numel(hit_cols)
        col = hit_cols{k};
        mid = col(length('atlas_') + 1 : end - length(suffix_str));
        lc_mid = lower(mid);
        lc_name = lower(atlas_name);
        idx = strfind(lc_mid, lc_name);
        if ~isempty(idx)
            % Strip the matched substring plus a trailing underscore.
            start_after = idx(1) + length(atlas_name);
            if start_after <= length(mid) && mid(start_after) == '_'
                start_after = start_after + 1;
            end
            region = mid(start_after:end);
            if isempty(region), region = mid; end
        else
            region = mid;
        end
        matched_cols{end + 1} = col; %#ok<AGROW>
        matched_labels{end + 1} = region; %#ok<AGROW>
    end
    return
end

% Tier C: any atlas_*_<suffix> column. Region name = the full middle.
is_match = startsWith(cols, 'atlas_') & endsWith(cols, suffix_str);
hit_cols = cols(is_match);
for k = 1:numel(hit_cols)
    col = hit_cols{k};
    matched_cols{end + 1} = col; %#ok<AGROW>
    matched_labels{end + 1} = col(length('atlas_') + 1 : end - length(suffix_str)); %#ok<AGROW>
end
end


function pooled = local_pool_subjects(long)
% Collapse multiple runs per subject first, then pool across subjects.
% Uses findgroups + splitapply (function handles) rather than groupsummary
% (which only accepts a fixed set of method-name strings, not 'numel').
% Output: condition, region, lag_seconds, mean, sem, sd, n.

% Per (subject, condition, region, lag): average across runs.
[G1, subj_u, cond_u, region_u, lag_u] = findgroups( ...
    long.subject, long.condition, long.region, long.lag_seconds);
v_per_run = splitapply(@(x) mean(x, 'omitnan'), long.value, G1);
g1 = table(subj_u, cond_u, region_u, lag_u, v_per_run, ...
    'VariableNames', {'subject','condition','region','lag_seconds','value'});

% Per (condition, region, lag): mean + SD + N (finite) across subjects.
[G2, cond_u2, region_u2, lag_u2] = findgroups( ...
    g1.condition, g1.region, g1.lag_seconds);
m  = splitapply(@(x) mean(x, 'omitnan'), g1.value, G2);
sd = splitapply(@(x) std(x, 'omitnan'),  g1.value, G2);
nn = splitapply(@(x) sum(isfinite(x)),   g1.value, G2);

pooled = table(cond_u2, region_u2, lag_u2, m, sd, nn, sd ./ sqrt(max(nn, 1)), ...
    'VariableNames', {'condition','region','lag_seconds','mean','sd','n','sem'});
end


function top_regions = local_rank_regions(pooled, rank_by, top_n)
% Rank by a per-region scalar collapsed across (condition, lag), then keep
% the top top_n region names in descending order.
regions = unique(pooled.region, 'stable');
scores = zeros(numel(regions), 1);

for r = 1:numel(regions)
    mask = pooled.region == regions(r);
    vals = pooled.mean(mask);
    vals = vals(isfinite(vals));
    if isempty(vals), scores(r) = -Inf; continue; end
    switch lower(rank_by)
        case 'peak_abs'
            scores(r) = max(abs(vals));
        case 'peak'
            scores(r) = max(vals);
        case 'auc_abs'
            scores(r) = sum(abs(vals));
        case 'sd'
            scores(r) = std(vals);
        otherwise
            error('plot_hrf_atlas_curves:UnknownRankBy', ...
                'Unknown RankBy: %s. Use peak_abs, peak, auc_abs, or sd.', rank_by);
    end
end

[~, order] = sort(scores, 'descend');
keep = order(1:min(top_n, numel(order)));
top_regions = cellstr(regions(keep));
end


% =========================================================================
% Plot grid
% =========================================================================
function fig = local_plot_grid(pooled, top_regions, opts)
n = numel(top_regions);
if ~isempty(opts.Layout) && numel(opts.Layout) == 2
    nrows = opts.Layout(1); ncols = opts.Layout(2);
else
    [nrows, ncols] = local_auto_grid(n);
end

if isempty(opts.FigureSize)
    fig_w = 280 * ncols;
    fig_h = 220 * nrows + 60;
else
    fig_w = opts.FigureSize(1);
    fig_h = opts.FigureSize(2);
end

fig = figure('Color', 'w', 'Position', [80, 80, fig_w, fig_h]);
tl = tiledlayout(fig, nrows, ncols, 'Padding', 'compact', 'TileSpacing', 'compact');

cond_list = cellstr(unique(pooled.condition, 'stable'));
colors = local_color_palette(numel(cond_list));

% Y-axis: share scale across panels so amplitudes are visually comparable.
y_min = Inf; y_max = -Inf;
for r = 1:n
    mask_r = pooled.region == string(top_regions{r});
    band_lo = pooled.mean(mask_r) - local_err(pooled, mask_r, opts.ErrorBand);
    band_hi = pooled.mean(mask_r) + local_err(pooled, mask_r, opts.ErrorBand);
    band_lo = band_lo(isfinite(band_lo));
    band_hi = band_hi(isfinite(band_hi));
    if ~isempty(band_lo), y_min = min(y_min, min(band_lo)); end
    if ~isempty(band_hi), y_max = max(y_max, max(band_hi)); end
end
if ~isfinite(y_min) || ~isfinite(y_max)
    y_min = -1; y_max = 1;
end
y_pad = 0.05 * (y_max - y_min + eps);

for r = 1:n
    ax = nexttile(tl);
    hold(ax, 'on');
    if logical(opts.YZero)
        yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.75, 'HandleVisibility', 'off');
    end

    for c = 1:numel(cond_list)
        cond = string(cond_list{c});
        mask = pooled.region == string(top_regions{r}) & pooled.condition == cond;
        if ~any(mask), continue; end
        sub = sortrows(pooled(mask, :), 'lag_seconds');
        x = sub.lag_seconds(:);
        y = sub.mean(:);
        e = local_err(pooled, mask, opts.ErrorBand);
        if any(e > 0)
            xx = [x; flipud(x)];
            yy = [y + e; flipud(y - e)];
            fill(ax, xx, yy, colors(c, :), 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
        plot(ax, x, y, '-', 'Color', colors(c, :), 'LineWidth', 1.6, ...
            'DisplayName', char(cond));
    end

    title(ax, char(top_regions{r}), 'Interpreter', 'none', 'FontWeight', 'normal');
    xlabel(ax, 'Lag (s)'); ylabel(ax, 'Mean signal');
    grid(ax, 'on'); box(ax, 'on');
    ylim(ax, [y_min - y_pad, y_max + y_pad]);

    if r == 1
        legend(ax, 'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
    end
end

if isempty(char(opts.Title))
    tstr = sprintf('Top %d %s regions  (model=%s, object=%s, %s ± %s)', ...
        n, char(opts.AtlasName), char(opts.Model), char(opts.Object), ...
        char(opts.Normalize), char(opts.ErrorBand));
else
    tstr = char(opts.Title);
end
title(tl, tstr, 'Interpreter', 'none');
end


function e = local_err(pooled, mask, band)
switch lower(char(band))
    case 'sem'
        e = pooled.sem(mask);
    case 'sd'
        e = pooled.sd(mask);
    case 'none'
        e = zeros(sum(mask), 1);
    otherwise
        e = pooled.sem(mask);
end
e(~isfinite(e)) = 0;
e = e(:);
% match the sort the caller will apply
[~, ord] = sort(pooled.lag_seconds(mask));
e = e(ord);
end


% =========================================================================
% Utilities
% =========================================================================
function [nrows, ncols] = local_auto_grid(n)
ncols = max(1, ceil(sqrt(n)));
nrows = ceil(n / ncols);
end


function s = local_normalize_suffix(mode)
switch lower(strtrim(char(mode)))
    case 'mean', s = 'mean';
    case 'l1',   s = 'meanL1';
    case 'none', s = 'sum';
    otherwise,   s = char(mode);
end
end


function s = local_get_string(t, i, col)
s = '';
if any(strcmp(col, t.Properties.VariableNames))
    val = t.(col)(i);
    if isstring(val) || ischar(val)
        s = char(val);
    elseif iscell(val)
        s = char(val{1});
    else
        try
            s = char(string(val));
        catch
        end
    end
end
end


function colors = local_color_palette(k)
% A reasonable diverging palette for up to ~8 conditions.
base = [
    0.215 0.494 0.722;  % blue
    0.894 0.102 0.110;  % red
    0.302 0.686 0.290;  % green
    0.596 0.306 0.639;  % purple
    1.000 0.498 0.000;  % orange
    1.000 1.000 0.200;  % yellow
    0.651 0.337 0.157;  % brown
    0.969 0.506 0.749]; % pink
if k <= size(base, 1)
    colors = base(1:k, :);
else
    colors = lines(k);
end
end
