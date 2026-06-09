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
%   'Source'      - which score-column family to plot:
%                     'atlas'     (default) -> atlas_<name>_<region>_<suffix>
%                     'signature' -> sig_<set>_<name>  (e.g. sig_all_NPS)
%                     'imageset'  -> map_<set>_<name>  (e.g. map_bucknerlab_*)
%                   For 'signature'/'imageset', use 'Set' to pick a set and
%                   AtlasObj/AtlasName/Normalize are ignored. Despite the
%                   function name, all three families plot identically (one
%                   curve per unit: region / signature / network).
%   'Set'         - for Source 'signature'/'imageset': case-insensitive
%                   substring of the set token to select, e.g. 'bucknerlab'
%                   or 'all'. Default '' (every sig_*/map_* column). The
%                   curve label is the column with the source_<set>_ prefix
%                   stripped.
%   'Model'       - which row's model column to use. Default 'sfir'.
%   'Object'      - 'beta' (default) or 't'.
%   'Conditions'  - cellstr; subset of conditions to plot. Supports glob
%                   wildcards (e.g. '*_heat_start_ttl_1'). Default [] (all).
%   'CollapseConditions' - false (default) keeps each matched condition as
%                   its own curve. true relabels every condition matching
%                   a given pattern to that pattern's canonical name, so
%                   all matches pool into ONE averaged curve. E.g. with
%                   8 body-site conditions all ending in _heat_start_ttl_1,
%                   'Conditions', {'*_heat_start_ttl_1'} +
%                   'CollapseConditions', true gives a single
%                   'heat_start_ttl_1' curve (per study_label) instead of 8.
%   'ConditionLabels' - cellstr parallel to Conditions; explicit collapsed
%                   names. Default derives the name from each pattern
%                   (strip '*' and surrounding underscores).
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
% First arg accepts:
%   - a single input_table (current behavior; single-study mode)
%   - a struct array with fields .label and .table (or .input_table),
%     e.g. struct('label', {'lf_acc','lf_exp'}, 'table', {it1, it2}).
%     Each source's CSVs are tagged with the supplied label so per-region
%     panels can compare across studies/contexts where BIDS condition
%     names collide.
p.addRequired('input_table', @(x) istable(x) || local_is_source_struct(x));
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
p.addParameter('CollapseConditions', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('ConditionLabels', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('BalancedNesting', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Source', 'atlas', @(x) ischar(x) || isstring(x));
p.addParameter('Set', '', @(x) ischar(x) || isstring(x));
p.parse(input_table, varargin{:});
opts = p.Results;

model = lower(char(opts.Model));
object = lower(char(opts.Object));
suffix = local_normalize_suffix(opts.Normalize);
source_family = lower(strtrim(char(opts.Source)));
if ~ismember(source_family, {'atlas', 'signature', 'imageset'})
    error('plot_hrf_atlas_curves:UnknownSource', ...
        'Source must be ''atlas'', ''signature'', or ''imageset''. Got ''%s''.', source_family);
end

% 1. Walk the input_table(s), load the matching score CSV per row, extract
%    the requested source columns. Output is a long table per row, tagged
%    with a 'study_label' string column when multiple sources were passed.
match_spec = struct('family', source_family, 'set', char(opts.Set), ...
    'atlas_obj', opts.AtlasObj, 'atlas_name', char(opts.AtlasName), 'suffix', suffix);
sources = local_normalize_sources(input_table);
long_chunks = cell(numel(sources), 1);
for s = 1:numel(sources)
    chunk = local_collect_source_long(sources(s).table, model, object, ...
        match_spec, logical(opts.Verbose));
    if isempty(chunk) || height(chunk) == 0, continue; end
    chunk.study_label = repmat(string(sources(s).label), height(chunk), 1);
    % source_id is unique per input table (distinct from study_label, since
    % several sources can share a label). Used by balanced nesting to give
    % each source equal weight within a shared label.
    chunk.source_id = repmat(string(sprintf('src%03d', s)), height(chunk), 1);
    long_chunks{s} = chunk;
end
keep = ~cellfun(@isempty, long_chunks);
if ~any(keep)
    pfx = local_source_prefix(source_family);
    if strcmp(source_family, 'atlas')
        hint = sprintf(['Pass ''AtlasObj'' (preferred) or ''AtlasName'' substring ' ...
            'to disambiguate, or check that scoring wrote %s columns (suffix _%s) ' ...
            'for this model/object.'], pfx, suffix);
    else
        hint = sprintf(['Check that scoring wrote %s columns for this ' ...
            'model/object, and that your ''Set'' substring (''%s'') matches ' ...
            'the set token in the column names.'], pfx, char(opts.Set));
    end
    error('plot_hrf_atlas_curves:NoData', ...
        'No %s* score columns found for source=''%s'', model=%s, object=%s. %s', ...
        pfx, source_family, model, object, hint);
end
long = vertcat(long_chunks{keep});

% Preserve the original (pre-collapse) condition so balanced nesting can
% give each original condition (e.g. each body site) equal weight even
% after they're relabeled to a shared collapsed name. When not collapsing,
% orig_condition == condition (a no-op level).
long.orig_condition = long.condition;

if ~isempty(opts.Conditions)
    requested = cellstr(string(opts.Conditions));
    if logical(opts.CollapseConditions)
        % Collapse mode: every condition matching a pattern is RELABELED
        % to that pattern's canonical name, so all matches pool into one
        % averaged series instead of one series per distinct condition.
        % orig_condition keeps the body-site identity for balanced nesting.
        labels = local_resolve_condition_labels(requested, opts.ConditionLabels);
        [keep_cond, new_cond] = local_collapse_conditions(long.condition, requested, labels);
        long = long(keep_cond, :);
        long.condition = new_cond(keep_cond);
    else
        long = long(local_condition_pattern_match(long.condition, requested), :);
    end
end
if height(long) == 0
    error('plot_hrf_atlas_curves:NoMatchingConditions', ...
        'After Conditions filter, no rows remain.');
end

% 2. Pool across (subject, run) per (condition, region, lag).
pooled = local_pool_subjects(long, logical(opts.BalancedNesting));
if logical(opts.Verbose)
    if isempty(pooled) || height(pooled) == 0
        fprintf('  pooled: 0 rows  <-- this will cause "No regions to plot"\n');
    else
        fprintf('  pooled: %d rows; %d unique regions; ranges: mean=[%.3g, %.3g], n=[%d, %d]\n', ...
            height(pooled), numel(unique(pooled.region)), ...
            min(pooled.mean, [], 'omitnan'), max(pooled.mean, [], 'omitnan'), ...
            min(pooled.n), max(pooled.n));
    end
end

% 3. Pick regions.
if ~isempty(opts.Regions)
    % Allow either exact match against pooled.region (when AtlasObj was
    % passed and labels are clean) or case-insensitive substring match
    % (when AtlasName/no-arg paths leave pooled.region holding the full
    % atlas_<name>_<region> token). Substring match handles e.g. user
    % passing 'Cblm_VIIIa_L' when pooled has
    % 'CANLab2024_MNI152NLin2009cAsym_coarse_2mm_Cblm_VIIIa_L'.
    requested = cellstr(string(opts.Regions));
    pool_regions = cellstr(unique(pooled.region));
    top_regions = {};
    for r = 1:numel(requested)
        req = requested{r};
        exact = pool_regions(strcmp(pool_regions, req));
        if ~isempty(exact)
            top_regions = [top_regions; exact(:)]; %#ok<AGROW>
            continue
        end
        substr = pool_regions(contains(lower(pool_regions), lower(req)));
        if ~isempty(substr)
            top_regions = [top_regions; substr(:)]; %#ok<AGROW>
        end
    end
    top_regions = unique(top_regions, 'stable');
else
    top_regions = local_rank_regions(pooled, char(opts.RankBy), opts.TopN);
end
if isempty(top_regions)
    avail = cellstr(unique(pooled.region));
    n_show = min(30, numel(avail));
    avail_str = strjoin(avail(1:n_show), ', ');
    if numel(avail) > n_show
        avail_str = sprintf('%s, ... (%d total)', avail_str, numel(avail));
    end
    if ~isempty(opts.Regions)
        % Regions were requested but none matched -- almost always means
        % the score CSVs were scored with a DIFFERENT atlas than the one
        % whose region names are being requested. List what IS available
        % so the mismatch is obvious.
        requested_str = strjoin(cellstr(string(opts.Regions)), ', ');
        error('plot_hrf_atlas_curves:NoMatchingRegions', ...
            ['None of the requested Regions {%s} matched any atlas column ' ...
             'in the score CSVs. This usually means the CSVs were scored ' ...
             'with a different atlas than the one you are plotting -- e.g. ' ...
             'scored with canlab2024 but plotting painpathways2024 region ' ...
             'names. Re-score with the desired AtlasObj (and a distinct ' ...
             '''AtlasName'' if its atlas_name collides with another atlas).\n' ...
             'Available atlas regions in the CSVs (%d): %s'], ...
            requested_str, numel(avail), avail_str);
    else
        error('plot_hrf_atlas_curves:NoRegions', ...
            ['No regions to plot. pooled table has %d rows / %d unique ' ...
             'regions. If pooled is non-empty but ranking returned empty, ' ...
             'the RankBy metric (%s) may have produced all-NaN scores -- ' ...
             'try ''peak_abs''.\nAvailable regions: %s'], ...
            height(pooled), numel(avail), char(opts.RankBy), avail_str);
    end
end

% 4. Plot.
fig = local_plot_grid(pooled, top_regions, opts);
end


% =========================================================================
% Data collection
% =========================================================================
function long = local_collect_source_long(input_table, model, object, match_spec, verbose)
% Returns one row per (subject, run_label, condition, region, lag_seconds,
% value), where 'region' holds the unit label (atlas region / signature /
% network) regardless of source family. Column matching is dispatched by
% match_spec.family ('atlas' | 'signature' | 'imageset').
long = table();
v = input_table.Properties.VariableNames;
file_col = sprintf('%s_scores_file', object);
if ~any(strcmp(file_col, v))
    error('plot_hrf_atlas_curves:NoFileColumn', ...
        'input_table is missing %s.', file_col);
end

suffix_str = ['_' match_spec.suffix];
labels = local_labels_from_atlas(match_spec.atlas_obj);

if verbose
    fprintf('plot_hrf_atlas_curves: input_table has %d rows; source=''%s'', model=''%s'', object=''%s''\n', ...
        height(input_table), match_spec.family, model, object);
    switch match_spec.family
        case 'atlas'
            if ~isempty(labels)
                fprintf('  matcher: AtlasObj.labels (%d labels)\n', numel(labels));
            elseif ~isempty(match_spec.atlas_name)
                fprintf('  matcher: AtlasName substring ''%s''\n', match_spec.atlas_name);
            else
                fprintf('  matcher: auto (any atlas_*_%s)\n', match_spec.suffix);
            end
        case {'signature', 'imageset'}
            pfx = local_source_prefix(match_spec.family);
            if ~isempty(match_spec.set)
                fprintf('  matcher: %s* columns, Set substring ''%s''\n', pfx, match_spec.set);
            else
                fprintf('  matcher: any %s* column\n', pfx);
            end
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

    [region_cols, region_labels] = local_match_source_columns(cols, match_spec, labels, suffix_str);
    if isempty(region_cols)
        n_skipped_no_atlas = n_skipped_no_atlas + 1;
        if verbose && n_skipped_no_atlas <= 2
            pfx = local_source_prefix(match_spec.family);
            like = cols(startsWith(cols, pfx));
            fprintf('  row %d (model=%s): no matching %s columns. %s* present (first 3): %s\n', ...
                i, row_model, match_spec.family, pfx, strjoin(like(1:min(3, numel(like))), ', '));
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
        % Cheap counts only -- avoid unique() / splitapply on the full
        % (potentially 10M+) row table, which can be slow or run into
        % memory pressure and is not worth it for diagnostics.
        fprintf('  long table: %d rows; unique regions: %d; unique conditions: %d\n', ...
            height(long), numel(unique(long.region)), numel(unique(long.condition)));
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


function pfx = local_source_prefix(family)
switch family
    case 'atlas',     pfx = 'atlas_';
    case 'signature', pfx = 'sig_';
    case 'imageset',  pfx = 'map_';
    otherwise,        pfx = 'atlas_';
end
end


function [matched_cols, matched_labels] = local_match_source_columns(cols, match_spec, labels, suffix_str)
% Dispatch column matching by source family.
switch match_spec.family
    case 'atlas'
        [matched_cols, matched_labels] = local_match_atlas_columns( ...
            cols, labels, match_spec.atlas_name, suffix_str);
    case {'signature', 'imageset'}
        [matched_cols, matched_labels] = local_match_prefixed_columns( ...
            cols, local_source_prefix(match_spec.family), match_spec.set);
    otherwise
        matched_cols = {}; matched_labels = {};
end
end


function [matched_cols, matched_labels] = local_match_prefixed_columns(cols, prefix, set_name)
% Match sig_/map_ columns. Exclude *_se uncertainty columns. When set_name
% is non-empty, restrict to columns whose name contains the set token
% (case-insensitive) and strip 'prefix<set>_' to get the unit label;
% otherwise strip just 'prefix' (label keeps the set, e.g. 'all_NPS').
matched_cols = {};
matched_labels = {};
set_name = char(set_name);
for k = 1:numel(cols)
    col = cols{k};
    if ~startsWith(col, prefix), continue; end
    if endsWith(col, '_se'), continue; end
    if ~isempty(set_name) && ~contains(lower(col), lower(set_name))
        continue
    end
    mid = col(length(prefix) + 1 : end);   % drop 'sig_' / 'map_'
    if ~isempty(set_name)
        lc_mid = lower(mid);
        idx = strfind(lc_mid, lower(set_name));
        if ~isempty(idx)
            start_after = idx(1) + length(set_name);
            if start_after <= length(mid) && mid(start_after) == '_'
                start_after = start_after + 1;
            end
            label = mid(start_after:end);
            if isempty(label), label = mid; end
        else
            label = mid;
        end
    else
        label = mid;
    end
    matched_cols{end + 1} = col; %#ok<AGROW>
    matched_labels{end + 1} = label; %#ok<AGROW>
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


function pooled = local_pool_subjects(long, balanced)
% Two-stage pooling: (1) collapse everything within each subject down to a
% single value per (study_label, condition, region, lag); (2) average those
% per-subject values across subjects, with SEM = sd/sqrt(n_subjects).
%
% Stage 1 nesting:
%   balanced=true (default) -- hierarchical, equal weight at each level.
%     Average AWAY the nuisance factors one at a time, innermost first:
%       run_label  -> average runs within each (subject, source, body site)
%       orig_condition -> average body sites equally within each (subject, source)
%       source_id  -> average sources equally within each (subject, label)
%     This makes the result independent of how many runs/body-sites/source
%     rows each cell happens to have -- no assumption of balance.
%   balanced=false -- legacy row-weighted: a single flat mean over all
%     nuisance rows at once (cells with more rows get more weight).
%
% In both cases the curve-identity columns kept through to stage 2 are
% (subject, study_label[if present], condition, region, lag_seconds).
if nargin < 2, balanced = true; end

has_study = any(strcmp('study_label', long.Properties.VariableNames));

% Columns that define a per-subject curve cell (kept through stage 1).
keep_cols = {'subject', 'condition', 'region', 'lag_seconds'};
if has_study
    keep_cols = [keep_cols, {'study_label'}];
end
keep_cols = intersect(keep_cols, long.Properties.VariableNames, 'stable');

% Nuisance factors to average away, innermost -> outermost.
nuisance = {'run_label', 'orig_condition', 'source_id'};
nuisance = intersect(nuisance, long.Properties.VariableNames, 'stable');

work = long;
if balanced
    % Peel one nuisance factor at a time so each level is equal-weighted.
    for i = 1:numel(nuisance)
        f = nuisance{i};
        group_cols = setdiff(work.Properties.VariableNames, {f, 'value'}, 'stable');
        work = local_group_mean(work, group_cols);
    end
    % Any remaining nuisance columns (none expected) collapse here.
    g1 = local_group_mean(work, keep_cols);
else
    % Legacy: single flat mean over all nuisance rows at once.
    g1 = local_group_mean(work, keep_cols);
end

% Stage 2: across-subject mean + SEM. n is the count of FINITE subjects
% (not GroupCount, which would include subjects whose value is NaN at this
% lag, e.g. an all-NaN region/volume), so the SEM denominator is honest.
final_group = intersect({'study_label', 'condition', 'region', 'lag_seconds'}, ...
    g1.Properties.VariableNames, 'stable');
g2 = groupsummary(g1, final_group, ...
    {@(x) mean(x, 'omitnan'), @(x) std(x, 'omitnan'), @(x) sum(isfinite(x))}, 'value');

pooled = table();
for c = 1:numel(final_group)
    pooled.(final_group{c}) = g2.(final_group{c});
end
pooled.mean = g2.fun1_value;
pooled.sd = g2.fun2_value;
pooled.n = g2.fun3_value;
pooled.sem = g2.fun2_value ./ sqrt(max(g2.fun3_value, 1));
end


function T = local_group_mean(T, group_cols)
% Group T by group_cols and replace 'value' with the per-group mean,
% dropping all other (nuisance) columns. NaN-omitting.
g = groupsummary(T, group_cols, 'mean', 'value');
out = table();
for c = 1:numel(group_cols)
    out.(group_cols{c}) = g.(group_cols{c});
end
out.value = g.mean_value;
T = out;
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

% Each curve is one (study_label, condition) series. In single-source
% mode pooled has no study_label column, so we degenerate to condition-only.
has_study = any(strcmp('study_label', pooled.Properties.VariableNames));
if has_study
    series_keys = unique(pooled(:, {'study_label', 'condition'}), 'stable');
else
    series_keys = unique(pooled(:, {'condition'}), 'stable');
end
colors = local_color_palette(height(series_keys));

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

    for k = 1:height(series_keys)
        cond = series_keys.condition(k);
        mask = pooled.region == string(top_regions{r}) & pooled.condition == cond;
        if has_study
            study = series_keys.study_label(k);
            mask = mask & pooled.study_label == study;
            % Series legend: "label | condition" when label is non-empty,
            % else just the condition.
            if strlength(study) > 0
                display_name = sprintf('%s | %s', char(study), char(cond));
            else
                display_name = char(cond);
            end
        else
            display_name = char(cond);
        end
        if ~any(mask), continue; end
        sub = sortrows(pooled(mask, :), 'lag_seconds');
        x = sub.lag_seconds(:);
        y = sub.mean(:);
        e = local_err(pooled, mask, opts.ErrorBand);
        if any(e > 0)
            xx = [x; flipud(x)];
            yy = [y + e; flipud(y - e)];
            fill(ax, xx, yy, colors(k, :), 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
        plot(ax, x, y, '-', 'Color', colors(k, :), 'LineWidth', 1.6, ...
            'DisplayName', display_name);
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


function mask = local_condition_pattern_match(cond_vec, patterns)
% Glob-style match: '*foo' matches anything ending in 'foo', 'foo*'
% matches anything starting with 'foo', '*foo*' matches anything
% containing 'foo'. Patterns without '*' use exact (case-sensitive)
% match. Returns a column logical mask the same height as cond_vec.
cond_str = string(cond_vec);
mask = false(numel(cond_str), 1);
for i = 1:numel(patterns)
    p = char(patterns{i});
    if contains(p, '*')
        % regexptranslate('wildcard', ...) handles MATLAB's glob -> regex
        % conversion (* -> .*, escaping of regex metacharacters).
        rx = ['^', regexptranslate('wildcard', p), '$'];
        hit = ~cellfun('isempty', regexp(cellstr(cond_str), rx, 'once'));
    else
        hit = cond_str == string(p);
    end
    mask = mask | hit(:);
end
end


function labels = local_resolve_condition_labels(patterns, user_labels)
% Per-pattern collapsed label. If the user supplied ConditionLabels
% (parallel to Conditions) use those; otherwise derive a clean label
% from each pattern by stripping '*' and surrounding underscores:
%   '*_heat_start_ttl_1'  -> 'heat_start_ttl_1'
%   'leftface*'           -> 'leftface'
user_labels = cellstr(string(user_labels));
labels = cell(1, numel(patterns));
for i = 1:numel(patterns)
    if numel(user_labels) >= i && ~isempty(user_labels{i})
        labels{i} = user_labels{i};
    else
        lab = strrep(patterns{i}, '*', '');
        lab = regexprep(lab, '^_+', '');   % strip leading underscores
        lab = regexprep(lab, '_+$', '');   % strip trailing underscores
        if isempty(lab), lab = patterns{i}; end
        labels{i} = lab;
    end
end
end


function [keep, new_cond] = local_collapse_conditions(cond_vec, patterns, labels)
% For each row, find the FIRST pattern it matches and relabel the
% condition to that pattern's collapsed label. Rows matching no pattern
% are dropped (keep=false). Returns a full-length keep mask and a
% full-length relabeled condition vector (only meaningful where keep).
cond_str = string(cond_vec);
n = numel(cond_str);
keep = false(n, 1);
new_cond = cond_str;
for i = 1:numel(patterns)
    p = char(patterns{i});
    if contains(p, '*')
        rx = ['^', regexptranslate('wildcard', p), '$'];
        hit = ~cellfun('isempty', regexp(cellstr(cond_str), rx, 'once'));
    else
        hit = cond_str == string(p);
    end
    hit = hit(:);
    % Only assign rows not already claimed by an earlier pattern.
    assign = hit & ~keep;
    new_cond(assign) = string(labels{i});
    keep = keep | hit;
end
end


function tf = local_is_source_struct(x)
% True when x is a struct array carrying labeled input_tables -- accepts
% either a .table field or .input_table field name. The label field is
% required.
tf = isstruct(x) && isfield(x, 'label') && (isfield(x, 'table') || isfield(x, 'input_table'));
end


function sources = local_normalize_sources(input_data)
% Returns a non-empty struct array with fields .label (string) and
% .table (table). A bare table comes back as a single entry with empty
% label, preserving single-source plot behavior.
if istable(input_data)
    sources = struct('label', {""}, 'table', {input_data});
    return
end
n = numel(input_data);
sources = struct('label', cell(1, n), 'table', cell(1, n));
for i = 1:n
    if isfield(input_data(i), 'table')
        sources(i).table = input_data(i).table;
    else
        sources(i).table = input_data(i).input_table;
    end
    sources(i).label = string(input_data(i).label);
    if ~istable(sources(i).table)
        error('plot_hrf_atlas_curves:BadSource', ...
            'sources(%d).table must be a table.', i);
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
