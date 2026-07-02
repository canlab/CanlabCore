function R = hrf_causality(out_dir, varargin)
%HRF_CAUSALITY Systems-level lead/lag from an HRF output dir (IO driver).
%
% End-to-end HRF-informed causality for a study output directory: pulls each
% node's estimated sFIR kernel from the score CSVs, extracts the matching
% node TIMESERIES from the 4-D preprocessed BOLD (signatures via
% pattern-expression, atlas regions via apply_parcellation), deconvolves with
% the per-node kernels, and runs Granger causality with group statistics
% (hrf_causality_analyze). This removes the regional-HRF latency confound that
% otherwise reverses BOLD effective-connectivity directions.
%
% Reads hrf_study_config.mat (written by the SLURM study pipeline) from
% out_dir for the BOLD/events files, run/subject structure, TR, and the
% signature/imageset/atlas specs.
%
% Usage
% -----
%   R = hrf_causality(out_dir)                                  % signatures
%   R = hrf_causality(out_dir, 'Unit','atlas', 'Atlas','ppat')  % regions
%   R = hrf_causality({dir1, dir2})                             % POOL two dirs
%   R = hrf_causality({lf,obs}, 'Condition','rest_stim')        % GC in blocks
%   R = hrf_causality(out_dir, 'MaxRuns',2)                     % quick smoke test
%
% To CONTRAST two task states, run once per condition (or per dir) and pass
% both results to hrf_causality_contrast:
%   Rr = hrf_causality({lf,obs}, 'Condition','rest_stim');
%   Rn = hrf_causality({lf,obs}, 'Condition','nback-stimblock');
%   C  = hrf_causality_contrast(Rr, Rn);   hrf_plot_causality(C);
%
% Inputs
% ------
%   out_dir - a study output directory (with hrf_study_config.mat and
%             *_map_scores.csv), OR a cell of such directories whose subjects
%             are POOLED (e.g. {lf_distractmap, obs_distractmap} bodysites; the
%             same subject id across dirs combines that subject's runs).
%
% Optional (name-value)
% ---------------------
%   'Unit'        - 'signature' (default) or 'atlas'.
%   'Condition'   - restrict Granger to the event blocks of this trial_type
%                   (glob ok; e.g. 'rest_stim', 'nback-stimblock'). Each block
%                   becomes a separate GC realization. Default '' = whole run.
%   'MinSegLen'   - drop condition blocks shorter than this many TRs. Default 20.
%   'Atlas'       - for Unit='atlas': which atlas token (e.g. 'canlab2024',
%                   'ppat'); default = first atlas in the config.
%   'Nodes'       - restrict to these node names (region/signature), glob ok.
%                   Default {} = all nodes of the unit.
%   'KernelModel' - model whose sFIR curve is the deconvolution kernel.
%                   Default 'sfir'. 'KernelObject' default 'beta'.
%   'KernelCondition' - condition (glob) whose HRF is the kernel; '' (default)
%                   averages the HRF across conditions (a cleaner, more
%                   condition-general hemodynamic kernel).
%   'EvokedMode'  - 'both' (default) | 'remove' | 'keep' (see analyze).
%   'Conditional' - conditional MVGC (default false = pairwise).
%   'Order','MaxOrder','Nperm','DeconvMethod' - passed through.
%   'MaxRuns'     - cap number of runs processed (smoke testing). Default Inf.
%   'Verbose'/'doverbose' - progress chatter (default true).
%
% Output
% ------
%   R - the hrf_causality_analyze struct (per-mode net_group/t/p/p_fdr +
%       net_subj), plus .unit, .nodes, .kernel_lags, .run_files.
%
% See also: hrf_causality_analyze, hrf_deconv_timeseries,
%           hrf_granger_causality, hrf_apply_maps_to_wholebrain.

p = inputParser;
p.addRequired('out_dir', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('Unit', 'signature', @(x) ischar(x) || isstring(x));
p.addParameter('Atlas', '', @(x) ischar(x) || isstring(x));
p.addParameter('Nodes', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('KernelModel', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('KernelObject', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('KernelCondition', '', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('EvokedMode', 'both', @(x) ischar(x) || isstring(x));
p.addParameter('Conditional', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Order', 'bic', @(x) (ischar(x) || isstring(x)) || isscalar(x));
p.addParameter('MaxOrder', 10, @(x) isscalar(x) && x >= 1);
p.addParameter('Nperm', 0, @(x) isscalar(x) && x >= 0);
p.addParameter('Correction', 'fdr', @(x) ischar(x) || isstring(x));
p.addParameter('GroupNperm', 5000, @(x) isscalar(x) && x >= 100);
p.addParameter('DeconvMethod', 'ridge', @(x) ischar(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('MinSegLen', 20, @(x) isscalar(x) && x >= 4);
p.addParameter('MaxRuns', Inf, @(x) isscalar(x) && x >= 1);
p.addParameter('ReturnData', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(out_dir, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);
if ~isempty(opts.doverbose), verbose = logical(opts.doverbose); end
od_list = local_dir_list(out_dir);
unit = lower(char(opts.Unit));

cfgs = cell(1, numel(od_list));
for di = 1:numel(od_list), cfgs{di} = local_load_config(od_list{di}); end
TR = local_pa(cfgs{1}.pipeline_args, 'TR', NaN);
if ~isfinite(TR), error('hrf_causality:NoTR', 'Could not read TR from config.pipeline_args.'); end

% ---- 1. kernels from the score CSVs (group-mean sFIR curve per node), pooled across dirs
[kernels, nodes_all, lags] = local_build_kernels(od_list, unit, opts);
[kernels, nodes] = local_select_kernels(kernels, nodes_all, opts.Nodes);
if isempty(nodes), error('hrf_causality:NoNodes', 'No %s nodes after Nodes filter.', unit); end
cond_txt = local_cond_text(opts.Condition);
if verbose
    fprintf('hrf_causality: %d %s node(s) over %d dir(s); kernel=%s/%s, %d lags, TR=%.3g; condition=%s\n', ...
        numel(nodes), unit, numel(od_list), char(opts.KernelModel), char(opts.KernelObject), numel(lags), TR, cond_txt);
end

% Fast path: build ONE image_vector holding only the requested node maps, so
% each run applies just those instead of the whole signature+imageset set (44+
% maps). Huge win when few nodes are needed (e.g. mediation M='NPS' -> 1 map).
% Built once; falls back to the full apply if it can't be assembled.
sig_subset = [];
if strcmp(unit, 'signature')
    try
        sig_subset = local_node_map_subset(cfgs{1}.signature_sets, cfgs{1}.image_sets, nodes);
    catch
        sig_subset = [];
    end
    if verbose && ~isempty(sig_subset)
        fprintf('  fast extract: applying %d requested map(s)/run instead of the full set\n', size(sig_subset.dat, 2));
    end
end

% ---- 2. per-run node timeseries from the 4-D BOLD (across all dirs) ------
tsRuns = {}; confRuns = {}; segRuns = {}; subjUsed = {}; usedFiles = {}; usedEvents = {};
total = 0;
for di = 1:numel(od_list)
    cfg = cfgs{di};
    fmri_files = cellstr(string(cfg.fmri_files));
    events_files = local_cellstr_or_empty(cfg, 'events_files', numel(fmri_files));
    subjects = cellstr(string(cfg.subject_ids));
    atlas_obj = [];
    if strcmp(unit, 'atlas'), atlas_obj = local_pick_atlas(cfg, opts.Atlas); end
    for r = 1:numel(fmri_files)
        if total >= opts.MaxRuns, break; end
        f = fmri_files{r};
        if exist(f, 'file') ~= 2
            warning('hrf_causality:MissingBOLD', 'BOLD not found, skipping: %s', f); continue
        end
        if verbose, fprintf('  [dir %d run %d] %s\n', di, r, local_short(f)); end
        bold = fmri_data(f, 'noverbose');
        ts = local_extract_timeseries(bold, unit, nodes, cfg, atlas_obj, sig_subset);
        if isempty(ts), warning('hrf_causality:NoTS', 'No timeseries for %s; skipping.', local_short(f)); continue; end
        T = size(ts, 1);
        ef = ''; if r <= numel(events_files), ef = events_files{r}; end
        tsRuns{end + 1} = ts; %#ok<AGROW>
        confRuns{end + 1} = local_events_design(events_files, r, T, TR); %#ok<AGROW>
        segRuns{end + 1} = local_condition_mask(events_files, r, T, TR, opts.Condition); %#ok<AGROW>
        subjUsed{end + 1} = subjects{min(r, numel(subjects))}; %#ok<AGROW>
        usedFiles{end + 1} = f; %#ok<AGROW>
        usedEvents{end + 1} = ef; %#ok<AGROW>
        total = total + 1;
    end
end
if isempty(tsRuns), error('hrf_causality:NoRuns', 'No runs produced timeseries.'); end

% Prune nodes missing/constant in ANY run so kernels + ts stay column-aligned.
valid = true(1, numel(nodes));
for r = 1:numel(tsRuns)
    valid = valid & (all(isfinite(tsRuns{r}), 1) & std(tsRuns{r}, 0, 1) > 0);
end
if ~all(valid) && verbose, fprintf('  pruning %d node(s) missing/constant in some run\n', sum(~valid)); end
kernels = kernels(:, valid); nodes = nodes(valid);
for r = 1:numel(tsRuns), tsRuns{r} = tsRuns{r}(:, valid); end
if isempty(nodes), error('hrf_causality:NoValidNodes', 'No nodes valid across all runs.'); end

% Early exit: hand back the extracted per-run timeseries + kernels + events so
% other methods (e.g. hrf_causality_mediation) can reuse the expensive BOLD
% extraction without re-loading the 4-D data.
if logical(opts.ReturnData)
    R = struct('tsRuns', {tsRuns}, 'kernels', kernels, 'nodes', {nodes}, ...
        'subjects', {subjUsed}, 'events_files', {usedEvents}, 'run_files', {usedFiles(:)'}, ...
        'TR', TR, 'dirs', {od_list(:)'}, 'unit', unit, 'kernel_lags', lags);
    return
end

% Warn if a condition was requested but matched no event blocks anywhere.
if ~strcmp(cond_txt, '(whole run)') && ~any(cellfun(@(m) ~isempty(m) && any(m), segRuns))
    warning('hrf_causality:NoConditionBlocks', ...
        'Condition ''%s'' matched no event blocks in any run; check trial_type labels.', cond_txt);
end

R = hrf_causality_analyze(tsRuns, kernels, ...
    'Subjects', subjUsed, 'Nodes', nodes, 'EvokedMode', opts.EvokedMode, ...
    'Confounds', confRuns, 'Segments', segRuns, 'MinSegLen', opts.MinSegLen, ...
    'Conditional', opts.Conditional, 'Order', opts.Order, 'MaxOrder', opts.MaxOrder, ...
    'Nperm', opts.Nperm, 'DeconvMethod', opts.DeconvMethod, ...
    'Correction', opts.Correction, 'GroupNperm', opts.GroupNperm, 'doverbose', verbose);
R.unit = unit; R.kernel_lags = lags; R.run_files = usedFiles(:);
R.dirs = od_list(:)'; R.condition = cond_txt;
end


% =========================================================================
function cfg = local_load_config(od)
cf = fullfile(od, 'hrf_study_config.mat');
if exist(cf, 'file') ~= 2
    error('hrf_causality:NoConfig', 'hrf_study_config.mat not found in %s', od);
end
S = load(cf, 'config'); cfg = S.config;
end


function L = local_dir_list(x)
% Normalize the out_dir argument to a cellstr list of directories to pool.
if iscell(x)
    L = cellfun(@(c) char(string(c)), x, 'uni', 0);
elseif isstring(x) && ~isscalar(x)
    L = cellstr(x);
else
    L = {char(string(x))};
end
end


function s = local_cond_text(c)
cc = cellstr(string(c));
cc = cc(~cellfun(@(z) isempty(strtrim(z)), cc));
if isempty(cc), s = '(whole run)'; else, s = strjoin(cc, '|'); end
end


function mask = local_condition_mask(events_files, r, T, TR, condition)
% Logical [T x 1] marking the run frames inside the requested condition's
% event blocks (trial_type glob-matched). [] if no condition / no events ->
% the analyze step then treats the whole run as one segment.
mask = [];
cond = cellstr(string(condition));
cond = cond(~cellfun(@(z) isempty(strtrim(z)), cond));
if isempty(cond), return; end
if isempty(events_files) || r > numel(events_files), return; end
ef = events_files{r};
if isempty(ef) || exist(ef, 'file') ~= 2, return; end
try
    E = readtable(ef, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
catch
    return
end
v = E.Properties.VariableNames;
if ~all(ismember({'onset', 'trial_type'}, v)), return; end
onset = double(E.onset);
dur = zeros(size(onset)); if any(strcmp('duration', v)), dur = double(E.duration); end
tt = cellstr(string(E.trial_type));
sel = local_node_mask(tt, cond);          % glob trial-type match
mask = false(T, 1);
rows = find(sel);
for j = rows(:)'
    on = max(1, floor(onset(j) / TR) + 1);
    off = min(T, on + max(1, round(dur(j) / TR)) - 1);
    if on <= T, mask(on:off) = true; end
end
end


function [K, nodes, lags] = local_build_kernels(od, unit, opts)
% Group-mean sFIR curve per node, averaged across the kernel condition(s),
% pooled across all subjects' score CSVs (and all dirs) for the model/object.
prefix = local_unit_prefix(unit);
model = lower(char(opts.KernelModel)); object = lower(char(opts.KernelObject));
dirs = od; if ~iscell(dirs), dirs = {char(dirs)}; end
csvs = [];
for di = 1:numel(dirs)
    cc = dir(fullfile(dirs{di}, sprintf('*_hrf_%s_%s_map_scores.csv', model, object)));
    if isempty(cc), cc = dir(fullfile(dirs{di}, sprintf('*_hrf_%s_map_scores.csv', object))); end
    csvs = [csvs; cc]; %#ok<AGROW>
end
if isempty(csvs)
    error('hrf_causality:NoKernelCsv', 'No %s/%s score CSVs in the given dir(s).', model, object);
end
acc = containers.Map('KeyType', 'char', 'ValueType', 'any');   % node -> running [L x k] curves
lags = [];
for c = 1:numel(csvs)
    T = readtable(fullfile(csvs(c).folder, csvs(c).name), 'TextType', 'string');
    v = T.Properties.VariableNames;
    if ~any(strcmp('lag_index', v)) && ~any(strcmp('lag_seconds', v)), continue; end
    [lagkey, this_lags] = local_lagkey(T);
    cond = string(T.condition);
    cmask = local_cond_mask(cond, opts.KernelCondition);
    cols = v(startsWith(v, prefix) & ~endsWith(v, '_se'));
    for i = 1:numel(cols)
        nm = local_node_from_col(cols{i}, unit);
        curve = local_curve_for_node(T, cols{i}, cmask, lagkey, this_lags);
        if isempty(curve), continue; end
        if isempty(lags) || numel(this_lags) > numel(lags), lags = this_lags(:); end
        if isKey(acc, nm), acc(nm) = [acc(nm), curve(:)]; else, acc(nm) = curve(:); end
    end
end
nodes = acc.keys; nodes = nodes(:)';
L = numel(lags);
K = zeros(L, numel(nodes));
for n = 1:numel(nodes)
    cv = acc(nodes{n});
    cv = local_pad_rows(cv, L);
    K(:, n) = mean(cv, 2, 'omitnan');
end
end


function curve = local_curve_for_node(T, col, cmask, lagkey, lags)
% Average the HRF curve over the selected conditions: group by lag, mean.
y = double(T.(col));
lg = double(T.(lagkey));
y = y(cmask); lg = lg(cmask);
[ul, ~, gi] = unique(lg, 'stable');
m = accumarray(gi, y, [], @(x) mean(x, 'omitnan'));
% reorder to ascending lag matching `lags`
[~, ord] = sort(ul);
m = m(ord); uls = ul(ord);
curve = nan(numel(lags), 1);
[tf, loc] = ismember(round(lags(:), 6), round(uls(:), 6));
curve(tf) = m(loc(tf));
end


function ts = local_extract_timeseries(bold, unit, nodes, cfg, atlas_obj, sig_subset)
% Return [T x numel(nodes)] timeseries in `nodes` order (NaN column if a node
% is unavailable in this run). Columns stay aligned with the kernels; the
% driver prunes any node that is invalid across runs after the loop.
T = size(bold.dat, 2);
ts = nan(T, numel(nodes));
if strcmp(unit, 'signature')
    metric = local_field_or(cfg, 'similarity_metric', 'dotproduct');
    if ~isempty(sig_subset)
        % fast path: apply ONLY the requested node maps (named by node)
        sc = hrf_apply_maps_to_wholebrain(bold, 'SignatureSets', {}, ...
            'ImageSets', sig_subset, 'SimilarityMetric', metric);
    else
        sc = hrf_apply_maps_to_wholebrain(bold, ...
            'SignatureSets', cfg.signature_sets, 'ImageSets', cfg.image_sets, ...
            'SimilarityMetric', metric);
    end
    v = sc.Properties.VariableNames;
    for n = 1:numel(nodes)
        col = local_match_score_col(v, nodes{n});
        if ~isempty(col), ts(:, n) = double(sc.(col)); end
    end
else
    parc = apply_parcellation(bold, atlas_obj);
    if isstruct(parc) && isfield(parc, 'dat'), P = parc.dat; else, P = parc; end
    if size(P, 1) ~= T && size(P, 2) == T, P = P'; end   % want [T x regions]
    labels = cellstr(string(atlas_obj.labels));
    for n = 1:numel(nodes)
        li = find(strcmp(labels, nodes{n}), 1);
        if ~isempty(li) && li <= size(P, 2), ts(:, n) = P(:, li); end
    end
end
end


function obj = local_node_map_subset(sig_sets, image_sets, want)
% Build one image_vector holding just the requested node maps -- searched
% across the signature sets AND imagesets and named by node -- so the extractor
% applies ONLY those maps per run instead of the whole set. Returns [] if none
% of the wanted nodes is a map (caller falls back to the full apply).
obj = []; got = {};
want = cellstr(string(want));
setspecs = [cellstr(string(sig_sets)); cellstr(string(image_sets))];
setspecs = setspecs(~cellfun(@(s) isempty(strtrim(s)), setspecs));
for s = 1:numel(setspecs)
    if all(ismember(lower(want), lower(got))), break; end     % already have all
    try
        [o, nm] = load_image_set(setspecs{s}, 'noverbose');
    catch
        continue
    end
    nm = cellstr(string(nm));
    for w = 1:numel(want)
        if any(strcmpi(got, want{w})), continue; end
        hit = find(strcmpi(nm, want{w}), 1);
        if isempty(hit), continue; end
        sub = get_wh_image(o, hit);
        if isempty(obj)
            obj = sub;
        else
            try, sub = resample_space(sub, obj); catch, end
            obj.dat = [obj.dat, sub.dat];
        end
        got{end + 1} = want{w}; %#ok<AGROW>
    end
end
if isempty(obj), return; end
try, obj.metadata_table = table(); catch, end
obj.image_names = char(cellstr(string(got)));
end


function conf = local_events_design(events_files, r, T, TR)
% Build an unconvolved [T x q] task design (one boxcar per trial_type) for the
% 'remove' evoked mode. Best-effort BIDS events.tsv parsing; [] if unavailable.
conf = [];
if isempty(events_files) || r > numel(events_files), return; end
ef = events_files{r};
if isempty(ef) || exist(ef, 'file') ~= 2, return; end
try
    E = readtable(ef, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
catch
    return
end
v = E.Properties.VariableNames;
if ~all(ismember({'onset', 'trial_type'}, v)), return; end
onset = double(E.onset);
dur = zeros(size(onset)); if any(strcmp('duration', v)), dur = double(E.duration); end
tt = string(E.trial_type);
utt = unique(tt, 'stable');
conf = zeros(T, numel(utt));
for k = 1:numel(utt)
    rows = find(tt == utt(k));
    for j = rows(:)'
        on = max(1, floor(onset(j) / TR) + 1);
        off = min(T, on + max(1, round(dur(j) / TR)) - 1);
        if on <= T, conf(on:off, k) = 1; end
    end
end
conf = conf(:, any(conf ~= 0, 1));
end


% =========================================================================
function pfx = local_unit_prefix(unit)
switch unit
    case 'signature', pfx = 'sig_';
    case 'atlas',     pfx = 'atlas_';
    case 'imageset',  pfx = 'map_';
    otherwise, error('hrf_causality:Unit', 'Unit must be signature|atlas|imageset.');
end
end

function nm = local_node_from_col(col, unit)
parts = strsplit(col, '_');
if numel(parts) < 3, nm = col; return; end
switch unit
    case 'atlas'
        if ismember(parts{end}, {'mean', 'meanL1', 'sum'}), nm = strjoin(parts(3:end-1), '_');
        else, nm = strjoin(parts(3:end), '_'); end
    otherwise
        nm = strjoin(parts(3:end), '_');
end
end

function col = local_match_score_col(v, node)
% Find the sig_/map_ column whose parsed name equals `node`.
col = '';
cand = v(startsWith(v, 'sig_') | startsWith(v, 'map_'));
cand = cand(~endsWith(cand, '_se'));
for i = 1:numel(cand)
    parts = strsplit(cand{i}, '_');
    if numel(parts) >= 3 && strcmp(strjoin(parts(3:end), '_'), node), col = cand{i}; return; end
end
end

function [lagkey, lags] = local_lagkey(T)
v = T.Properties.VariableNames;
if any(strcmp('lag_seconds', v)), lagkey = 'lag_seconds'; else, lagkey = 'lag_index'; end
lags = unique(double(T.(lagkey)), 'stable');
lags = sort(lags);
end

function mask = local_cond_mask(cond, want)
want = cellstr(string(want));
want = want(~cellfun(@(s) isempty(strtrim(s)), want));
if isempty(want), mask = true(numel(cond), 1); return; end
mask = false(numel(cond), 1);
for i = 1:numel(want)
    pat = strtrim(want{i});
    if any(pat == '*' | pat == '?')
        rx = ['^', regexptranslate('wildcard', pat), '$'];
        mask = mask | ~cellfun('isempty', regexp(cellstr(cond), rx, 'once'));
    else
        mask = mask | (cond == string(pat));
    end
end
end

function M = local_pad_rows(M, L)
if size(M, 1) < L, M = [M; nan(L - size(M, 1), size(M, 2))]; elseif size(M, 1) > L, M = M(1:L, :); end
end

function [K, nodes] = local_select_kernels(K, nodes_all, want)
% Keep kernel columns whose node name matches `want` (glob ok); empty => all.
keep = local_node_mask(nodes_all, want);
K = K(:, keep); nodes = nodes_all(keep);
end

function keep = local_node_mask(nodes, want)
want = cellstr(string(want));
want = want(~cellfun(@(s) isempty(strtrim(s)), want));
if isempty(want), keep = true(1, numel(nodes)); return; end
keep = false(1, numel(nodes));
for i = 1:numel(nodes)
    for j = 1:numel(want)
        pat = strtrim(want{j});
        if any(pat == '*' | pat == '?')
            hit = ~isempty(regexp(nodes{i}, ['^', regexptranslate('wildcard', pat), '$'], 'once'));
        else
            hit = strcmpi(nodes{i}, pat);
        end
        if hit, keep(i) = true; end
    end
end
end

function atl = local_pick_atlas(cfg, want)
objs = local_field_or(cfg, 'atlas_objs', {});
names = cellstr(string(local_field_or(cfg, 'atlas_names', {})));
if isempty(objs) && ~isempty(local_field_or(cfg, 'atlas_obj', []))
    atl = cfg.atlas_obj; return
end
if isempty(objs), error('hrf_causality:NoAtlas', 'No atlas in config for Unit=atlas.'); end
idx = 1;
if ~isempty(want)
    hit = find(strcmpi(names, char(want)), 1);
    if ~isempty(hit), idx = hit; else, warning('hrf_causality:AtlasName', 'Atlas ''%s'' not in config; using %s.', char(want), names{1}); end
end
atl = objs{idx};
end

function val = local_pa(pa, key, default)
val = default;
for i = 1:2:numel(pa) - 1
    if ischar(pa{i}) && strcmpi(pa{i}, key), val = pa{i + 1}; return; end
end
end

function val = local_field_or(s, f, default)
if isfield(s, f) && ~isempty(s.(f)), val = s.(f); else, val = default; end
end

function c = local_cellstr_or_empty(cfg, f, n)
if isfield(cfg, f) && ~isempty(cfg.(f)), c = cellstr(string(cfg.(f))); else, c = repmat({''}, 1, n); end
end

function s = local_short(f)
[~, n, e] = fileparts(char(f)); s = [n e];
end
