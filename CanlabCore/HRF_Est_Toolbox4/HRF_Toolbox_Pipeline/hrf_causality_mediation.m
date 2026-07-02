function R = hrf_causality_mediation(source, varargin)
%HRF_CAUSALITY_MEDIATION Trial-level mediation X -> M -> Y on HRF-deconvolved nodes.
%
% Driver for the mediation method: extracts per-trial response amplitudes for
% each node from the HRF-DECONVOLVED BOLD (the same deconvolution used by the
% Granger pipeline), builds the trial-level X, M, Y variables from those node
% amplitudes and/or the BIDS events, and runs multilevel mediation
% (hrf_mediation_analyze). Answers e.g. "does signature/region M mediate the
% effect of stimulus on rating", "does A mediate the effect of stimulus on B",
% or "does A mediate the rest-vs-nback condition effect".
%
% Any of X / M / Y can be:
%   * a NODE name (signature/region)  -> that node's per-trial amplitude
%   * an EVENTS column ('temp','rating',...) -> that column's per-trial value
%   * 'condition:<glob>'              -> 0/1 indicator of trial_type match
% so all outcome types are available through one API.
%
% :Usage:
% ::
%     % stimulus intensity -> NPS -> pain rating, pooled over bodysites
%     R = hrf_causality_mediation({lf,obs}, 'X','temp', 'M','NPS', 'Y','rating', ...
%             'TrialType','*stimblock*');
%     % does NPS mediate the effect of stim on SIIPS (node -> node)?
%     R = hrf_causality_mediation({lf,obs}, 'X','temp','M','NPS','Y','SIIPS');
%     % reuse an extraction prep instead of re-loading BOLD:
%     prep = hrf_causality({lf,obs}, 'Unit','signature', 'ReturnData',true);
%     R = hrf_causality_mediation(prep, 'X','temp','M','NPS','Y','rating');
%
% :Inputs:
%   **source:** a dir / cell of dirs (passed to hrf_causality for extraction),
%             OR a prep struct from hrf_causality(...,'ReturnData',true).
%
% :Required name-value:
%   **'X','M','Y':** the three mediation variables (specs as above).
%
% :Optional Inputs:
%   **'TrialType':** glob on trial_type selecting which events are the trials
%             (the M / amplitude anchor). Use the stimulus BLOCK, e.g.
%             'nback-stimblock' or 'rest_stim' (NOT '*stimblock*', which also
%             grabs the zero-duration *_ttl_* markers). Default '' => rows with
%             a finite 'rating'.
%   **'TrialWindow':** seconds post-onset to average the proxy over for a
%             trial's amplitude. Default [] => the event's own duration.
%   **'PairWindow':** seconds. When an events column used for X/Y is NaN on the
%             trial row (e.g. 'rating', which lives on a separate rating event),
%             its value is paired from the nearest event that has it -- the next
%             one within PairWindow, else the nearest overall. Default 90.
%   **'Amplitude':** 'mean' (default) or 'peak' of the proxy over the window.
%   **'DeconvMethod','Lambda':** passed to hrf_deconv_timeseries.
%   **'Standardize','Nboot':** passed to hrf_mediation_analyze.
%   Extraction passthrough (when source is dirs): 'Unit','Nodes','Atlas',
%   'KernelModel','KernelObject','KernelCondition','MaxRuns'.
%   **'Verbose'/'doverbose':** default true.
%
% :Output:
%   **R:** the hrf_mediation_analyze struct (a/b/cp/c/ab paths + stats), plus
%          .x/.m/.y specs, .trialtype, .dirs.
%
% See also: hrf_mediation_analyze, hrf_causality, hrf_deconv_timeseries.

p = inputParser;
p.addRequired('source');
p.addParameter('X', '', @(x) ischar(x) || isstring(x));
p.addParameter('M', '', @(x) ischar(x) || isstring(x));
p.addParameter('Y', '', @(x) ischar(x) || isstring(x));
p.addParameter('TrialType', '', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('TrialWindow', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('PairWindow', 90, @(x) isscalar(x) && x > 0);
p.addParameter('Amplitude', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('DeconvMethod', 'ridge', @(x) ischar(x) || isstring(x));
p.addParameter('Lambda', [], @(x) isempty(x) || isscalar(x));
p.addParameter('Standardize', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Nboot', 5000, @(x) isscalar(x) && x >= 100);
% extraction passthrough
p.addParameter('Unit', 'signature', @(x) ischar(x) || isstring(x));
p.addParameter('Nodes', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Atlas', '', @(x) ischar(x) || isstring(x));
p.addParameter('KernelModel', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('KernelObject', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('KernelCondition', '', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('MaxRuns', Inf, @(x) isscalar(x) && x >= 1);
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(source, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);
if ~isempty(opts.doverbose), verbose = logical(opts.doverbose); end
if isempty(char(opts.X)) || isempty(char(opts.M)) || isempty(char(opts.Y))
    error('hrf_causality_mediation:Specs', 'X, M and Y must all be specified.');
end

% ---- 1. obtain the extraction prep (proxy-ready timeseries) --------------
if isstruct(source) && isfield(source, 'tsRuns')
    prep = source;
else
    % Restrict extraction to only the node-valued X/M/Y specs (events columns
    % like temp/rating and 'condition:' specs are resolved later, not extracted)
    % so the BOLD apply is done for just those maps -- a large speedup.
    ext_nodes = {};
    for spc = {char(opts.X), char(opts.M), char(opts.Y)}
        if ~startsWith(spc{1}, 'condition:'), ext_nodes{end + 1} = spc{1}; end %#ok<AGROW>
    end
    if ~isempty(opts.Nodes), ext_nodes = [cellstr(string(opts.Nodes)); ext_nodes(:)]; end
    ext_nodes = unique(ext_nodes, 'stable');
    prep = hrf_causality(source, 'ReturnData', true, 'Unit', opts.Unit, 'Nodes', ext_nodes, ...
        'Atlas', opts.Atlas, 'KernelModel', opts.KernelModel, 'KernelObject', opts.KernelObject, ...
        'KernelCondition', opts.KernelCondition, 'MaxRuns', opts.MaxRuns, 'doverbose', verbose);
end
nodes = cellstr(string(prep.nodes));
TR = prep.TR;
nrun = numel(prep.tsRuns);

% ---- 2. per-run: deconvolve, then per-trial node amplitudes -------------
subj_amp = containers.Map('KeyType', 'char', 'ValueType', 'any');  % subject -> bundle
for r = 1:nrun
    proxy = hrf_deconv_timeseries(prep.tsRuns{r}, prep.kernels, ...
        'Method', opts.DeconvMethod, 'Lambda', opts.Lambda);
    b = local_trial_bundle(prep.events_files{r}, proxy, nodes, TR, opts);
    if isempty(b.amp), continue; end
    sid = char(string(prep.subjects{r}));
    if isKey(subj_amp, sid), subj_amp(sid) = local_cat_bundle(subj_amp(sid), b);
    else, subj_amp(sid) = b; end
end
sids = subj_amp.keys;
if isempty(sids), error('hrf_causality_mediation:NoTrials', 'No trials extracted (check TrialType / events).'); end

% ---- 3. resolve X/M/Y to per-subject trial vectors ----------------------
Xc = {}; Mc = {}; Yc = {};
for s = 1:numel(sids)
    b = subj_amp(sids{s});
    x = local_resolve(opts.X, b, nodes);
    m = local_resolve(opts.M, b, nodes);
    y = local_resolve(opts.Y, b, nodes);
    Xc{end + 1} = x; Mc{end + 1} = m; Yc{end + 1} = y; %#ok<AGROW>
end

% ---- 4. mediation -------------------------------------------------------
R = hrf_mediation_analyze(Xc, Mc, Yc, 'Names', {char(opts.X), char(opts.M), char(opts.Y)}, ...
    'Standardize', opts.Standardize, 'Nboot', opts.Nboot, 'doverbose', verbose);
R.x = char(opts.X); R.m = char(opts.M); R.y = char(opts.Y);
R.trialtype = char(string(opts.TrialType));
if isfield(prep, 'dirs'), R.dirs = prep.dirs; end
end


% =========================================================================
function b = local_trial_bundle(ef, proxy, nodes, TR, opts)
% Per-trial node amplitudes + passthrough events columns for one run.
b = struct('amp', [], 'trial_type', strings(0, 1), 'cols', struct());
if isempty(ef) || exist(ef, 'file') ~= 2, return; end
try
    E = readtable(ef, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
catch
    return
end
v = E.Properties.VariableNames;
if ~any(strcmp('onset', v)), return; end
onset = double(E.onset);
dur = zeros(size(onset)); if any(strcmp('duration', v)), dur = double(E.duration); end
tt = strings(numel(onset), 1); if any(strcmp('trial_type', v)), tt = string(E.trial_type); end

sel = local_trial_select(E, tt, opts.TrialType);
rows = find(sel);
if isempty(rows), return; end

T = size(proxy, 1); N = numel(nodes);
amp = nan(numel(rows), N);
use_peak = strcmpi(char(opts.Amplitude), 'peak');
for i = 1:numel(rows)
    j = rows(i);
    on = max(1, floor(onset(j) / TR) + 1);
    if isempty(opts.TrialWindow), w = max(dur(j), TR); else, w = opts.TrialWindow; end
    off = min(T, on + max(1, round(w / TR)) - 1);
    seg = proxy(on:off, :);
    if isempty(seg), continue; end
    if use_peak
        [~, pk] = max(abs(seg), [], 1);
        amp(i, :) = seg(sub2ind(size(seg), pk, 1:N));
    else
        amp(i, :) = mean(seg, 1, 'omitnan');
    end
end
b.amp = amp;
b.trial_type = tt(rows);
b.cols = struct();
% Per-trial value of every numeric events column, PAIRED: use the trial row's
% own value if finite; otherwise take the nearest event (preferring the next
% one within PairWindow s) that has a finite value. This is what links a heat
% trial to its rating, which is recorded on a separate 'rating' event.
for k = 1:numel(v)
    nm = v{k};
    if ismember(nm, {'onset', 'duration', 'trial_type'}), continue; end
    col = E.(nm);
    if ~isnumeric(col), continue; end
    b.cols.(matlab.lang.makeValidName(nm)) = local_paired_col(double(col), onset, rows, opts.PairWindow);
end
end


function vals = local_paired_col(colall, onset_all, rows, pairwin)
% For each anchor trial: its own value if finite, else the value of the
% nearest event with a finite value -- preferring the next event within
% pairwin seconds (e.g. the rating that follows the stimulus), else the
% nearest one overall.
fin = find(isfinite(colall));
vals = nan(numel(rows), 1);
for i = 1:numel(rows)
    j = rows(i);
    if isfinite(colall(j)), vals(i) = colall(j); continue; end
    if isempty(fin), continue; end
    d = onset_all(fin) - onset_all(j);
    foll = fin(d >= -1e-6 & d <= pairwin);
    if ~isempty(foll)
        [~, mi] = min(onset_all(foll) - onset_all(j));
        vals(i) = colall(foll(mi));
    else
        [~, mi] = min(abs(onset_all(fin) - onset_all(j)));
        vals(i) = colall(fin(mi));
    end
end
end


function a = local_cat_bundle(a, b)
a.amp = [a.amp; b.amp];
a.trial_type = [a.trial_type; b.trial_type];
fn = union(fieldnames(a.cols), fieldnames(b.cols));
for k = 1:numel(fn)
    va = local_field_or_nan(a.cols, fn{k}, size(a.amp, 1) - size(b.amp, 1));
    vb = local_field_or_nan(b.cols, fn{k}, size(b.amp, 1));
    a.cols.(fn{k}) = [va; vb];
end
end

function v = local_field_or_nan(s, f, n)
if isfield(s, f), v = s.(f)(:); else, v = nan(n, 1); end
end


function v = local_resolve(spec, b, nodes)
spec = char(string(spec));
ni = find(strcmpi(nodes, spec), 1);
if ~isempty(ni), v = b.amp(:, ni); return; end
if startsWith(spec, 'condition:')
    v = double(local_glob(b.trial_type, strtrim(spec(11:end)))); return
end
key = matlab.lang.makeValidName(spec);
if isfield(b.cols, key), v = b.cols.(key)(:); return; end
if isfield(b.cols, spec), v = b.cols.(spec)(:); return; end
if any(local_glob(b.trial_type, spec))
    v = double(local_glob(b.trial_type, spec)); return
end
error('hrf_causality_mediation:Unresolved', ...
    'Could not resolve ''%s'' as a node, events column, or trial_type.', spec);
end


function sel = local_trial_select(E, tt, trialtype)
pats = cellstr(string(trialtype));
pats = pats(~cellfun(@(s) isempty(strtrim(s)), pats));
if isempty(pats)
    % default: rows that have a finite rating (the behavioural trials)
    if any(strcmp('rating', E.Properties.VariableNames))
        sel = isfinite(double(E.rating));
    else
        sel = true(numel(tt), 1);
    end
    return
end
sel = false(numel(tt), 1);
for i = 1:numel(pats), sel = sel | local_glob(tt, pats{i}); end
end


function tf = local_glob(tt, pat)
tt = cellstr(string(tt));
pat = strtrim(pat);
if any(pat == '*' | pat == '?')
    rx = ['^', regexptranslate('wildcard', pat), '$'];
    tf = ~cellfun('isempty', regexp(tt, rx, 'once'));
else
    tf = strcmpi(tt, pat);
end
tf = tf(:);
end
