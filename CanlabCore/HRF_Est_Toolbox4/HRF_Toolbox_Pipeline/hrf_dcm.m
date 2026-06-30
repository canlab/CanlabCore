function R = hrf_dcm(source, varargin)
%HRF_DCM Bilinear DCM directed connectivity via SPM (EXPERIMENTAL).
%
% Estimates a deterministic bilinear Dynamic Causal Model (SPM) on a SMALL,
% user-specified set of nodes and returns the posterior endogenous
% connectivity as a directed net-flow result compatible with
% hrf_plot_causality / hrf_causality_contrast. DCM is generative: it fits the
% hemodynamic response itself, so it is given the RAW node timeseries (NOT the
% deconvolved proxies that Granger/mediation use), with the task conditions as
% driving inputs.
%
% Status: the DCM struct construction and direction recovery are VALIDATED on
% synthetic ground truth -- spm_dcm_estimate converges on the built struct and
% recovers an injected 1->2 influence (net>0, correct sign). It has not yet
% been run at full study scale (each estimate is slow and needs the 4-D BOLD,
% and node selection matters), so treat first real-data results as provisional
% and prefer PEB on the returned .gcm for group inference (see below).
%
% :Usage:
% ::
%     R = hrf_dcm({lf,obs}, 'Unit','atlas', 'Atlas','ppat', ...
%             'Nodes',{'Thal_VPLM_L','dpIns_L','dACC'}, ...
%             'Conditions',{'rest_stim','nback-stimblock'});
%     hrf_plot_causality(R, 'Mode','dcm');
%     % group Bayesian inference (recommended):  PEB = spm_dcm_peb(R.gcm(:));
%
% :Inputs:
%   **source:** dirs / cell of dirs (passed to hrf_causality for extraction),
%             OR a prep struct from hrf_causality(...,'ReturnData',true).
%
% :Required name-value:
%   **'Nodes':**       the DCM node list (<= MaxNodes; DCM does not scale).
%   **'Conditions':**  trial_type(s) (glob) used as the model inputs (driving,
%                      and modulatory if 'B' is set). Built from the BIDS events.
%
% :Optional Inputs:
%   **'A':** endogenous structure -- 'full' (default, all-to-all + self) or an
%            [n x n] 0/1 matrix.
%   **'C':** driving inputs -- 'all' (default, every input drives every node)
%            or an [n x m] 0/1 matrix.
%   **'B':** modulatory -- 'none' (default), 'all', a condition name (that
%            input modulates all connections), or an [n x n x m] 0/1 array.
%   **'TE':** echo time, seconds. Default 0.04.
%   **'Bins':** microtime bins per TR for the inputs. Default 16.
%   **'MaxNodes':** guard; error if more nodes requested. Default 8.
%   **'MaxRuns':** cap runs (smoke test). Default Inf.
%   Extraction passthrough: 'Unit','Atlas','KernelModel','KernelObject',
%   'KernelCondition'. (Kernels are built by the extractor but unused by DCM.)
%   **'Verbose'/'doverbose':** default true.
%
% :Output:
%   **R:** struct with .modes={'dcm'} and .dcm.net_group/.t/.p/.p_fdr/.net_subj
%          (directed i->j net flow, group one-sample t across subjects), .nodes,
%          .subjects, .A_group (posterior A, SPM target<-source convention),
%          and .gcm (cell of per-run estimated DCM structs, for spm_dcm_peb).
%
% See also: hrf_causality, hrf_plot_causality, hrf_causality_contrast,
%           spm_dcm_estimate, spm_dcm_peb.

p = inputParser;
p.addRequired('source');
p.addParameter('Nodes', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Conditions', {}, @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('A', 'full', @(x) (ischar(x) || isstring(x)) || isnumeric(x));
p.addParameter('C', 'all', @(x) (ischar(x) || isstring(x)) || isnumeric(x));
p.addParameter('B', 'none', @(x) (ischar(x) || isstring(x)) || isnumeric(x));
p.addParameter('TE', 0.04, @(x) isscalar(x) && x > 0);
p.addParameter('Bins', 16, @(x) isscalar(x) && x >= 1);
p.addParameter('MaxNodes', 8, @(x) isscalar(x) && x >= 2);
p.addParameter('MaxRuns', Inf, @(x) isscalar(x) && x >= 1);
p.addParameter('Unit', 'atlas', @(x) ischar(x) || isstring(x));
p.addParameter('Atlas', '', @(x) ischar(x) || isstring(x));
p.addParameter('KernelModel', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('KernelObject', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('KernelCondition', '', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(source, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);
if ~isempty(opts.doverbose), verbose = logical(opts.doverbose); end

if exist('spm_dcm_estimate', 'file') ~= 2
    error('hrf_dcm:NoSPMDCM', 'spm_dcm_estimate (SPM12 DCM) is required on the path.');
end
if isempty(opts.Nodes), error('hrf_dcm:NoNodes', 'DCM needs an explicit small ''Nodes'' list.'); end
if isempty(local_clean(opts.Conditions)), error('hrf_dcm:NoConditions', 'DCM needs ''Conditions'' (the task inputs).'); end
try, spm('defaults', 'fMRI'); catch, end

% ---- extraction prep (RAW node timeseries) ------------------------------
if isstruct(source) && isfield(source, 'tsRuns')
    prep = source;
else
    prep = hrf_causality(source, 'ReturnData', true, 'Unit', opts.Unit, 'Nodes', opts.Nodes, ...
        'Atlas', opts.Atlas, 'KernelModel', opts.KernelModel, 'KernelObject', opts.KernelObject, ...
        'KernelCondition', opts.KernelCondition, 'MaxRuns', opts.MaxRuns, 'doverbose', verbose);
end
nodes = cellstr(string(prep.nodes));
n = numel(nodes);
if n > opts.MaxNodes
    error('hrf_dcm:TooManyNodes', '%d nodes requested; DCM caps at MaxNodes=%d. Narrow ''Nodes''.', n, opts.MaxNodes);
end
TR = prep.TR;

% ---- per-run DCM estimation --------------------------------------------
gcm = {}; runSubj = {};
for r = 1:numel(prep.tsRuns)
    y = prep.tsRuns{r};
    if any(~isfinite(y(:))) || any(std(y, 0, 1) == 0), continue; end
    U = local_build_inputs(prep.events_files{r}, size(y, 1), TR, opts.Conditions, opts.Bins);
    if isempty(U.u) || size(U.u, 2) == 0
        warning('hrf_dcm:NoInputs', 'No condition inputs for run %d; skipping.', r); continue
    end
    DCM = local_build_dcm(y, nodes, U, TR, opts);
    if verbose, fprintf('  DCM est run %d/%d (%d nodes, %d inputs, T=%d)...\n', r, numel(prep.tsRuns), n, size(U.u, 2), size(y, 1)); end
    try
        DCMe = spm_dcm_estimate(DCM);
    catch est_err
        warning('hrf_dcm:EstimateFailed', 'spm_dcm_estimate failed on run %d: %s', r, est_err.message); continue
    end
    gcm{end + 1} = DCMe; %#ok<AGROW>
    runSubj{end + 1} = char(string(prep.subjects{r})); %#ok<AGROW>
end
if isempty(gcm), error('hrf_dcm:NoDCM', 'No DCM was successfully estimated.'); end

% ---- collect Ep.A -> directed net, group across subjects ----------------
% SPM convention A(target,source); convert to i->j (A_mine = A_spm.').
[usubj, ~, sidx] = unique(runSubj, 'stable');
nSubj = numel(usubj);
net_subj = nan(n, n, nSubj);
A_subj = nan(n, n, nSubj);
for s = 1:nSubj
    runs_s = find(sidx == s);
    As = nan(n, n, numel(runs_s));
    for q = 1:numel(runs_s)
        As(:, :, q) = gcm{runs_s(q)}.Ep.A.';   % -> i->j convention
    end
    Am = mean(As, 3, 'omitnan');
    A_subj(:, :, s) = Am;
    net_subj(:, :, s) = Am - Am.';
end

R = struct('modes', {{'dcm'}}, 'nodes', string(nodes), 'subjects', {usubj(:)'}, ...
    'nsubj', nSubj, 'ndcm', numel(gcm), 'A_group', mean(A_subj, 3, 'omitnan'), 'gcm', {gcm(:)'});
R.dcm = local_group_stats(net_subj, nodes);
R.dcm.net_subj = net_subj;
if isfield(prep, 'dirs'), R.dirs = prep.dirs; end

if verbose
    [~, ix] = max(R.dcm.net_group(:) .* (R.dcm.p(:) < 0.05));
    [si, di] = ind2sub([n n], ix);
    fprintf('hrf_dcm: %d DCMs, %d subjects, %d nodes.\n', numel(gcm), nSubj, n);
    if isfinite(R.dcm.net_group(si, di)) && R.dcm.p(si, di) < 0.05
        fprintf('  strongest sig net flow: %s -> %s (net=%.3f Hz, p=%.2g)\n', nodes{si}, nodes{di}, R.dcm.net_group(si, di), R.dcm.p(si, di));
    else
        fprintf('  no net flow significant at p<.05 (group t, n=%d; consider PEB on R.gcm)\n', nSubj);
    end
end
end


% =========================================================================
function DCM = local_build_dcm(y, names, U, TR, opts)
[T, n] = size(y);
m = numel(U.name);
DCM = struct();
DCM.a = local_A(opts.A, n);
DCM.b = local_B(opts.B, n, m, U.name);
DCM.c = local_C(opts.C, n, m);
DCM.d = zeros(n, n, 0);
DCM.U = U;
DCM.Y.y = y;
DCM.Y.dt = TR;
DCM.Y.name = names(:)';
DCM.Y.X0 = ones(T, 1);
DCM.Y.Q = spm_Ce(repmat(T, 1, n));
DCM.v = T;
DCM.n = n;
DCM.TE = opts.TE;
DCM.delays = repmat(TR / 2, n, 1);
DCM.options = struct('nonlinear', 0, 'two_state', 0, 'stochastic', 0, ...
    'centre', 1, 'induced', 0, 'maxnodes', max(n, 8), 'nograph', 1);
end


function U = local_build_inputs(events_file, T, TR, conditions, bins)
% Box-function inputs at microtime dt = TR/bins, length T*bins, one column per
% requested condition (glob), built from the BIDS events onsets/durations.
U = struct('u', zeros(T * bins, 0), 'dt', TR / bins, 'name', {{}});
pats = local_clean(conditions);
if isempty(pats) || isempty(events_file) || exist(events_file, 'file') ~= 2, return; end
try
    E = readtable(events_file, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
catch
    return
end
v = E.Properties.VariableNames;
if ~all(ismember({'onset', 'trial_type'}, v)), return; end
onset = double(E.onset);
dur = zeros(size(onset)); if any(strcmp('duration', v)), dur = double(E.duration); end
tt = string(E.trial_type);
dt = TR / bins; nmt = T * bins;
u = []; names = {};
for k = 1:numel(pats)
    sel = local_glob(tt, pats{k});
    if ~any(sel), continue; end
    col = zeros(nmt, 1);
    rows = find(sel);
    for j = rows(:)'
        i0 = max(1, floor(onset(j) / dt) + 1);
        i1 = min(nmt, i0 + max(1, round(max(dur(j), dt) / dt)) - 1);
        col(i0:i1) = 1;
    end
    u = [u, col]; names{end + 1} = char(pats{k}); %#ok<AGROW>
end
U.u = u; U.name = names;
end


function A = local_A(spec, n)
if isnumeric(spec), A = double(spec ~= 0); return; end
switch lower(char(spec))
    case 'full', A = ones(n);
    case {'none', 'diagonal'}, A = eye(n);
    otherwise, A = ones(n);
end
end

function C = local_C(spec, n, m)
if isnumeric(spec), C = double(spec ~= 0); return; end
switch lower(char(spec))
    case 'all', C = ones(n, m);
    case 'none', C = zeros(n, m);
    otherwise, C = ones(n, m);
end
end

function B = local_B(spec, n, m, input_names)
if isnumeric(spec)
    if ndims(spec) == 3, B = double(spec ~= 0); else, B = repmat(double(spec ~= 0), 1, 1, m); end
    return
end
s = lower(char(spec));
switch s
    case {'none', ''}
        B = zeros(n, n, m);
    case 'all'
        B = ones(n, n, m);
    otherwise
        % a condition name: that input modulates all connections
        B = zeros(n, n, m);
        hit = find(strcmpi(input_names, s), 1);
        if ~isempty(hit), B(:, :, hit) = 1; end
end
end


function S = local_group_stats(net_subj, nodes)
[N, ~, nSubj] = size(net_subj);
mu = mean(net_subj, 3, 'omitnan');
t = nan(N); pv = nan(N);
if nSubj >= 2
    se = std(net_subj, 0, 3, 'omitnan') / sqrt(nSubj);
    t = mu ./ se; t(se == 0) = 0;
    pv = 2 * (1 - local_tcdf(abs(t), nSubj - 1));
end
pv(logical(eye(N))) = NaN;
S = struct('net_group', mu, 't', t, 'p', pv, 'p_fdr', local_fdr(pv), 'nodes', {nodes});
end


function pfdr = local_fdr(pv)
pfdr = nan(size(pv));
mask = isfinite(pv); ps = pv(mask); m = numel(ps);
if m == 0, return; end
[sp, ord] = sort(ps(:));
adj = min(1, flipud(cummin(flipud(sp .* m ./ (1:m)'))));
out = nan(m, 1); out(ord) = adj; pfdr(mask) = out;
end

function pcdf = local_tcdf(tval, df)
x = df ./ (df + tval .^ 2);
pcdf = 1 - 0.5 * betainc(x, df / 2, 0.5);
end

function tf = local_glob(tt, pat)
tt = cellstr(string(tt)); pat = strtrim(char(pat));
if any(pat == '*' | pat == '?')
    tf = ~cellfun('isempty', regexp(tt, ['^', regexptranslate('wildcard', pat), '$'], 'once'));
else
    tf = strcmpi(tt, pat);
end
tf = tf(:);
end

function c = local_clean(x)
c = cellstr(string(x));
c = c(~cellfun(@(s) isempty(strtrim(s)), c));
end
