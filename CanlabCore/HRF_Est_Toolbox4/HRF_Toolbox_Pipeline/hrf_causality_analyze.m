function R = hrf_causality_analyze(tsRuns, kernels, varargin)
%HRF_CAUSALITY_ANALYZE Deconvolve -> Granger -> group stats (the compute core).
%
% Computation half of the HRF-informed causality pipeline: takes per-run node
% timeseries plus per-node HRF kernels, deconvolves each run, optionally
% removes the task-evoked component, runs Granger causality per subject
% (pooling that subject's runs), and aggregates the directed net-flow across
% subjects with a one-sample t-test (BH-FDR over the off-diagonal).
%
% The IO half -- extracting signature/region timeseries from 4-D BOLD and
% pulling the sFIR kernels from the score CSVs -- lives in hrf_causality,
% which calls this. Kept separate so the inference is unit-testable without
% touching disk.
%
% Usage
% -----
%   R = hrf_causality_analyze(tsRuns, kernels, 'Nodes', names)
%   R = hrf_causality_analyze(tsRuns, kernels, 'Subjects', subjvec, ...
%           'EvokedMode','both', 'Confounds', confRuns)
%
% Inputs
% ------
%   tsRuns  - cell{1..nRun} of [T x N] node timeseries (one matrix per run;
%             same N columns/order in every run).
%   kernels - [L x N] per-node HRF kernels (column n deconvolves node n), or
%             [L x 1] to use one kernel for all nodes. Sampled at the data TR.
%
% Optional (name-value)
% ---------------------
%   'Subjects'   - [nRun x 1] subject id per run (numeric or cellstr). Runs
%                  with the same id are pooled for that subject's GC. Default:
%                  all runs = one subject.
%   'Nodes'      - N node names. Default n1..nN.
%   'EvokedMode' - 'both' (default) | 'remove' | 'keep'. 'remove' regresses
%                  Confounds out of the neural proxy before GC (endogenous
%                  coupling); 'keep' uses the full proxy; 'both' computes each
%                  and returns them side by side.
%   'Confounds'  - cell{1..nRun} of [T x q] regressors removed when the mode
%                  is 'remove'/'both'. For proxy-space removal these should be
%                  NEURAL-level (UNconvolved) task regressors. Default {} =>
%                  'remove' falls back to removing the per-run mean only.
%   'Segments'   - cell{1..nRun} of logical [T x 1] masks. When given, each
%                  run's proxy is split into the contiguous TRUE blocks of its
%                  mask and each block becomes a separate GC realization, so
%                  Granger is computed ONLY within that condition's blocks
%                  (no autoregression across block gaps). This is how
%                  condition-specific causality is done (e.g. rest_stim vs
%                  nback blocks). Default {} => whole run.
%   'MinSegLen'  - drop condition blocks shorter than this many samples.
%                  Default 20.
%   'DeconvMethod','Lambda' - passed to hrf_deconv_timeseries. Default ridge.
%   'Conditional','Order','MaxOrder','Nperm' - passed to hrf_granger_causality.
%   'Verbose'/'doverbose' - print a summary (default true).
%
% Output
% ------
%   R - struct:
%     .modes        cellstr of evoked modes computed
%     .(mode)       per mode: struct with
%         .net_group [N x N]  mean across subjects of net flow (gc-gc')
%         .t,.p      [N x N]  one-sample t and p across subjects (diag NaN)
%         .p_fdr     [N x N]  BH-FDR-corrected p over the off-diagonal
%         .net_subj  [N x N x nSubj]  per-subject net-flow matrices
%     .nodes, .subjects, .nsubj, .nrun, .order
%
% See also: hrf_deconv_timeseries, hrf_granger_causality, hrf_causality.

p = inputParser;
p.addRequired('tsRuns', @(x) iscell(x) && ~isempty(x));
p.addRequired('kernels', @(x) isnumeric(x) && ~isempty(x));
p.addParameter('Subjects', [], @(x) isempty(x) || isvector(x) || iscell(x));
p.addParameter('Nodes', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('EvokedMode', 'both', @(x) ischar(x) || isstring(x));
p.addParameter('Confounds', {}, @(x) iscell(x));
p.addParameter('Segments', {}, @(x) iscell(x));
p.addParameter('MinSegLen', 20, @(x) isscalar(x) && x >= 4);
p.addParameter('DeconvMethod', 'ridge', @(x) ischar(x) || isstring(x));
p.addParameter('Lambda', [], @(x) isempty(x) || isscalar(x));
p.addParameter('Conditional', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Order', 'bic', @(x) (ischar(x) || isstring(x)) || isscalar(x));
p.addParameter('MaxOrder', 10, @(x) isscalar(x) && x >= 1);
p.addParameter('Nperm', 0, @(x) isscalar(x) && x >= 0);
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(tsRuns, kernels, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);
if ~isempty(opts.doverbose), verbose = logical(opts.doverbose); end

nRun = numel(tsRuns);
N = size(tsRuns{1}, 2);
nodes = local_node_names(opts.Nodes, N);

subj = opts.Subjects;
if isempty(subj), subj = ones(nRun, 1); end
subj = local_to_subject_ids(subj, nRun);
[usubj, ~, sidx] = unique(subj, 'stable');
nSubj = numel(usubj);

modes = local_modes(opts.EvokedMode);
have_conf = ~isempty(opts.Confounds);
have_seg = ~isempty(opts.Segments);
minseglen = opts.MinSegLen;

% Deconvolve every run once (kernel is mode-independent).
proxy = cell(1, nRun);
for r = 1:nRun
    proxy{r} = hrf_deconv_timeseries(tsRuns{r}, kernels, ...
        'Method', opts.DeconvMethod, 'Lambda', opts.Lambda);
end

gargs = {'Nodes', nodes, 'Conditional', opts.Conditional, 'Order', opts.Order, ...
    'MaxOrder', opts.MaxOrder, 'Nperm', opts.Nperm, 'doverbose', false};

R = struct('modes', {modes}, 'nodes', string(nodes), 'nsubj', nSubj, 'nrun', nRun, 'order', NaN);
R.subjects = usubj(:)';
for mi = 1:numel(modes)
    mode = modes{mi};
    net_subj = nan(N, N, nSubj);
    order_used = nan(1, nSubj);
    for s = 1:nSubj
        runs_s = find(sidx == s);
        Xs = {};
        for q = 1:numel(runs_s)
            rr = runs_s(q);
            Xr = proxy{rr};
            if ~strcmp(mode, 'keep')
                conf = [];
                if have_conf && rr <= numel(opts.Confounds), conf = opts.Confounds{rr}; end
                Xr = local_remove_evoked(Xr, conf);
            end
            segmask = [];
            if have_seg && rr <= numel(opts.Segments), segmask = opts.Segments{rr}; end
            Xs = [Xs, local_segment(Xr, segmask, minseglen)]; %#ok<AGROW>
        end
        if isempty(Xs), continue; end       % no usable blocks for this subject
        try
            G = hrf_granger_causality(Xs, gargs{:});
            net_subj(:, :, s) = G.net;
            order_used(s) = G.order;
        catch ge
            warning('hrf_causality_analyze:SubjectGCFailed', ...
                'Subject %s GC failed (%s); skipping.', usubj{s}, ge.message);
        end
    end
    R.(mode) = local_group_stats(net_subj, nodes);
    R.(mode).net_subj = net_subj;
    if mi == 1, R.order = round(median(order_used, 'omitnan')); end
end

if verbose
    fprintf('hrf_causality_analyze: N=%d nodes, %d run(s), %d subject(s); modes: %s\n', ...
        N, nRun, nSubj, strjoin(modes, ', '));
    for mi = 1:numel(modes)
        S = R.(modes{mi});
        [~, ix] = max(S.net_group(:) .* (S.p(:) < 0.05));
        [si, di] = ind2sub([N N], ix);
        if isfinite(S.net_group(si, di)) && S.p(si, di) < 0.05
            fprintf('  [%s] strongest sig net flow: %s -> %s (net=%.3f, p=%.2g, p_fdr=%.2g)\n', ...
                modes{mi}, nodes{si}, nodes{di}, S.net_group(si, di), S.p(si, di), S.p_fdr(si, di));
        else
            fprintf('  [%s] no net flow significant at p<.05\n', modes{mi});
        end
    end
end
end


% =========================================================================
function S = local_group_stats(net_subj, nodes)
[N, ~, nSubj] = size(net_subj);
mu = mean(net_subj, 3, 'omitnan');
t = nan(N); pv = nan(N);
if nSubj >= 2
    sd = std(net_subj, 0, 3, 'omitnan');
    se = sd / sqrt(nSubj);
    t = mu ./ se;
    t(se == 0) = 0;
    pv = 2 * (1 - local_tcdf(abs(t), nSubj - 1));
end
pv(logical(eye(N))) = NaN;
S = struct('net_group', mu, 't', t, 'p', pv, 'p_fdr', local_fdr(pv), 'nodes', {nodes});
end


function segs = local_segment(X, mask, minlen)
% Split X [T x N] into the contiguous blocks where mask is true (each >=
% minlen rows). Empty mask => the whole run as a single segment. Each block
% becomes a separate GC realization so no autoregression crosses a gap.
if isempty(mask)
    segs = {X};
    return
end
mask = logical(mask(:));
T = size(X, 1);
if numel(mask) ~= T
    L = min(numel(mask), T); mask = mask(1:L); X = X(1:L, :);
end
segs = {};
d = diff([false; mask; false]);
starts = find(d == 1);
stops = find(d == -1) - 1;
for b = 1:numel(starts)
    if stops(b) - starts(b) + 1 >= minlen
        segs{end + 1} = X(starts(b):stops(b), :); %#ok<AGROW>
    end
end
end


function Xc = local_remove_evoked(X, conf)
% Remove the evoked component from each column: regress out [1, conf] (or
% just the mean if no confounds), keep the residual.
T = size(X, 1);
D = ones(T, 1);
if ~isempty(conf)
    if size(conf, 1) ~= T
        error('hrf_causality_analyze:ConfRows', 'Confound rows (%d) ~= run length (%d).', size(conf, 1), T);
    end
    D = [D, conf];
end
Xc = X - D * (D \ X);
end


function pfdr = local_fdr(pv)
% Benjamini-Hochberg over the finite off-diagonal entries.
pfdr = nan(size(pv));
mask = isfinite(pv);
ps = pv(mask);
m = numel(ps);
if m == 0, return; end
[sp, ord] = sort(ps(:));
adj = sp .* m ./ (1:m)';
adj = min(1, flipud(cummin(flipud(adj))));
out = nan(m, 1); out(ord) = adj;
pfdr(mask) = out;
end


function modes = local_modes(em)
switch lower(strtrim(char(em)))
    case 'both',   modes = {'remove', 'keep'};
    case 'remove', modes = {'remove'};
    case 'keep',   modes = {'keep'};
    otherwise, error('hrf_causality_analyze:EvokedMode', 'EvokedMode must be both|remove|keep.');
end
end


function names = local_node_names(nodes, N)
if isempty(nodes)
    names = arrayfun(@(i) sprintf('n%d', i), 1:N, 'uni', 0);
else
    names = cellstr(string(nodes));
    if numel(names) ~= N
        error('hrf_causality_analyze:Nodes', 'Nodes has %d names but N=%d.', numel(names), N);
    end
end
end


function ids = local_to_subject_ids(subj, nRun)
if iscell(subj) || isstring(subj)
    ids = cellstr(string(subj));
else
    ids = cellstr(string(subj(:)));
end
if numel(ids) ~= nRun
    error('hrf_causality_analyze:Subjects', 'Subjects has %d entries but %d runs.', numel(ids), nRun);
end
end


function pcdf = local_tcdf(tval, df)
% Student-t CDF via regularized incomplete beta (no Stats toolbox needed).
x = df ./ (df + tval .^ 2);
ib = 0.5 * betainc(x, df / 2, 0.5);
pcdf = 1 - ib;             % for tval >= 0
end
