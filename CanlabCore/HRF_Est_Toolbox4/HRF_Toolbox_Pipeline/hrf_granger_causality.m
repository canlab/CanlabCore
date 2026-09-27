function G = hrf_granger_causality(X, varargin)
%HRF_GRANGER_CAUSALITY Directed Granger causality on (deconvolved) timeseries.
%
% Pooled-VAR Granger causality for a set of signals (regions/signatures),
% intended to run on the HRF-DECONVOLVED proxies from hrf_deconv_timeseries
% so the directionality reflects neural lead/lag rather than regional
% hemodynamic latency. Handles multiple runs/subjects as separate
% realizations of a common VAR (no autoregression across run boundaries),
% chooses model order by BIC, and reports both the directed influence matrix
% and its asymmetry (net flow).
%
% Granger logic: i "Granger-causes" j if the past of i improves prediction
% of j beyond j's own past (pairwise), or beyond j's own past AND every
% other node's past (conditional / multivariate -- removes influence routed
% through or shared with third nodes).
%
% Usage
% -----
%   G = hrf_granger_causality(X)                       % one run, [T x N]
%   G = hrf_granger_causality({X1, X2, ...})           % multi-run/subject
%   G = hrf_granger_causality(X, 'Conditional', true, 'Nodes', names)
%   G = hrf_granger_causality(X, 'Order', 4, 'Nperm', 1000)
%
% Inputs
% ------
%   X - [T x N] timeseries (columns = nodes), OR a cell array of such
%       matrices (one per run/subject), each with the SAME N columns.
%
% Optional (name-value)
% ---------------------
%   'Nodes'       - cellstr/string of N node names (for labeling). Default {}.
%   'Order'       - VAR model order: integer, or 'bic' (default) / 'aic' to
%                   pick over 1..MaxOrder.
%   'MaxOrder'    - max order searched when Order is 'bic'/'aic'. Default 10.
%   'Conditional' - true => conditional MVGC (condition on all other nodes);
%                   false (default) => pairwise GC.
%   'Detrend'     - linearly detrend each run column first. Default true.
%   'Zscore'      - z-score each run column first (recommended; makes runs
%                   commensurate before pooling). Default true.
%   'Nperm'       - 0 (default) => parametric F-test p-values. >0 => add a
%                   null by circularly shifting the source within each run
%                   (preserves each series' autocorrelation), p_perm column.
%   'Verbose'/'doverbose' - print a short summary. Default true.
%
% Outputs
% -------
%   G - struct:
%     .gc        [N x N]  directed influence; gc(i,j) = i -> j (log-ratio of
%                         restricted/full residual variance, >= 0).
%     .fstat     [N x N]  F statistic for each i -> j.
%     .pval      [N x N]  parametric F-test p-value (diag = NaN).
%     .pval_perm [N x N]  permutation p-value (only if Nperm > 0).
%     .net       [N x N]  gc - gc' (positive = net flow i -> j).
%     .order      scalar  VAR order used.
%     .conditional logical
%     .n_eff      scalar  pooled samples used.
%     .nodes      string  node names.
%
% Notes / caveats
% ---------------
% * Run this on DECONVOLVED signals. On raw BOLD, regional HRF differences
%   can reverse the inferred direction (David et al. 2008). See
%   hrf_deconv_timeseries.
% * Task-evoked common drive can manufacture spurious GC. Remove the evoked
%   task mean (or analyze within a condition) before calling this if your
%   question is about endogenous coupling.
% * GC is a statistical (predictive) notion of directed influence, not proof
%   of a physical mechanism. Treat results as hypotheses; corroborate with a
%   generative model (DCM) for strong claims.
%
% Example
% -------
%   T = 600; a = randn(T,1); b = zeros(T,1);
%   for t = 3:T, b(t) = 0.5*b(t-1) + 0.6*a(t-2) + 0.3*randn; end  % a -> b
%   G = hrf_granger_causality([a b], 'Nodes', {'A','B'});
%   G.net   % G.net(1,2) > 0  => net flow A -> B
%
% See also: hrf_deconv_timeseries, hrf_causality, hrf_misspec_metrics.

p = inputParser;
p.addRequired('X', @(x) isnumeric(x) || iscell(x));
p.addParameter('Nodes', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Order', 'bic', @(x) (ischar(x) || isstring(x)) || (isscalar(x) && x >= 1));
p.addParameter('MaxOrder', 10, @(x) isscalar(x) && x >= 1);
p.addParameter('Conditional', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Detrend', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Zscore', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Nperm', 0, @(x) isscalar(x) && x >= 0);
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(X, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);
if ~isempty(opts.doverbose), verbose = logical(opts.doverbose); end

% ---- normalize to a cell of runs, preprocess ----------------------------
if isnumeric(X), runs = {X}; else, runs = X(:)'; end
N = size(runs{1}, 2);
for r = 1:numel(runs)
    if size(runs{r}, 2) ~= N
        error('hrf_granger_causality:Ncols', 'All runs must have the same N columns.');
    end
    Xr = double(runs{r});
    if logical(opts.Detrend), Xr = detrend(Xr, 1); end
    if logical(opts.Zscore),  Xr = local_zscore(Xr); end
    runs{r} = Xr;
end

nodes = local_node_names(opts.Nodes, N);
conditional = logical(opts.Conditional);

% ---- model order --------------------------------------------------------
if ischar(opts.Order) || isstring(opts.Order)
    porder = local_select_order(runs, opts.MaxOrder, lower(char(opts.Order)));
else
    porder = round(opts.Order);
end

% ---- pooled target + lag tensor -----------------------------------------
[Y, Lag, run_rows] = local_build_lags(runs, porder);   % Y [M x N], Lag [M x N x p]
M = size(Y, 1);
if M <= N * porder + 2
    error('hrf_granger_causality:TooFewSamples', ...
        'Only %d pooled samples for order %d, N=%d -- increase data or lower Order.', M, porder, N);
end

% ---- GC matrix ----------------------------------------------------------
[gc, fstat, pval, df2] = local_gc_matrix(Y, Lag, porder, conditional);

G = struct();
G.gc = gc; G.fstat = fstat; G.pval = pval;
G.net = gc - gc';
G.order = porder; G.conditional = conditional;
G.n_eff = M; G.nodes = string(nodes);
G.df = [porder, df2];

% ---- permutation null (optional) ----------------------------------------
if opts.Nperm > 0
    G.pval_perm = local_perm_pvals(runs, porder, conditional, gc, opts.Nperm, run_rows);
end

if verbose
    fprintf('hrf_granger_causality: N=%d nodes, %d run(s), order=%d (%s), %s, M=%d samples\n', ...
        N, numel(runs), porder, local_order_src(opts.Order), ...
        local_tern(conditional, 'conditional MVGC', 'pairwise'), M);
    [~, ix] = max(G.net(:));
    [si, di] = ind2sub([N N], ix);
    fprintf('  strongest net flow: %s -> %s (net=%.3f, p=%.2g)\n', ...
        nodes{si}, nodes{di}, G.net(si, di), pval(si, di));
end
end


% =========================================================================
function [gc, fstat, pval, df2] = local_gc_matrix(Y, Lag, p, conditional)
N = size(Y, 2);
M = size(Y, 1);
gc = zeros(N); fstat = zeros(N); pval = nan(N);
intercept = ones(M, 1);

for j = 1:N
    yj = Y(:, j);
    ownj = local_lagblock(Lag, j);           % M x p

    if conditional
        % full = intercept + own + ALL other nodes' lags
        others = setdiff(1:N, j);
        fullblk = [intercept, ownj, local_lagblock(Lag, others)];
        rss_full = local_rss(yj, fullblk);
        df_full = size(fullblk, 2);
        df2 = M - df_full;
        for i = others
            % restricted = full minus source i's lag columns
            keep = setdiff(others, i);
            restr = [intercept, ownj, local_lagblock(Lag, keep)];
            rss_r = local_rss(yj, restr);
            [gc(i, j), fstat(i, j), pval(i, j)] = local_gc_stat(rss_r, rss_full, p, df2);
        end
    else
        restr = [intercept, ownj];
        rss_r = local_rss(yj, restr);
        df_full = size(restr, 2) + p;
        df2 = M - df_full;
        for i = 1:N
            if i == j, continue; end
            fullblk = [intercept, ownj, local_lagblock(Lag, i)];
            rss_f = local_rss(yj, fullblk);
            [gc(i, j), fstat(i, j), pval(i, j)] = local_gc_stat(rss_r, rss_f, p, df2);
        end
    end
end
end


function [g, F, pv] = local_gc_stat(rss_r, rss_f, p, df2)
rss_f = max(rss_f, realmin);
g = max(log(rss_r / rss_f), 0);
F = ((rss_r - rss_f) / p) / (rss_f / df2);
F = max(F, 0);
pv = 1 - local_fcdf(F, p, df2);
end


function B = local_lagblock(Lag, cols)
% Flatten Lag(:, cols, :) into [M x (numel(cols)*p)].
[M, ~, p] = size(Lag);
cols = cols(:)';
B = zeros(M, numel(cols) * p);
c = 0;
for ci = cols
    B(:, c + (1:p)) = reshape(Lag(:, ci, :), M, p);
    c = c + p;
end
end


function rss = local_rss(y, D)
b = D \ y;
res = y - D * b;
rss = sum(res .^ 2);
end


function [Y, Lag, run_rows] = local_build_lags(runs, p)
% Build pooled target Y [M x N] and lag tensor Lag [M x N x p], aligned so
% row m of Y corresponds to lags Lag(m, :, 1..p). Lags never cross run edges.
N = size(runs{1}, 2);
blocks = cellfun(@(r) max(size(r, 1) - p, 0), runs);
M = sum(blocks);
Y = zeros(M, N);
Lag = zeros(M, N, p);
run_rows = cell(1, numel(runs));
off = 0;
for r = 1:numel(runs)
    Xr = runs{r}; Tr = size(Xr, 1);
    if Tr <= p, run_rows{r} = []; continue; end
    rows = (off + 1):(off + (Tr - p));
    Y(rows, :) = Xr(p + 1:Tr, :);
    for k = 1:p
        Lag(rows, :, k) = Xr(p + 1 - k:Tr - k, :);
    end
    run_rows{r} = rows;
    off = off + (Tr - p);
end
end


function porder = local_select_order(runs, maxord, crit)
% Pick VAR order by AIC/BIC on the pooled N-variate model. Fit each candidate
% on the SAME sample (trimmed to maxord) so scores are comparable.
N = size(runs{1}, 2);
best = inf; porder = 1;
for p = 1:maxord
    [Y, Lag] = local_build_lags(runs, p);
    M = size(Y, 1);
    if M <= N * p + 2, break; end
    intercept = ones(M, 1);
    resid = zeros(M, N);
    for j = 1:N
        D = [intercept, local_lagblock(Lag, 1:N)];
        b = D \ Y(:, j);
        resid(:, j) = Y(:, j) - D * b;
    end
    Sigma = (resid' * resid) / M;
    ld = local_logdet(Sigma);
    nparams = N * (N * p + 1);
    switch crit
        case 'aic', score = M * ld + 2 * nparams;
        otherwise,  score = M * ld + nparams * log(M);   % bic
    end
    if score < best, best = score; porder = p; end
end
end


function pv = local_perm_pvals(runs, p, conditional, gc_obs, nperm, ~)
% Null by circularly shifting each SOURCE series within each run (preserves
% per-series autocorrelation, breaks directed coupling). Source-specific
% shift => recompute that source's contribution.
N = size(runs{1}, 2);
ge = zeros(N);   % count of null >= observed
for b = 1:nperm
    sruns = runs;
    for r = 1:numel(sruns)
        Tr = size(sruns{r}, 1);
        if Tr < 4, continue; end
        sh = 1 + mod(b * 7 + r * 13, Tr - 1);   % deterministic, varied shift
        sruns{r} = circshift(sruns{r}, sh, 1);
    end
    [Yp, Lagp] = local_build_lags(sruns, p);
    gcp = local_gc_matrix(Yp, Lagp, p, conditional);
    ge = ge + (gcp >= gc_obs);
end
pv = (1 + ge) / (1 + nperm);
pv(logical(eye(N))) = NaN;
end


% =========================================================================
function Z = local_zscore(X)
mu = mean(X, 1); sd = std(X, 0, 1); sd(sd == 0) = 1;
Z = (X - mu) ./ sd;
end

function names = local_node_names(nodes, N)
if isempty(nodes)
    names = arrayfun(@(i) sprintf('n%d', i), 1:N, 'uni', 0);
else
    names = cellstr(string(nodes));
    if numel(names) ~= N
        error('hrf_granger_causality:Nodes', 'Nodes has %d names but N=%d.', numel(names), N);
    end
end
end

function ld = local_logdet(S)
[Rc, fl] = chol(S);
if fl == 0, ld = 2 * sum(log(diag(Rc))); else, ld = log(max(det(S), realmin)); end
end

function pcdf = local_fcdf(F, d1, d2)
% F CDF via regularized incomplete beta (avoids a Stats-toolbox dependency).
x = d1 * F / (d1 * F + d2);
pcdf = betainc(x, d1 / 2, d2 / 2);
end

function s = local_order_src(o)
if ischar(o) || isstring(o), s = char(o); else, s = 'fixed'; end
end

function s = local_tern(c, a, b)
if c, s = a; else, s = b; end
end
