function G = hrf_group_stats(data, varargin)
%HRF_GROUP_STATS Group inference with permutation correction (shared engine).
%
% One canonical group-level test used across the HRF causality/mediation/curve
% API, so every summary shares the same small-n-appropriate inference instead
% of ad-hoc parametric t + FDR. Given a per-subject stack it returns the group
% estimate, t, uncorrected p, and a CORRECTED p / significance mask by the
% chosen method:
%
%   'permutation' (default) - sign-flip (one-sample) or label (two-sample)
%                  permutation, FWER-controlled over the tested family via the
%                  maximum |t| statistic (Nichols & Holmes 2002). Exact when
%                  2^n is enumerated (one-sample, n<=13). The whole array is
%                  permuted per subject, so correlations across cells are
%                  preserved under the null -- the right choice at n~7-11.
%   'cluster'    - temporal cluster-mass permutation (Maris & Oostenveld 2007)
%                  along 'ClusterDim' (e.g. lag); needs a 2-D cell array.
%   'fdr'        - Benjamini-Hochberg across the tested family.
%   'none'       - uncorrected.
%
% :Usage:
% ::
%     G = hrf_group_stats(net_subj, 'Mask', ~eye(N))          % connectivity edges
%     G = hrf_group_stats(curve_subj, 'ClusterDim', 1, 'Correction','cluster')
%     G = hrf_group_stats(x, 'Design','twosample', 'Group', grp)
%
% :Inputs:
%   **data:** numeric array [d1 x ... x dk x nSubj]; the LAST dimension indexes
%             subjects. Cells are the tested family (e.g. edges, term x lag).
%
% :Optional Inputs:
%   **'Design':**   'onesample' (default; test mean~=0 by sign-flip) or
%                   'twosample' (test group difference by label permutation).
%   **'Group':**    for twosample: [nSubj] vector of two group labels.
%   **'Correction':** 'permutation'|'cluster'|'fdr'|'fdr_all'|'none'. Default
%                   'permutation'.
%   **'Threshold':** corrected-p / FWER threshold. Default 0.05.
%   **'Nperm':**    permutations when the exact set is too large. Default 5000.
%   **'Mask':**     logical over [d1..dk] of cells to include in the family
%                   (e.g. off-diagonal). Default: cells finite in >=1 subject.
%   **'ClusterDim':** dimension for 'cluster' contiguity (2-D data only).
%   **'ClusterFormP':** cluster-forming p. Default 0.05.
%
% :Output:
%   **G:** struct with .est, .t, .p, .p_corr, .sig (all size [d1..dk]), plus
%          .n, .correction, .design, .mask.
%
% See also: hrf_causality_analyze, hrf_causality_contrast, hrf_dcm,
%           hrf_mediation_analyze, hrf_animate_wordcloud.

p = inputParser;
p.addRequired('data', @isnumeric);
p.addParameter('Design', 'onesample', @(x) ischar(x) || isstring(x));
p.addParameter('Group', [], @(x) isempty(x) || isvector(x) || iscell(x));
p.addParameter('Correction', 'permutation', @(x) ischar(x) || isstring(x));
p.addParameter('Threshold', 0.05, @(x) isscalar(x) && x > 0 && x <= 1);
p.addParameter('Nperm', 5000, @(x) isscalar(x) && x >= 100);
p.addParameter('Mask', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.addParameter('ClusterDim', [], @(x) isempty(x) || isscalar(x));
p.addParameter('ClusterFormP', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.parse(data, varargin{:});
opts = p.Results;

sz = size(data);
nsub = sz(end);
dims = sz(1:end - 1);
if isscalar(dims), dims = [dims, 1]; end
ncell = prod(dims);
X = reshape(data, ncell, nsub);              % [cells x nSubj]
design = lower(char(opts.Design));
corr = lower(char(opts.Correction));

% family mask (over cells)
if isempty(opts.Mask)
    mask = any(isfinite(X), 2);
else
    mask = logical(opts.Mask(:));
end

% two-sample group labels
if strcmp(design, 'twosample')
    g = local_two_groups(opts.Group, nsub);
else
    g = [];
end

% observed statistic
[t_obs, est] = local_stat(X, design, g);
[pobs] = local_param_p(t_obs, X, design, g);

% correction
switch corr
    case {'permutation', 'perm', 'maxt'}
        pcorr = local_perm_fwe(X, mask, design, g, t_obs, opts.Nperm);
        sig = pcorr <= opts.Threshold;
    case {'cluster', 'tcluster'}
        [pcorr, sig] = local_perm_cluster(X, mask, dims, design, g, opts);
    case {'fdr_all', 'fdrall', 'fdr'}
        pcorr = nan(ncell, 1);
        pcorr(mask) = local_bh(pobs(mask));
        sig = pcorr <= opts.Threshold;
    otherwise                                    % 'none'
        pcorr = pobs;
        sig = pobs < opts.Threshold;
end
sig(~mask) = false;
pobs(~mask) = NaN; pcorr(~mask) = NaN;

G = struct();
G.est = reshape(est, dims);
G.t = reshape(t_obs, dims);
G.p = reshape(pobs, dims);
G.p_corr = reshape(pcorr, dims);
G.sig = reshape(sig, dims);
G.n = nsub;
G.correction = corr;
G.design = design;
G.mask = reshape(mask, dims);
end


% =========================================================================
function [t, est] = local_stat(X, design, g)
% t and estimate per cell (rows of X = cells, cols = subjects).
if strcmp(design, 'twosample')
    A = X(:, g == 1); B = X(:, g == 2);
    mA = mean(A, 2, 'omitnan'); mB = mean(B, 2, 'omitnan');
    vA = var(A, 0, 2, 'omitnan'); vB = var(B, 0, 2, 'omitnan');
    nA = sum(isfinite(A), 2); nB = sum(isfinite(B), 2);
    se = sqrt(vA ./ max(nA, 1) + vB ./ max(nB, 1));
    t = (mA - mB) ./ se; t(se == 0) = 0;
    est = mA - mB;
else
    m = mean(X, 2, 'omitnan');
    se = std(X, 0, 2, 'omitnan') ./ sqrt(sum(isfinite(X), 2));
    t = m ./ se; t(se == 0) = 0;
    est = m;
end
t(~isfinite(t)) = 0;
end


function pp = local_param_p(t, X, design, g)
if strcmp(design, 'twosample')
    A = X(:, g == 1); B = X(:, g == 2);
    vA = var(A, 0, 2, 'omitnan'); vB = var(B, 0, 2, 'omitnan');
    nA = sum(isfinite(A), 2); nB = sum(isfinite(B), 2);
    sa = vA ./ max(nA, 1); sb = vB ./ max(nB, 1);
    df = (sa + sb) .^ 2 ./ (sa .^ 2 ./ max(nA - 1, 1) + sb .^ 2 ./ max(nB - 1, 1));
    df(~isfinite(df) | df < 1) = 1;
else
    df = sum(isfinite(X), 2) - 1; df(df < 1) = 1;
end
pp = 2 * (1 - local_tcdf(abs(t), df));
end


function pcorr = local_perm_fwe(X, mask, design, g, t_obs, nperm)
% FWER over the masked family via max|t| null.
absobs = abs(t_obs);
ge = zeros(size(t_obs));
perms = local_perm_set(X, design, g, nperm);
for k = 1:numel(perms)
    tp = abs(local_stat_perm(X, design, g, perms{k}));
    mx = max(tp(mask), [], 'omitnan');
    if isempty(mx) || ~isfinite(mx), mx = 0; end
    ge = ge + (mx >= absobs);
end
pcorr = ge / numel(perms);
end


function [pcorr, sig] = local_perm_cluster(X, mask, dims, design, g, opts)
% Temporal cluster-mass along ClusterDim (2-D data). Falls back to FWER if no
% valid cluster dim.
if isempty(opts.ClusterDim) || numel(dims) ~= 2
    pcorr = local_perm_fwe(X, mask, design, g, local_stat(X, design, g), opts.Nperm);
    sig = pcorr <= opts.Threshold; return
end
cdim = opts.ClusterDim;
n = size(X, 2);
tcrit = local_tinv(1 - opts.ClusterFormP / 2, max(n - 1, 1));
maskM = reshape(mask, dims);
tobsM = reshape(local_stat(X, design, g), dims);
if cdim == 2, tobsM = tobsM'; maskM = maskM'; end   % put cluster dim as rows
oc = local_clusters(tobsM, maskM, tcrit);
perms = local_perm_set(X, design, g, opts.Nperm);
maxmass = zeros(numel(perms), 1);
for k = 1:numel(perms)
    tm = reshape(local_stat_perm(X, design, g, perms{k}), dims);
    if cdim == 2, tm = tm'; end
    cl = local_clusters(tm, maskM, tcrit);
    if ~isempty(cl), maxmass(k) = max([cl.mass]); end
end
pM = nan(size(tobsM)); sM = false(size(tobsM));
for c = 1:numel(oc)
    pc = mean(maxmass >= oc(c).mass);
    pM(oc(c).cells) = pc; sM(oc(c).cells) = pc <= opts.Threshold;
end
if cdim == 2, pM = pM'; sM = sM'; end
pcorr = pM(:); sig = sM(:);
end


function cl = local_clusters(tmat, maskM, tcrit)
[L, T] = size(tmat);
cl = struct('cells', {}, 'mass', {});
for j = 1:T
    col = tmat(:, j); mk = maskM(:, j);
    for s = [1 -1]
        supra = (s * col > tcrit) & isfinite(col) & mk;
        d = diff([false; supra; false]);
        starts = find(d == 1); stops = find(d == -1) - 1;
        for b = 1:numel(starts)
            idx = (starts(b):stops(b))';
            cl(end + 1) = struct('cells', sub2ind([L T], idx, repmat(j, numel(idx), 1)), ...
                'mass', sum(abs(col(idx)))); %#ok<AGROW>
        end
    end
end
end


function perms = local_perm_set(X, design, g, nperm)
n = size(X, 2);
if strcmp(design, 'twosample')
    perms = cell(1, nperm);
    perms{1} = g;                                 % identity first
    for k = 2:nperm, perms{k} = g(randperm(n)); end
else
    if n <= 13
        m = 2 ^ n; perms = cell(1, m);
        for f = 1:m
            s = 1 - 2 * bitget(f - 1, 1:n);
            perms{f} = s;
        end
    else
        perms = cell(1, nperm);
        perms{1} = ones(1, n);
        for k = 2:nperm, perms{k} = 2 * (rand(1, n) > 0.5) - 1; end
    end
end
end


function t = local_stat_perm(X, design, g, perm)
if strcmp(design, 'twosample')
    t = local_stat(X, design, perm);              % perm = shuffled labels
else
    t = local_stat(X .* perm, design, g);          % perm = sign vector (1 x n)
end
end


function padj = local_bh(pv)
pv = pv(:); m = numel(pv);
[sp, ord] = sort(pv);
adj = min(1, flipud(cummin(flipud(sp .* m ./ (1:m)'))));
padj = nan(m, 1); padj(ord) = adj;
end


function g = local_two_groups(grp, nsub)
if isempty(grp), error('hrf_group_stats:NoGroup', 'Design=twosample needs ''Group'' labels.'); end
if iscell(grp) || isstring(grp)
    [~, ~, g] = unique(cellstr(string(grp)), 'stable');
else
    [~, ~, g] = unique(grp(:), 'stable');
end
if numel(g) ~= nsub, error('hrf_group_stats:GroupLen', 'Group has %d labels, need %d.', numel(g), nsub); end
if max(g) ~= 2, error('hrf_group_stats:TwoGroups', 'Design=twosample needs exactly 2 groups.'); end
g = g(:)';
end


function pcdf = local_tcdf(tval, df)
x = df ./ (df + tval .^ 2);
pcdf = 1 - 0.5 * betainc(x, df ./ 2, 0.5);
end


function t = local_tinv(pp, df)
lo = 0; hi = 100;
for it = 1:60
    mid = (lo + hi) / 2;
    if local_tcdf(mid, df) < pp, lo = mid; else, hi = mid; end
end
t = (lo + hi) / 2;
end
