function C = hrf_causality_contrast(R1, R2, varargin)
%HRF_CAUSALITY_CONTRAST Contrast directed net-flow between two conditions.
%
% Compares the per-subject directed net-flow of two hrf_causality results and
% returns the group-level difference (R1 - R2) with statistics, in the same
% struct shape hrf_plot_causality consumes. Use it for either kind of
% condition contrast:
%   * within-data task states  -- rest_stim vs nback-stimblock (each computed
%     with hrf_causality(..., 'Condition', ...));
%   * between-directory states -- acceptance vs experience (each computed from
%     its own dir set, e.g. {lf_acc,obs_acc} vs {lf_exp,obs_exp}).
%
% Pairing is automatic: if the two results share >=2 subject ids it does a
% PAIRED (within-subject) one-sample t on the per-subject differences (the
% right test when the same people did both conditions); otherwise it falls
% back to an unpaired two-sample (Welch) t.
%
% :Usage:
% ::
%     Rr = hrf_causality({lf,obs}, 'Condition','rest_stim');
%     Rn = hrf_causality({lf,obs}, 'Condition','nback-stimblock');
%     C  = hrf_causality_contrast(Rr, Rn, 'Label1','rest', 'Label2','nback');
%     hrf_plot_causality(C);     % net_group = rest - nback directed flow
%
% :Inputs:
%   **R1, R2:** hrf_causality / hrf_causality_analyze result structs (must
%             carry .net_subj and .subjects; share at least one node).
%
% :Optional Inputs:
%   **'Mode':**    evoked mode to contrast ('remove'/'keep'); default = first
%                  mode present in both.
%   **'Paired':**  [] auto (default), true (force paired on shared subjects),
%                  or false (force unpaired two-sample).
%   **'Label1'/'Label2':** names for the two conditions (for the title).
%
% :Output:
%   **C:** struct with .modes={Mode}, .nodes, .subjects (paired) and
%          .(Mode).net_group (R1-R2) / .t / .p / .p_fdr / .net_subj, plus
%          .paired, .n, .label1, .label2 -- ready for hrf_plot_causality.
%
% See also: hrf_causality, hrf_causality_analyze, hrf_plot_causality.

p = inputParser;
p.addRequired('R1', @isstruct);
p.addRequired('R2', @isstruct);
p.addParameter('Mode', '', @(x) ischar(x) || isstring(x));
p.addParameter('Paired', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.addParameter('Label1', 'cond1', @(x) ischar(x) || isstring(x));
p.addParameter('Label2', 'cond2', @(x) ischar(x) || isstring(x));
p.addParameter('Correction', 'fdr', @(x) ischar(x) || isstring(x));
p.addParameter('GroupNperm', 5000, @(x) isscalar(x) && x >= 100);
p.parse(R1, R2, varargin{:});
opts = p.Results;

mode = local_common_mode(R1, R2, opts.Mode);
[A, B, nodes] = local_align(R1, R2, mode);     % A,B: [N x N x nSubj_i] on shared nodes
N = numel(nodes);

s1 = local_subjects(R1, size(A, 3));
s2 = local_subjects(R2, size(B, 3));
shared = intersect(s1, s2, 'stable');
paired = opts.Paired;
if isempty(paired), paired = numel(shared) >= 2; else, paired = logical(paired); end

if paired
    [~, i1] = ismember(shared, s1);
    [~, i2] = ismember(shared, s2);
    D = A(:, :, i1) - B(:, :, i2);
    G = hrf_group_stats(D, 'Mask', ~eye(N), 'Correction', opts.Correction, 'Nperm', opts.GroupNperm);
    stats = local_from_G(G, nodes); stats.net_subj = D;
    C = local_pack(mode, nodes, stats, true, numel(shared), opts);
    C.subjects = shared(:)';
else
    AB = cat(3, A, B);
    grp = [ones(1, size(A, 3)), 2 * ones(1, size(B, 3))];
    G = hrf_group_stats(AB, 'Mask', ~eye(N), 'Design', 'twosample', 'Group', grp, ...
        'Correction', opts.Correction, 'Nperm', opts.GroupNperm);
    stats = local_from_G(G, nodes); stats.net_subj = [];
    C = local_pack(mode, nodes, stats, false, [size(A, 3) size(B, 3)], opts);
end
fprintf('hrf_causality_contrast: %s - %s | %s | %s test | %s\n', ...
    char(opts.Label1), char(opts.Label2), mode, ...
    local_tern(paired, sprintf('PAIRED n=%d', numel(shared)), ...
        sprintf('unpaired n=%d/%d', size(A, 3), size(B, 3))), ...
    sprintf('%d nodes', N));
end


% =========================================================================
function mode = local_common_mode(R1, R2, want)
m1 = local_modes(R1); m2 = local_modes(R2);
common = intersect(m1, m2, 'stable');
if isempty(common), error('hrf_causality_contrast:NoMode', 'R1 and R2 share no evoked mode.'); end
if ~isempty(char(want))
    if ~ismember(char(want), common)
        error('hrf_causality_contrast:Mode', 'Mode ''%s'' not present in both results.', char(want));
    end
    mode = char(want);
else
    mode = common{1};
end
end

function m = local_modes(R)
if isfield(R, 'modes'), m = cellstr(string(R.modes));
else, m = intersect({'remove', 'keep'}, fieldnames(R)'); end
end

function [A, B, nodes] = local_align(R1, R2, mode)
n1 = cellstr(string(R1.nodes)); n2 = cellstr(string(R2.nodes));
nodes = intersect(n1, n2, 'stable');
if numel(nodes) < 2, error('hrf_causality_contrast:Nodes', 'R1 and R2 share <2 nodes.'); end
A = local_subnet(R1.(mode).net_subj, n1, nodes);
B = local_subnet(R2.(mode).net_subj, n2, nodes);
end

function S = local_subnet(net_subj, nodes_in, nodes_keep)
[~, loc] = ismember(nodes_keep, nodes_in);
S = net_subj(loc, loc, :);
end

function subj = local_subjects(R, n)
if isfield(R, 'subjects') && ~isempty(R.subjects)
    subj = cellstr(string(R.subjects));
else
    subj = arrayfun(@(i) sprintf('_s%d', i), 1:n, 'uni', 0);   % unnamed -> unpaired
end
end

function stats = local_onesample(D, nodes)
[N, ~, n] = size(D);
mu = mean(D, 3, 'omitnan');
sd = std(D, 0, 3, 'omitnan');
se = sd / sqrt(n);
t = mu ./ se; t(se == 0) = 0;
pv = 2 * (1 - local_tcdf(abs(t), n - 1));
pv(logical(eye(N))) = NaN;
stats = struct('net_group', mu, 't', t, 'p', pv, 'p_fdr', local_fdr(pv), 'nodes', {nodes});
end

function stats = local_twosample(A, B, nodes)
N = size(A, 1);
mA = mean(A, 3, 'omitnan'); mB = mean(B, 3, 'omitnan');
vA = var(A, 0, 3, 'omitnan'); vB = var(B, 0, 3, 'omitnan');
nA = sum(isfinite(A), 3); nB = sum(isfinite(B), 3);
se = sqrt(vA ./ max(nA, 1) + vB ./ max(nB, 1));
t = (mA - mB) ./ se; t(se == 0) = 0;
df = (vA ./ max(nA, 1) + vB ./ max(nB, 1)) .^ 2 ./ ...
    ((vA ./ max(nA, 1)) .^ 2 ./ max(nA - 1, 1) + (vB ./ max(nB, 1)) .^ 2 ./ max(nB - 1, 1));
df(~isfinite(df) | df < 1) = 1;
pv = 2 * (1 - local_tcdf(abs(t), df));
pv(logical(eye(N))) = NaN;
stats = struct('net_group', mA - mB, 't', t, 'p', pv, 'p_fdr', local_fdr(pv), 'nodes', {nodes});
end

function stats = local_from_G(G, nodes)
stats = struct('net_group', G.est, 't', G.t, 'p', G.p, 'p_fdr', G.p_corr, ...
    'sig', G.sig, 'nodes', {nodes});
end

function C = local_pack(mode, nodes, stats, paired, n, opts)
C = struct('modes', {{mode}}, 'nodes', string(nodes), 'paired', paired, 'n', n, ...
    'label1', char(opts.Label1), 'label2', char(opts.Label2), ...
    'contrast', sprintf('%s - %s', char(opts.Label1), char(opts.Label2)));
C.(mode) = stats;
end

function pfdr = local_fdr(pv)
pfdr = nan(size(pv));
mask = isfinite(pv); ps = pv(mask); m = numel(ps);
if m == 0, return; end
[sp, ord] = sort(ps(:));
adj = sp .* m ./ (1:m)';
adj = min(1, flipud(cummin(flipud(adj))));
out = nan(m, 1); out(ord) = adj; pfdr(mask) = out;
end

function pcdf = local_tcdf(tval, df)
x = df ./ (df + tval .^ 2);
pcdf = 1 - 0.5 * betainc(x, df / 2, 0.5);
end

function s = local_tern(c, a, b)
if c, s = a; else, s = b; end
end
