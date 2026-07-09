function out = count_map(obj, varargin)
% count_map  Count-map / count-table across subjects for an RSM/RDM.
%
% For a subject-stacked rsm (.dat is k x k x N, one slice per subject) counts,
% in every RSM cell / grouping-block / contrast, HOW MANY subjects meet a
% criterion -- the "consistency" summary you report in a paper (e.g. "in 9/11
% subjects hot patterns were more similar to each other than to warm"). Returns
% both a count-MAP (a matrix you can imagesc) and a count-TABLE (one row per
% cell), the table also carrying the group effect and a group p-value.
%
% Usage
% -----
%   out = R.count_map()                                   % block map, sign criterion
%   out = R.count_map('Granularity','full')               % k x k condition map
%   out = R.count_map('Criterion','threshold','Threshold',0.1)
%   out = R.count_map('Criterion','permutation','Nperm',2000)   % per-subject test
%   out = R.count_map('Granularity','contrasts', 'Contrasts', ...
%             { {'hot','hot'}, {'hot','warm'}, {'hot','imagine'} })
%
% Optional name-value
% -------------------
%   'Granularity' 'blocks' (default; uses .groupings -> G x G map) |
%                 'full' (k x k condition pairs) | 'contrasts' (a list).
%   'Contrasts'   for 'contrasts': cell array of {a,b} pairs (a,b are grouping
%                 names, index vectors, or logical masks; a==b => within-block).
%   'Criterion'   how a subject is counted in a cell:
%                 'sign' (default)      -- oriented value beats 0 (more similar);
%                 'threshold'           -- oriented value beats 'Threshold';
%                 'permutation'         -- the subject's cell is significant vs a
%                                          within-subject label-permutation null.
%   'Threshold'   scalar for 'threshold'. Default 0.
%   'Tail'        'right' (default; more similar) | 'left' | 'both'.
%   'Alpha'       significance for 'permutation'. Default 0.05.
%   'Nperm'       within-subject permutations for 'permutation'. Default 1000.
%   'Transform'   'auto' (default; Fisher-z for correlation-like RSMs) |
%                 'fisherz' | 'none'.
%   'Reduction'   'mean' (default) | 'median' | 'sum' -- how each block reduces.
%   'Correction'  group-level p correction over cells (hrf_group_stats):
%                 'fdr' (default) | 'permutation' | 'none'.
%   'GroupNperm'  permutations for the group correction. Default 5000.
%
% Output struct `out`
% -------------------
%   .count_map   [G x G] / [k x k] matrix of counts (contrasts: [nc x 1]).
%   .prop_map    counts / N.
%   .n           number of subjects.
%   .labels      block / condition / contrast names (rows of the map).
%   .table       one row per cell: name_a, name_b, count, n, proportion,
%                group_mean, group_t, group_p (uncorrected), group_p_corr, sig.
%   .criterion, .granularity, .tail, .is_dissimilarity.
%
% Higher = more similar: for an RDM (.is_dissimilarity true) the value is
% oriented (negated) before the criterion so 'right'/'more similar' always
% means the same thing; the reported group_mean is the raw (dis)similarity.
%
% See also: rsm/cells, rsm/cells_table, rsm/ttest_contrasts, hrf_group_stats.

if builtin('numel', obj) > 1, error('rsm:count_map:nonScalar', 'count_map expects a scalar rsm.'); end

p = inputParser;
p.addParameter('Granularity', 'blocks', @(x) ischar(x) || isstring(x));
p.addParameter('Contrasts', {}, @iscell);
p.addParameter('Criterion', 'sign', @(x) ischar(x) || isstring(x));
p.addParameter('Threshold', 0, @isscalar);
p.addParameter('Tail', 'right', @(x) ischar(x) || isstring(x));
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('Nperm', 1000, @(x) isscalar(x) && x >= 100);
p.addParameter('Transform', 'auto', @(x) ischar(x) || isstring(x));
p.addParameter('Reduction', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('Correction', 'fdr', @(x) ischar(x) || isstring(x));
p.addParameter('GroupNperm', 5000, @(x) isscalar(x) && x >= 100);
p.parse(varargin{:});
o = p.Results;

D = obj.dat;
if ismatrix(D), D = reshape(D, size(D, 1), size(D, 2), 1); end
[k, ~, N] = size(D);
tf = local_resolve_transform(o.Transform, obj.metric);
gran = lower(char(o.Granularity));

% ---- assemble the list of cells to count -------------------------------
[specs, mapsz, rowlab, collab] = local_specs(obj, k, gran, o.Contrasts);
nc = numel(specs);
if nc == 0, error('rsm:count_map:noCells', 'No cells to count (check groupings / Contrasts).'); end

% ---- per-subject value in every cell -----------------------------------
V = nan(nc, N);
for c = 1:nc
    for s = 1:N
        V(c, s) = local_spec_value(D(:, :, s), specs(c), tf, o.Reduction);
    end
end
oriented = V; if obj.is_dissimilarity, oriented = -V; end

% ---- per-subject "hit" -> count ----------------------------------------
switch lower(char(o.Criterion))
    case 'threshold', hit = local_tail(oriented, o.Threshold, o.Tail);
    case 'permutation'
        hit = false(nc, N);
        for c = 1:nc
            for s = 1:N
                hit(c, s) = local_subject_perm(D(:, :, s), specs(c), oriented(c, s), obj.is_dissimilarity, ...
                    tf, o.Reduction, o.Tail, o.Alpha, o.Nperm);
            end
        end
    otherwise, hit = local_tail(oriented, 0, o.Tail);       % 'sign'
end
count = sum(hit, 2, 'omitnan');
prop = count / N;

% ---- group-level effect + p over cells (shared engine) -----------------
G = hrf_group_stats(V, 'Correction', o.Correction, 'Nperm', o.GroupNperm);

% ---- pack map + table --------------------------------------------------
name_a = {specs.name_a}'; name_b = {specs.name_b}';
tbl = table(string(name_a), string(name_b), count, repmat(N, nc, 1), prop, ...
    G.est(:), G.t(:), G.p(:), G.p_corr(:), G.sig(:), ...
    'VariableNames', {'name_a', 'name_b', 'count', 'n', 'proportion', ...
    'group_mean', 'group_t', 'group_p', 'group_p_corr', 'sig'});

out = struct();
out.n = N; out.criterion = lower(char(o.Criterion)); out.granularity = gran;
out.tail = lower(char(o.Tail)); out.is_dissimilarity = obj.is_dissimilarity;
out.table = tbl; out.labels = rowlab; out.collabels = collab;
if isempty(mapsz)
    out.count_map = count; out.prop_map = prop;                  % contrasts -> vector
else
    out.count_map = local_to_matrix(count, specs, mapsz);
    out.prop_map = local_to_matrix(prop, specs, mapsz);
end
end


% =========================================================================
function [specs, mapsz, rowlab, collab] = local_specs(obj, k, gran, contrasts)
specs = struct('name_a', {}, 'name_b', {}, 'terms', {}, 'i', {}, 'j', {});
mapsz = []; rowlab = {}; collab = {};
switch gran
    case 'blocks'
        gn = fieldnames(obj.groupings);
        if isempty(gn), error('rsm:count_map:noGroupings', 'Granularity=blocks needs obj.groupings.'); end
        G = numel(gn); mapsz = [G G]; rowlab = gn; collab = gn;
        for i = 1:G
            for j = i:G
                specs(end + 1) = local_spec(gn{i}, gn{j}, local_term(obj.groupings.(gn{i}), obj.groupings.(gn{j})), i, j); %#ok<AGROW>
            end
        end
    case 'full'
        lab = obj.labels; if isempty(lab), lab = arrayfun(@(x) sprintf('c%d', x), 1:k, 'uni', 0); end
        mapsz = [k k]; rowlab = lab; collab = lab;
        for i = 1:k
            for j = i + 1:k
                specs(end + 1) = local_spec(lab{i}, lab{j}, local_term(i, j), i, j); %#ok<AGROW>
            end
        end
    case 'contrasts'
        if isempty(contrasts), error('rsm:count_map:noContrasts', 'Granularity=contrasts needs ''Contrasts''.'); end
        for c = 1:numel(contrasts)
            e = contrasts{c};
            if isstruct(e) && isfield(e, 'within') && isfield(e, 'between')
                a = local_idx(obj, e.within); b = local_idx(obj, e.between);
                if isfield(e, 'name'), nm = e.name; else, nm = sprintf('%s>btw', local_nm(e.within)); end
                specs(end + 1) = local_spec(nm, local_nm(e.between), [local_term(a, a), local_term(a, b)], c, c); %#ok<AGROW>
            elseif iscell(e) && numel(e) == 2
                ai = local_idx(obj, e{1}); bi = local_idx(obj, e{2});
                specs(end + 1) = local_spec(local_nm(e{1}), local_nm(e{2}), local_term(ai, bi), c, c); %#ok<AGROW>
            else
                g = e; if iscell(g), g = g{1}; end
                ai = local_idx(obj, g);
                specs(end + 1) = local_spec(local_nm(g), local_nm(g), local_term(ai, ai), c, c); %#ok<AGROW>
            end
        end
        rowlab = {specs.name_a}';
    otherwise
        error('rsm:count_map:granularity', 'Granularity must be blocks | full | contrasts.');
end
end


function s = local_spec(na, nb, terms, i, j)
s = struct('name_a', char(string(na)), 'name_b', char(string(nb)), 'terms', terms, 'i', i, 'j', j);
end

function t = local_term(a, b)
t = struct('a', a(:)', 'b', b(:)', 'within', isempty(setxor(a(:)', b(:)')));
end

function v = local_spec_value(M, spec, tf, red)
t = spec.terms;
v = local_cell_value(M, t(1).a, t(1).b, t(1).within, tf, red);
if numel(t) > 1
    v = v - local_cell_value(M, t(2).a, t(2).b, t(2).within, tf, red);
end
end


function idx = local_idx(obj, g)
if ischar(g) || isstring(g)
    idx = obj.groupings.(char(g));
elseif islogical(g)
    idx = find(g);
else
    idx = g;
end
idx = idx(:)';
end

function nm = local_nm(g)
if ischar(g) || isstring(g), nm = char(g); else, nm = 'idx'; end
end


function val = local_cell_value(M, ai, bi, within, tf, red)
if within
    na = numel(ai);
    sub = M(ai, ai);
    x = sub(triu(true(na, na), 1));
else
    sub = M(ai, bi); x = sub(:);
end
x = local_transform(double(x), tf);
switch lower(char(red))
    case 'median', val = median(x, 'omitnan');
    case 'sum',    val = sum(x, 'omitnan');
    otherwise,     val = mean(x, 'omitnan');
end
end


function tf = local_subject_perm(M, spec, obs_oriented, is_diss, transform, red, tail, alpha, nperm)
% Within-subject significance: permute the k condition labels, recompute the
% cell, build a null of the oriented value, one/two-sided p, count if < alpha.
k = size(M, 1);
null = nan(nperm, 1);
for pp = 1:nperm
    q = randperm(k);
    v = local_spec_value(M(q, q), spec, transform, red);
    if is_diss, v = -v; end
    null(pp) = v;
end
switch lower(char(tail))
    case 'left',  pval = (1 + sum(null <= obs_oriented)) / (1 + nperm);
    case 'both',  pval = (1 + sum(abs(null) >= abs(obs_oriented))) / (1 + nperm);
    otherwise,    pval = (1 + sum(null >= obs_oriented)) / (1 + nperm);
end
tf = pval < alpha;
end


function h = local_tail(x, thr, tail)
switch lower(char(tail))
    case 'left', h = x < thr;
    case 'both', h = abs(x) > abs(thr);
    otherwise,   h = x > thr;
end
end


function x = local_transform(x, tf)
switch tf
    case 'fisherz', x = atanh(max(min(x, 1 - 1e-7), -1 + 1e-7));
    otherwise,      % 'none'
end
end


function tf = local_resolve_transform(spec, metric)
spec = lower(char(spec));
if strcmp(spec, 'fisherz'), tf = 'fisherz'; return; end
if strcmp(spec, 'none'), tf = 'none'; return; end
if any(strcmpi(char(metric), {'correlation', 'spearman', 'cosine'})), tf = 'fisherz'; else, tf = 'none'; end
end


function M = local_to_matrix(vals, specs, mapsz)
M = nan(mapsz);
for c = 1:numel(specs)
    M(specs(c).i, specs(c).j) = vals(c);
    M(specs(c).j, specs(c).i) = vals(c);
end
end
