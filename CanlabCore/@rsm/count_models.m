function out = count_models(obj, model_rdms, varargin)
% count_models  Count how many subjects each candidate model RDM explains.
%
% For a subject-stacked rsm (.dat is k x k x N) correlates every subject's
% RSM/RDM with one or more candidate MODEL RDMs and reports, across subjects,
% how many each model (a) WINS (is the best-fitting model for that subject) and
% (b) is significantly related to the data -- the model-selection analogue of
% count_map, and the "the categorical model was the best fit in 9/11 subjects"
% summary you put in a paper. Complements @rsm/compare (Nili 2014 group RFX):
% compare answers "is this model related / better on average", count_models
% answers "in how many individual subjects".
%
% Usage
% -----
%   R = compute_rsm(dat, 'group_by',{'condition','bodysite'}, 'subject_var','sub');
%   out = R.count_models({'condition','bodysite'});            % same-vs-different models
%   out = R.count_models(model_rsm_array, 'Criterion','permutation');
%   out.table                                                  % model | wins | n_sig | group stats
%
% Inputs
% ------
%   obj         a scalar per-subject-stacked rsm (k x k x N).
%   model_rdms  candidate models, as in @rsm/compare:
%                 - an rsm or array of rsm (each k x k)
%                 - a numeric k x k matrix (or cell of them)
%                 - a metadata column name / cellstr -> same-vs-different model
%                   RDMs built from obj.metadata_table (rsm.from_categorical).
%
% Optional name-value
% -------------------
%   'CorrelationType' 'spearman' (default) | 'kendall' | 'pearson'.
%   'Criterion'   how a subject counts as fit by a model:
%                 'sign' (default; oriented r>0) | 'threshold' | 'permutation'
%                 (per-subject RDM-label permutation null).
%   'Threshold'   scalar for 'threshold'. Default 0.
%   'Tail'        'right' (default) | 'left' | 'both'.
%   'Alpha'       significance for 'permutation'. Default 0.05.
%   'Nperm'       per-subject permutations. Default 1000.
%   'Winner'      'exclusive' (default; argmax over models, 1 win/subject) |
%                 'significant' (a subject can be "won" by every model it is
%                 significantly related to) | 'none'.
%   'Correction'  group-level p over models: 'fdr' (default) | 'permutation' | 'none'.
%   'GroupNperm'  permutations for the group correction. Default 5000.
%   'doplot'      logical (default false) -- bar chart of wins / n_sig.
%
% Output struct `out`
% -------------------
%   .r            [nModels x N] per-subject oriented correlations (similarity
%                 space; higher = better fit).
%   .model_names  {nModels x 1}.
%   .wins         [nModels x 1] # subjects for which the model is best.
%   .n_sig        [nModels x 1] # subjects the model significantly fits.
%   .n            number of subjects.
%   .table        model, wins, prop_wins, n_sig, prop_sig, group_mean_r,
%                 group_t, group_p, group_p_corr, sig.
%   .figure       handle if doplot.
%
% See also: rsm/compare, rsm/count_map, rsm/from_categorical.

if builtin('numel', obj) > 1, error('rsm:count_models:nonScalar', 'count_models expects a scalar rsm.'); end

p = inputParser;
p.addRequired('model_rdms');
p.addParameter('CorrelationType', 'spearman', @(x) ischar(x) || isstring(x));
p.addParameter('Criterion', 'sign', @(x) ischar(x) || isstring(x));
p.addParameter('Threshold', 0, @isscalar);
p.addParameter('Tail', 'right', @(x) ischar(x) || isstring(x));
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.addParameter('Nperm', 1000, @(x) isscalar(x) && x >= 100);
p.addParameter('Winner', 'exclusive', @(x) ischar(x) || isstring(x));
p.addParameter('Correction', 'fdr', @(x) ischar(x) || isstring(x));
p.addParameter('GroupNperm', 5000, @(x) isscalar(x) && x >= 100);
p.addParameter('doplot', false, @(x) islogical(x) || isnumeric(x));
p.parse(model_rdms, varargin{:});
o = p.Results;
ctype = lower(char(o.CorrelationType));

D = obj.dat;
if ismatrix(D), D = reshape(D, size(D, 1), size(D, 2), 1); end
[k, ~, N] = size(D);

% ---- build model vectors (oriented to similarity space) ----------------
[mnames, mvec] = local_models(obj, o.model_rdms, k);
nm = numel(mnames);
if nm == 0, error('rsm:count_models:noModels', 'No candidate models resolved.'); end
mask = triu(true(k, k), 1);

% ---- per-subject correlation with every model --------------------------
R = nan(nm, N);
data_sign = 1; if obj.is_dissimilarity, data_sign = -1; end   % -> similarity space
for s = 1:N
    Ms = D(:, :, s); dv = data_sign * Ms(mask);
    for m = 1:nm
        R(m, s) = local_corr(dv, mvec{m}, ctype);
    end
end

% ---- per-subject "hit" per model ---------------------------------------
switch lower(char(o.Criterion))
    case 'threshold', hit = local_tail(R, o.Threshold, o.Tail);
    case 'permutation'
        hit = false(nm, N);
        for s = 1:N
            Ms = D(:, :, s);
            for m = 1:nm
                hit(m, s) = local_model_perm(Ms, mvec{m}, mask, data_sign, R(m, s), ...
                    ctype, o.Tail, o.Alpha, o.Nperm);
            end
        end
    otherwise, hit = local_tail(R, 0, o.Tail);                 % 'sign'
end
n_sig = sum(hit, 2);

% ---- winners -----------------------------------------------------------
wins = zeros(nm, 1);
switch lower(char(o.Winner))
    case 'none'
        % no winner tally
    case 'significant'
        wins = n_sig;                                          % each sig model "wins" that subject
    otherwise                                                  % 'exclusive' argmax
        for s = 1:N
            [best, mi] = max(R(:, s));
            if isfinite(best), wins(mi) = wins(mi) + 1; end
        end
end

% ---- group effect over models (Fisher-z for the test) ------------------
G = rsm_group_stats(atanh(max(min(R, 1 - 1e-7), -1 + 1e-7)), o.Correction, o.GroupNperm);
group_mean_r = mean(R, 2, 'omitnan');

tbl = table(string(mnames(:)), wins, wins / N, n_sig, n_sig / N, group_mean_r, ...
    G.t(:), G.p(:), G.p_corr(:), G.sig(:), 'VariableNames', ...
    {'model', 'wins', 'prop_wins', 'n_sig', 'prop_sig', 'group_mean_r', ...
    'group_t', 'group_p', 'group_p_corr', 'sig'});

out = struct('r', R, 'model_names', {mnames(:)}, 'wins', wins, 'n_sig', n_sig, ...
    'n', N, 'criterion', lower(char(o.Criterion)), 'correlation_type', ctype, ...
    'winner', lower(char(o.Winner)), 'table', tbl);
if logical(o.doplot), out.figure = local_plot_models(out); end
end


% =========================================================================
function [names, vecs] = local_models(obj, models, k)
% Return {names} and {upper-tri vectors in SIMILARITY space} for each model.
names = {}; vecs = {};
if ischar(models) || isstring(models) || (iscell(models) && all(cellfun(@(x) ischar(x) || isstring(x), models)))
    % metadata column name(s) -> same-vs-different model RDM(s)
    cols = cellstr(string(models));
    for c = 1:numel(cols)
        M = local_categorical_rdm(obj, cols{c}, k);
        names{end + 1} = cols{c}; vecs{end + 1} = local_orient(M, true, k); %#ok<AGROW>
    end
    return
end
if ~iscell(models), models = {models}; end
for m = 1:numel(models)
    e = models{m};
    if isa(e, 'rsm')
        for q = 1:builtin('numel', e)
            eq = e(q);
            names{end + 1} = local_model_name(eq, numel(names) + 1); %#ok<AGROW>
            vecs{end + 1} = local_orient(eq.dat(:, :, 1), eq.is_dissimilarity, k); %#ok<AGROW>
        end
    elseif isnumeric(e)
        names{end + 1} = sprintf('model_%d', numel(names) + 1); %#ok<AGROW>
        vecs{end + 1} = local_orient(e, true, k);              % assume RDM %#ok<AGROW>
    end
end
end


function M = local_categorical_rdm(obj, col, k)
% Same-vs-different (0/1) RDM from a metadata column, k x k.
if isempty(obj.metadata_table) || ~ismember(col, obj.metadata_table.Properties.VariableNames)
    error('rsm:count_models:noMeta', 'Model column "%s" not found in obj.metadata_table.', col);
end
v = obj.metadata_table.(col);
if height(obj.metadata_table) ~= k
    error('rsm:count_models:metaSize', 'metadata_table has %d rows but RSM is %d-D.', height(obj.metadata_table), k);
end
key = local_key(v);
M = double(bsxfun(@ne, key, key'));                            % 1 = different (dissimilar)
end


function vec = local_orient(M, is_diss, k)
% Upper-tri vector oriented to similarity space (higher = more similar).
if size(M, 1) ~= k, error('rsm:count_models:modelSize', 'Model is %dx%d, expected %dx%d.', size(M, 1), size(M, 2), k, k); end
mask = triu(true(k, k), 1);
vec = M(mask);
if is_diss, vec = -vec; end
end


function tf = local_model_perm(Ms, mvec, mask, data_sign, obs, ctype, tail, alpha, nperm)
k = size(Ms, 1);
null = nan(nperm, 1);
for pp = 1:nperm
    q = randperm(k);
    Mp = Ms(q, q); dv = data_sign * Mp(mask);
    null(pp) = local_corr(dv, mvec, ctype);
end
switch lower(char(tail))
    case 'left',  pval = (1 + sum(null <= obs)) / (1 + nperm);
    case 'both',  pval = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
    otherwise,    pval = (1 + sum(null >= obs)) / (1 + nperm);
end
tf = pval < alpha;
end


function rho = local_corr(x, y, ctype)
x = x(:); y = y(:); good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);
if numel(x) < 3 || std(x) == 0 || std(y) == 0, rho = NaN; return; end
switch ctype
    case 'kendall', rho = corr(x, y, 'type', 'Kendall');
    case 'pearson', rho = corr(x, y, 'type', 'Pearson');
    otherwise,      rho = corr(x, y, 'type', 'Spearman');
end
end


function h = local_tail(x, thr, tail)
switch lower(char(tail))
    case 'left', h = x < thr;
    case 'both', h = abs(x) > abs(thr);
    otherwise,   h = x > thr;
end
end


function key = local_key(v)
% Integer key per unique value of a metadata column (numeric/cell/categorical).
if iscell(v) || isstring(v) || ischar(v) || iscategorical(v)
    [~, ~, key] = unique(cellstr(string(v)), 'stable');
else
    [~, ~, key] = unique(v, 'stable');
end
key = key(:);
end


function nm = local_model_name(r, idx)
if ~isempty(r.metric) && ~strcmpi(r.metric, 'none'), nm = r.metric; else, nm = sprintf('model_%d', idx); end
end


function h = local_plot_models(out)
h = figure('Color', 'w', 'Name', 'rsm count_models');
Y = [out.wins, out.n_sig];
bar(Y); ylim([0 out.n]); ylabel('# subjects'); box off;
set(gca, 'XTick', 1:numel(out.model_names), 'XTickLabel', out.model_names, 'XTickLabelRotation', 30);
legend({'wins', 'significant'}, 'Location', 'best'); legend boxoff;
title(sprintf('model fit across %d subjects (%s)', out.n, out.criterion), 'Interpreter', 'none');
end
