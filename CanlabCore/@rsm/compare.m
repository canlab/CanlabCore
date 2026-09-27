function result = compare(obj, model_rdms, varargin)
% compare  Formal comparison of a data RDM to candidate model RDMs.
%
% Implements the Nili et al. (2014) RSA inference framework: correlates the
% data RDM (per subject) with one or more candidate model RDMs, tests each
% model's relatedness to the data, tests differences between models, and
% estimates the noise ceiling (the range of correlations the true model
% could achieve given the noise in the data).
%
% Wraps the Kriegeskorte rsatoolbox `rsa.compareRefRDM2candRDMs` when it is
% on the path; otherwise uses an equivalent built-in implementation.
%
% Usage
% -----
%   R = compute_rsm(dat, 'group_by',{'condition','bodysite'}, 'subject_var','sub');
%   models = rsm.from_categorical(R.metadata_table, {'condition','bodysite'});
%   result = R.compare(models, 'correlation_type','kendall_taua');
%
%   result.r_mean            % group-mean correlation per model
%   result.relatedness_p     % is each model related to the data?
%   result.noise_ceiling     % [lower upper] bounds
%
% Inputs
% ------
%   obj         a data rsm. If 3D (per-subject slices), random-effects
%               inference is run across subjects. If 2D, only the
%               correlations + a permutation relatedness test are available.
%   model_rdms  candidate model RDMs. Accepts:
%                 - an rsm or array of rsm (each k x k)
%                 - a numeric k x k matrix
%                 - a metadata column name / cellstr (builds same-vs-different
%                   model RDMs from obj.metadata_table)
%
% Optional name-value
% -------------------
%   'correlation_type'  'kendall_taua' (default, recommended for categorical
%                       models per Nili 2014) | 'spearman' | 'pearson'
%   'relatedness_test'  'signrank' (default, subject RFX Wilcoxon signed-rank)
%                       | 'ttest' | 'none'
%   'differences_test'  'signrank' (default) | 'ttest' | 'none'
%   'multiple_testing'  'fdr' (default) | 'bonferroni' | 'none'
%   'noise_ceiling'     logical (default true). Estimate upper/lower bounds.
%   'tail'              'right' (default) | 'both'   (relatedness)
%   'engine'            'auto' (default) | 'rsatoolbox' | 'builtin'
%   'verbose'           logical (default true)
%
% Output
% ------
%   result struct:
%     .candidate_names    {1 x nModels}
%     .r                  [nSubjects x nModels] per-subject correlations
%     .r_mean             [1 x nModels] group-mean correlation
%     .r_se               [1 x nModels] standard error across subjects
%     .relatedness_p      [1 x nModels]
%     .relatedness_p_corr [1 x nModels]
%     .relatedness_sig    [1 x nModels] logical
%     .differences_p      [nModels x nModels] pairwise model comparison
%     .noise_ceiling      [lower upper]
%     .correlation_type, .engine, .history

if builtin('numel', obj) > 1
    error('rsm:compare:nonScalar', 'compare expects a scalar rsm.');
end

p = inputParser;
p.addParameter('correlation_type', 'kendall_taua', @(x) (ischar(x)||isstring(x)) && ismember(lower(char(x)), {'kendall_taua','kendall','spearman','pearson'}));
p.addParameter('relatedness_test', 'signrank', @(x) (ischar(x)||isstring(x)) && ismember(lower(char(x)), {'signrank','ttest','none'}));
p.addParameter('differences_test', 'signrank', @(x) (ischar(x)||isstring(x)) && ismember(lower(char(x)), {'signrank','ttest','none'}));
p.addParameter('multiple_testing', 'fdr', @(x) (ischar(x)||isstring(x)) && ismember(lower(char(x)), {'fdr','bonferroni','none'}));
p.addParameter('noise_ceiling',    true,  @(x) islogical(x) || isnumeric(x));
p.addParameter('tail',             'right', @(x) (ischar(x)||isstring(x)) && ismember(lower(char(x)), {'right','both'}));
p.addParameter('engine',           'auto', @(x) (ischar(x)||isstring(x)) && ismember(lower(char(x)), {'auto','rsatoolbox','builtin'}));
p.addParameter('verbose',          true,  @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

k = size(obj.dat, 1);

% =========================================================================
% Build data RDM stack (dissimilarity space, per subject)
% =========================================================================
% Convert similarity RSM -> dissimilarity RDM if needed
data_stack = obj.dat;
if ~obj.is_dissimilarity
    data_stack = 1 - data_stack;
end
n_subj = size(data_stack, 3);

% =========================================================================
% Normalize model RDMs to dissimilarity space
% =========================================================================
[models, model_names] = normalize_compare_models(model_rdms, k, obj);
n_models = numel(models);

% =========================================================================
% Engine selection
% =========================================================================
caps = probe_rsatoolbox(); %#ok<NASGU>
% The built-in engine is well-validated, side-effect-free (no figure/file
% saving), and exposes per-subject correlations, so 'auto' prefers it. Use
% 'engine','rsatoolbox' explicitly to route through rsa.compareRefRDM2candRDMs.
use_rsatoolbox = strcmpi(opt.engine, 'rsatoolbox');

result = struct();
result.candidate_names = model_names;
result.correlation_type = lower(char(opt.correlation_type));
result.history = {sprintf('%s: rsm.compare (%d models, %d subjects)', datestr(now), n_models, n_subj)};

if use_rsatoolbox
    try
        result = compare_via_rsatoolbox(data_stack, models, model_names, opt, result);
        result.engine = 'rsatoolbox';
        if opt.verbose, fprintf('rsm.compare: used rsa.compareRefRDM2candRDMs.\n'); end
        return
    catch ME
        if opt.verbose
            warning('rsm:compare:rsatoolboxFailed', ...
                'rsatoolbox engine failed (%s); falling back to builtin.', ME.message);
        end
    end
end

% =========================================================================
% Built-in implementation
% =========================================================================
result = compare_builtin(data_stack, models, model_names, opt, result);
result.engine = 'builtin';
if opt.verbose
    fprintf('rsm.compare: %d models, %d subjects, %s correlation.\n', ...
        n_models, n_subj, result.correlation_type);
end

end


% =========================================================================
function result = compare_builtin(data_stack, models, model_names, opt, result)

[k, ~, n_subj] = size(data_stack);
n_models = numel(models);
tril_mask = tril(true(k), -1);

% Vectorize model RDMs
model_vecs = cellfun(@(M) M(tril_mask), models, 'UniformOutput', false);

% Per-subject correlations to each model
r = nan(n_subj, n_models);
for s = 1:n_subj
    dvec = data_stack(:, :, s);
    dvec = dvec(tril_mask);
    for m = 1:n_models
        r(s, m) = rdm_corr(dvec, model_vecs{m}, opt.correlation_type);
    end
end

result.r      = r;
result.r_mean = mean(r, 1, 'omitnan');
result.r_se   = std(r, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(r), 1));

% Relatedness test (each model vs 0 across subjects)
rel_p = nan(1, n_models);
if n_subj >= 3 && ~strcmpi(opt.relatedness_test, 'none')
    for m = 1:n_models
        rel_p(m) = group_test(r(:, m), opt.relatedness_test, opt.tail);
    end
elseif n_subj < 3
    if opt.verbose
        warning('rsm:compare:fewSubjects', ...
            'Random-effects inference needs >= 3 subjects (have %d); relatedness_p = NaN.', n_subj);
    end
end
result.relatedness_p = rel_p;

% Multiple-comparison correction across models
[rel_p_corr, rel_sig] = correct_p(rel_p, opt.multiple_testing);
result.relatedness_p_corr = rel_p_corr;
result.relatedness_sig    = rel_sig;

% Pairwise differences test
diff_p = nan(n_models, n_models);
if n_subj >= 3 && n_models >= 2 && ~strcmpi(opt.differences_test, 'none')
    for a = 1:n_models
        for b = (a+1):n_models
            d = r(:, a) - r(:, b);
            pdiff = group_test(d, opt.differences_test, 'both');
            diff_p(a, b) = pdiff;
            diff_p(b, a) = pdiff;
        end
    end
end
result.differences_p = diff_p;

% Noise ceiling (Nili et al. 2014)
if logical(opt.noise_ceiling) && n_subj >= 2
    result.noise_ceiling = estimate_noise_ceiling(data_stack, opt.correlation_type);
else
    result.noise_ceiling = [NaN NaN];
end

end


% =========================================================================
function nc = estimate_noise_ceiling(data_stack, corr_type)
% Upper bound: mean correlation of each subject's RDM with the group-mean
% RDM (incl. that subject). Lower bound: leave-one-out (mean of others).
[k, ~, n_subj] = size(data_stack);
tril_mask = tril(true(k), -1);

vecs = zeros(nnz(tril_mask), n_subj);
for s = 1:n_subj
    slice = data_stack(:, :, s);
    vecs(:, s) = slice(tril_mask);
end
group_mean = mean(vecs, 2, 'omitnan');

upper_rs = nan(n_subj, 1);
lower_rs = nan(n_subj, 1);
for s = 1:n_subj
    upper_rs(s) = rdm_corr(vecs(:, s), group_mean, corr_type);
    others = mean(vecs(:, setdiff(1:n_subj, s)), 2, 'omitnan');
    lower_rs(s) = rdm_corr(vecs(:, s), others, corr_type);
end
nc = [mean(lower_rs, 'omitnan'), mean(upper_rs, 'omitnan')];
end


% =========================================================================
function r = rdm_corr(a, b, corr_type)
a = a(:); b = b(:);
mask = isfinite(a) & isfinite(b);
if nnz(mask) < 3, r = NaN; return; end
a = a(mask); b = b(mask);
switch lower(char(corr_type))
    case 'kendall_taua', r = kendall_taua(a, b);
    case 'kendall',      r = corr(a, b, 'Type', 'Kendall');
    case 'pearson',      r = corr(a, b, 'Type', 'Pearson');
    otherwise,           r = corr(a, b, 'Type', 'Spearman');
end
end


% =========================================================================
function tau = kendall_taua(x, y)
% Kendall's tau-a: (n_concordant - n_discordant) / (n*(n-1)/2).
% Unlike tau-b, tau-a does not correct for ties -- the RSA-recommended
% statistic when model RDMs predict many tied ranks (Nili et al. 2014).
n = numel(x);
if n < 2, tau = NaN; return; end
% Pairwise sign products summed via vectorized upper triangle
[I, J] = find(triu(true(n), 1));
dx = sign(x(I) - x(J));
dy = sign(y(I) - y(J));
tau = sum(dx .* dy) / (n*(n-1)/2);
end


% =========================================================================
function pval = group_test(vals, test_type, tail)
vals = vals(~isnan(vals));
if numel(vals) < 3, pval = NaN; return; end
switch lower(char(test_type))
    case 'signrank'
        if strcmpi(tail, 'right')
            pval = signrank(vals, 0, 'tail', 'right');
        else
            pval = signrank(vals, 0);
        end
    case 'ttest'
        if strcmpi(tail, 'right')
            [~, pval] = ttest(vals, 0, 'Tail', 'right');
        else
            [~, pval] = ttest(vals, 0);
        end
    otherwise
        pval = NaN;
end
end


% =========================================================================
function [adj, sig] = correct_p(ps, method)
valid = ~isnan(ps);
adj = ps; sig = false(size(ps));
if ~any(valid), return; end
pv = ps(valid);
n = numel(pv);
switch lower(char(method))
    case 'fdr'
        [~, order] = sort(pv); ranks = zeros(n,1); ranks(order) = 1:n;
        adj_v = min(pv(:) .* n ./ ranks(:), 1);
        sig_v = adj_v < 0.05;
    case 'bonferroni'
        adj_v = min(pv(:) * n, 1);
        sig_v = adj_v < 0.05;
    otherwise
        adj_v = pv(:);
        sig_v = pv(:) < 0.05;
end
adj(valid) = adj_v;
sig(valid) = sig_v;
end


% =========================================================================
function result = compare_via_rsatoolbox(data_stack, models, model_names, opt, result)
% Wrap rsa.compareRefRDM2candRDMs. Builds the ref RDM (per-subject stack)
% and candidate RDM structs in the rsatoolbox format.
[k, ~, n_subj] = size(data_stack);

% Reference RDMs: rsatoolbox accepts a [k x k x nSubj] stack
refRDMs = data_stack;

% Candidate RDMs: struct array with .RDM and .name
candRDMs = struct('RDM', {}, 'name', {}, 'color', {});
for m = 1:numel(models)
    candRDMs(m).RDM = models{m};
    candRDMs(m).name = model_names{m};
    candRDMs(m).color = [0 0 0];
end

% Map our correlation-type names to rsatoolbox's
switch lower(char(opt.correlation_type))
    case 'kendall_taua', corrtype = 'Kendall_taua';
    case 'spearman',     corrtype = 'Spearman';
    case 'pearson',      corrtype = 'Pearson';
    otherwise,           corrtype = 'Kendall_taua';
end

userOptions = struct();
% Required project-setup fields for rsatoolbox's file/figure machinery
userOptions.analysisName = 'canlab_rsm_compare';
userOptions.projectName  = 'canlab_rsm_compare';
userOptions.rootPath     = tempdir;
userOptions.RDMcorrelationType = corrtype;
userOptions.RDMrelatednessTest = 'subjectRFXsignedRank';
userOptions.RDMrelatednessThreshold = 0.05;
userOptions.RDMrelatednessMultipleTesting = upper(char(opt.multiple_testing));
userOptions.RDMcandRelatednessTest = 'subjectRFXsignedRank';
userOptions.plotpValues = '=';
userOptions.saveFigurePDF = false;
userOptions.saveFigurePNG = false;
userOptions.saveFigureFig = false;
userOptions.displayFigures = false;

stats_p_r = rsa.compareRefRDM2candRDMs(refRDMs, candRDMs, userOptions);

% Extract what we can into our standard result struct
result.r_mean = stats_p_r.candRelatedness_r(:)';
if isfield(stats_p_r, 'candRelatedness_p')
    result.relatedness_p = stats_p_r.candRelatedness_p(:)';
end
if isfield(stats_p_r, 'ceiling')
    result.noise_ceiling = stats_p_r.ceiling(:)';
end
result.rsatoolbox_raw = stats_p_r;
result.r = [];   % per-subject not always exposed by the toolbox call
[result.relatedness_p_corr, result.relatedness_sig] = ...
    correct_p(result.relatedness_p, opt.multiple_testing);
end


% =========================================================================
function [models, names] = normalize_compare_models(model_rdms, k, obj)
% Coerce candidate models to a cell of k x k dissimilarity matrices + names.
models = {}; names = {};

% Column-name form: build from obj.metadata_table
if ischar(model_rdms) || iscellstr(model_rdms) || isstring(model_rdms) %#ok<ISCLSTR>
    cols = cellstr(model_rdms);
    if isempty(obj.metadata_table)
        error('rsm:compare:noMetadata', ...
            'Column-name models require obj.metadata_table.');
    end
    for i = 1:numel(cols)
        Mrsm = rsm.from_categorical(obj.metadata_table.(cols{i}));
        models{end+1} = Mrsm.dat; %#ok<AGROW>
        names{end+1}  = cols{i};  %#ok<AGROW>
    end
elseif isa(model_rdms, 'rsm')
    for i = 1:numel(model_rdms)
        M = mean(model_rdms(i).dat, 3, 'omitnan');
        if ~model_rdms(i).is_dissimilarity, M = 1 - M; end
        models{end+1} = M; %#ok<AGROW>
        nm = regexprep(model_rdms(i).source, '^design:', '');
        nm = regexprep(nm, '^parcel:', '');
        if isempty(nm), nm = sprintf('model%d', i); end
        names{end+1} = nm; %#ok<AGROW>
    end
elseif isnumeric(model_rdms)
    models{1} = model_rdms; names{1} = 'model1';
elseif isstruct(model_rdms)
    for i = 1:numel(model_rdms)
        models{end+1} = model_rdms(i).RDM; %#ok<AGROW>
        if isfield(model_rdms(i), 'name') && ~isempty(model_rdms(i).name)
            names{end+1} = model_rdms(i).name; %#ok<AGROW>
        else
            names{end+1} = sprintf('model%d', i); %#ok<AGROW>
        end
    end
else
    error('rsm:compare:badModel', 'model_rdms must be a column name, rsm, numeric, or struct.');
end

for i = 1:numel(models)
    if ~isequal(size(models{i}), [k k])
        error('rsm:compare:modelSize', ...
            'Model %d is %dx%d but data RDM is %dx%d.', i, size(models{i},1), size(models{i},2), k, k);
    end
end
end
