function test_rsa_tools()
% test_rsa_tools  Unit test suite for the CanlabCore RSA tools.
%
% Runs a battery of assertions on synthetic data with planted ground truth.
% Each test prints PASS/FAIL. Throws at the end if any test failed.
%
% Usage:
%   test_rsa_tools
%
% Covers: rsm construction, compute_rsm metrics, cells/contrasts,
% reliability, drift, LME, parcelwise, searchlight, compare, and soft-dep
% fallbacks.

fprintf('\n======================================\n');
fprintf('  CanlabCore RSA tools test suite\n');
fprintf('======================================\n\n');

state = struct('n_pass', 0, 'n_fail', 0, 'failures', {{}});

state = run_test(@test_construction,        'rsm construction',        state);
state = run_test(@test_metrics,             'compute_rsm metrics',     state);
state = run_test(@test_crossnobis,          'crossnobis distance',     state);
state = run_test(@test_crossnobis_occurrence, 'crossnobis occurrence folds', state);
state = run_test(@test_cvcorr,              'cross-session cvcorr',    state);
state = run_test(@test_recode_reference,    'rsa_recode_reference',    state);
state = run_test(@test_cells_contrasts,     'cells + ttest_contrasts', state);
state = run_test(@test_reliability,         'reliability ICC',         state);
state = run_test(@test_drift,               'drift',                   state);
state = run_test(@test_lme,                 'rsa_lme',                 state);
state = run_test(@test_compare_models,      'rsa_compare_models',      state);
state = run_test(@test_compare,             'rsm.compare',             state);
state = run_test(@test_model_rdms,          'model RDM constructors',  state);

fprintf('\n--------------------------------------\n');
fprintf('  %d passed, %d failed\n', state.n_pass, state.n_fail);
fprintf('--------------------------------------\n');
if state.n_fail > 0
    error('test_rsa_tools:failures', 'Failed tests: %s', strjoin(state.failures, ', '));
end
end


% =========================================================================
function state = run_test(fn, name, state)
try
    fn();
    fprintf('  [PASS] %s\n', name);
    state.n_pass = state.n_pass + 1;
catch ME
    fprintf('  [FAIL] %s -- %s\n', name, ME.message);
    state.n_fail = state.n_fail + 1;
    state.failures{end+1} = name;
end
end


% =========================================================================
function dat = synth(varargin)
% Small planted dataset (k=12, structure = condition + bodysite).
p = inputParser;
p.addParameter('n_sub', 6); p.addParameter('n_ses', 2); p.addParameter('seed', 1);
p.parse(varargin{:}); o = p.Results;
rng(o.seed);
n_vox = 120;
conds = {'a','b','c'}; bs = {'w','x','y','z'};
cs = randn(3, n_vox); bsig = randn(4, n_vox);
P = zeros(12, n_vox); r = 0;
for c = 1:3, for b = 1:4, r=r+1; P(r,:) = 0.8*cs(c,:) + 0.4*bsig(b,:); end, end
X = []; sv={}; sev=[]; cv={}; bv={}; i=0;
for s = 1:o.n_sub, for se = 1:o.n_ses, for c = 1:3, for b = 1:4
    i=i+1; X(:,i) = P((c-1)*4+b,:)' + 0.4*randn(n_vox,1);
    sv{i,1}=sprintf('s%02d',s); sev(i,1)=se; cv{i,1}=conds{c}; bv{i,1}=bs{b};
end,end,end,end
dat = fmri_data; dat.dat = X;
dat.metadata_table = table(sv, sev, cv, bv, 'VariableNames', {'sub','sesno','condition','bodysite'});
end


% =========================================================================
function test_construction()
dat = synth();
R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'metric', 'spearman', 'verbose', false);
assert(isequal(size(R), [12 12 6]), 'expected 12x12x6');
assert(~isempty(fieldnames(R.groupings)), 'groupings should auto-attach');
assert(isfield(R.groupings, 'a'), 'condition grouping missing');
% Block structure: within-condition > between-condition
m = mean(R.dat, 3);
within = mean([m(1,2) m(1,3) m(2,3)]);   % a-block off diagonal (3 of the 4 bs)
between = mean([m(1,5) m(1,9)]);          % a vs b/c
assert(within > between, 'within-condition should exceed between');
end


% =========================================================================
function test_metrics()
dat = synth();
for metric = {'correlation','spearman','cosine','euclidean'}
    R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
        'metric', metric{1}, 'verbose', false);
    assert(isequal(size(R), [12 12 6]), sprintf('%s size', metric{1}));
    is_dis = ismember(metric{1}, {'euclidean'});
    assert(R.is_dissimilarity == is_dis, sprintf('%s is_dissimilarity', metric{1}));
end
end


% =========================================================================
function test_crossnobis()
dat = synth();
R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'metric', 'crossnobis', 'fold_var', 'sesno', 'verbose', false);
assert(R.is_dissimilarity, 'crossnobis is a dissimilarity');
m = mean(R.dat, 3);
% within-condition distance < between-condition distance
within = mean([m(1,2) m(1,3) m(2,3)]);
between = mean([m(1,5) m(1,9)]);
assert(within < between, 'crossnobis within < between');
end


% =========================================================================
function test_crossnobis_occurrence()
% Folds defined by occurrence rank should match an explicit fold column.
dat = synth();
mt = dat.metadata_table;
% Build an explicit occurrence-rank fold column over (sub, condition, bodysite)
key = strcat(string(mt.sub), '|', string(mt.condition), '|', string(mt.bodysite));
[~,~,cid] = unique(key); mt.fold = zeros(height(mt),1);
for c = unique(cid)', idx = find(cid==c); mt.fold(idx) = 1:numel(idx); end
dat.metadata_table = mt;
R_manual = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var','sub', ...
    'metric','crossnobis', 'fold_var','fold', 'verbose', false);
R_auto = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var','sub', ...
    'metric','crossnobis', 'fold_var','occurrence', 'verbose', false);
d = max(abs(R_manual.dat(:) - R_auto.dat(:)));
assert(d < 1e-9, sprintf('occurrence vs manual fold diff = %g', d));
end


% =========================================================================
function test_cvcorr()
% Cross-validated (cross-session) correlation: a SIMILARITY whose diagonal is
% the cross-fold reliability (not 1), and which requires a fold_var.
dat = synth();
R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var','sub', ...
    'metric','cvcorr', 'fold_var','sesno', 'verbose', false);
assert(~R.is_dissimilarity, 'cvcorr is a similarity, not a dissimilarity');
m = mean(R.dat, 3, 'omitnan');
dg = diag(m);
assert(~all(abs(dg(~isnan(dg)) - 1) < 1e-6), 'cvcorr diagonal must be cross-fold reliability, not all 1');
% fold_var is required
threw = false;
try
    compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var','sub', ...
        'metric','cvcorr', 'verbose', false);
catch
    threw = true;
end
assert(threw, 'cvcorr without fold_var should error');

% leave-one-fold-out scheme: runs, symmetric, similar structure to allpairs
Rlo = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var','sub', ...
    'metric','cvcorr', 'fold_var','sesno', 'cv_scheme','loo', 'verbose', false);
mlo = mean(Rlo.dat, 3, 'omitnan');
assert(max(abs(mlo - mlo'), [], 'all', 'omitnan') < 1e-9, 'cvcorr loo must be symmetric');
assert(corr(m(~isnan(m(:))), mlo(~isnan(mlo(:))), 'rows','pairwise') > 0.5, ...
    'cvcorr loo and allpairs should be broadly consistent');
end


% =========================================================================
function test_recode_reference()
v = {'Left Face','Right Arm','Left Face','Chest','Abdomen'};
out = rsa_recode_reference(v, 'Left Face', 'other_label', 'Other');
assert(isequal(out(:)', {'Left Face','Other','Left Face','Other','Other'}), 'recode mismatch');
% numeric + multi-reference
out2 = rsa_recode_reference([0 1 2 0 1], {'0','1'});
assert(isequal(out2(:)', {'0','1','Other','0','1'}), 'numeric recode mismatch');
end


% =========================================================================
function test_cells_contrasts()
dat = synth();
R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'metric', 'spearman', 'verbose', false);
within = R.cells('a', 'a');
between = R.cells('a', 'b');
assert(isequal(size(within), [6 1]), 'cells returns per-subject column');
assert(mean(within) > mean(between), 'within > between cells');
% ttest_contrasts
T = R.ttest_contrasts({'wa','a',[]; 'avsb','a','b'}, 'tail', 'right', 'correction', 'fdr');
assert(height(T) == 2, 'two contrasts');
assert(all(ismember({'Contrast','t','P','FDR_P','sig'}, T.Properties.VariableNames)), 'table cols');
% correction='none' must not crash (regression test for duplicate P column)
T2 = R.ttest_contrasts({'wa','a',[]}, 'correction', 'none');
assert(height(T2) == 1, 'none correction works');
end


% =========================================================================
function test_reliability()
dat = synth('n_ses', 4);
R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'level', 'session', ...
    'subject_var', 'sub', 'session_var', 'sesno', 'metric', 'spearman', 'verbose', false);
out = R.reliability('icc_type', '3-k');
assert(isstruct(out) && isfield(out, 'summary'), 'per-subject struct returned');
assert(out.summary.mean > 0.3, 'planted reliability should be substantial');
% replicate pool
icc = R.reliability('pool', 'replicate');
assert(isscalar(icc), 'replicate pool returns scalar');
end


% =========================================================================
function test_drift()
dat = synth('n_ses', 5);
R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'level', 'session', ...
    'subject_var', 'sub', 'session_var', 'sesno', 'metric', 'spearman', 'verbose', false);
% Single subject
idx = find(strcmp(R.replicate_table.sub, 's01'));
Rs = R; Rs.dat = R.dat(:,:,idx); Rs.replicate_table = R.replicate_table(idx,:);
[d_t, d_fit] = Rs.drift('reference', 'first', 'fit', 'linear');
assert(height(d_t) > 0, 'drift table populated');
assert(all(ismember({'condition','slope','P_value'}, d_fit.Properties.VariableNames)), 'fit cols');
% varargout with fit=none must not error
[d_t2, d_empty] = Rs.drift('reference', 'mean');
assert(istable(d_empty) && height(d_empty) == 0, 'empty fit table');
end


% =========================================================================
function test_lme()
dat = synth();
mdl = dat.rsa_lme('predictors', {'condition','bodysite','sesno'}, 'subject_var', 'sub', ...
    'verbose', false);
assert(isa(mdl, 'LinearMixedModel'), 'returns LME');
ce = mdl.Coefficients;
bc = ce.Estimate(strcmp(ce.Name, 'SameCondition'));
assert(bc > 0, 'SameCondition beta positive');
% icc + blups
icc = rsa_lme_icc(mdl);
assert(abs(sum(icc.summary.ICC) - 1) < 1e-6, 'ICC components sum to 1');
end


% =========================================================================
function test_compare_models()
dat = synth();
seq = rsa_model_sequence('Y ~ 1 + (1|Sub)');
seq = seq.add_term('SameCondition');
seq = seq.add_term('SameBodysite');
[T, best] = dat.rsa_compare_models(seq.formulas, ...
    'predictors', {'condition','bodysite'}, 'subject_var', 'sub', 'verbose', false);
assert(height(T) == 3, 'three models');
assert(best == 3, 'fullest model best by AIC');
assert(all(isfinite(T.lrt_p(2:end))), 'LRT p computed');
end


% =========================================================================
function test_compare()
dat = synth('n_sub', 9);
R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'metric', 'spearman', 'verbose', false);
result = R.compare({'condition','bodysite'}, 'correlation_type', 'kendall_taua', 'verbose', false);
assert(numel(result.candidate_names) == 2, 'two candidates');
assert(result.r_mean(1) > result.r_mean(2), 'condition (0.8) beats bodysite (0.4)');
assert(all(result.relatedness_sig), 'both planted models significant');
assert(result.noise_ceiling(2) >= result.noise_ceiling(1), 'upper >= lower NC');
end


% =========================================================================
function test_model_rdms()
dat = synth();
% from_categorical single + multi
M1 = rsm.from_categorical(dat.metadata_table, 'condition');
assert(numel(M1) == 1 && M1.is_dissimilarity, 'single categorical model');
M2 = rsm.from_categorical(dat.metadata_table, {'condition','bodysite'});
assert(numel(M2) == 2, 'two categorical models');
% from_metadata_distance
Md = rsm.from_metadata_distance(dat.metadata_table, 'sesno', 'metric', 'abs_diff');
assert(Md.is_dissimilarity, 'distance model is dissimilarity');
% from_design
X = [1 0 0; 0 1 0; 0 0 1];
Mx = rsm.from_design(X, 'names', {'p','q','r'});
assert(numel(Mx) == 3, 'three design models');
end
