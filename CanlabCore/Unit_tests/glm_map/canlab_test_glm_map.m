function tests = canlab_test_glm_map
%CANLAB_TEST_GLM_MAP Smoke/behavior tests for the glm_map object class.
%
% Covers construction, true Dependent property read-through, design
% diagnostics, direct-mode (2nd-level) fit over fmri_data.regress, event-mode
% build_design via the wrapped fmri_glm_design_matrix, import_SPM, and the
% main input-validation error paths.

tests = functiontests(localfunctions);
end


% =====================================================================
% Construction and Dependent properties
% =====================================================================
function test_empty_constructor(tc)
g = glm_map;
tc.verifyClass(g, 'glm_map');
tc.verifyEqual(g.level, 2);
tc.verifyFalse(g.is_timeseries);
tc.verifyFalse(g.is_fitted);
tc.verifyEqual(g.num_images, 0);
tc.verifyEqual(g.num_regressors, 0);
tc.verifyEqual(g.num_contrasts, 0);
end


function test_direct_design_dependent_properties(tc)
X = [ones(20, 1) zscore((1:20)') randn(20, 1)];
g = glm_map('X', X, 'level', 2, 'regressor_names', {'intercept', 'slope', 'nuis'});

tc.verifyEqual(g.num_images, 20);
tc.verifyEqual(g.num_regressors, 3);
tc.verifyEqual(g.X, X);                                  % read-through to Xdirect
tc.verifyEqual(g.regressor_names, {'intercept', 'slope', 'nuis'});
tc.verifyFalse(g.is_fitted);
end


function test_add_contrasts(tc)
X = [ones(20, 1) zscore((1:20)') randn(20, 1)];
g = glm_map('X', X, 'level', 2);
g = add_contrasts(g, [0 1 0; 0 0 1], {'slope_con', 'nuis_con'});

tc.verifyEqual(g.num_contrasts, 2);
tc.verifyEqual(size(g.contrasts), [3 2]);               % stored [regressors x contrasts]
tc.verifyEqual(g.contrast_names, {'slope_con', 'nuis_con'});
end


% =====================================================================
% Diagnostics (design only; no fit required)
% =====================================================================
function test_diagnostics_on_design(tc)
rng(0);
X = [ones(30, 1) zscore((1:30)') randn(30, 1)];
g = glm_map('X', X, 'level', 2);
g = add_contrasts(g, [0 1 0], {'slope'});
g = diagnostics(g, 'noverbose');

tc.verifyNumElements(g.diagnostics.Variance_inflation_factors, 3);
tc.verifyNumElements(g.diagnostics.Contrast_variance_inflation_factors, 1);
tc.verifyGreaterThan(g.diagnostics.condition_number, 0);
tc.verifyFalse(g.diagnostics.rank_deficient);
tc.verifyTrue(isstruct(g.diagnostics.collinearity_report));
end


function test_diagnostics_flags_rank_deficiency(tc)
X = [ones(20, 1) (1:20)' (1:20)'];                     % duplicate column
g = glm_map('X', X, 'level', 2);
w = warning('off', 'all');
c = onCleanup(@() warning(w));
g = diagnostics(g, 'noverbose');

tc.verifyTrue(g.diagnostics.rank_deficient);
tc.verifySize(g.diagnostics.collinearity_report.duplicate_column_pairs, [1 2]);
tc.verifyNotEmpty(g.warnings);
end


% =====================================================================
% Direct-mode (2nd-level / group) fit on real sample data
% =====================================================================
function test_direct_mode_fit(tc)
dat = canlab_get_sample_fmri_data();
n = size(dat.dat, 2);
g = glm_map('X', [ones(n, 1) zscore((1:n)')], 'level', 2, ...
    'regressor_names', {'intercept', 'cov'});
g = add_contrasts(g, [0 1], {'cov_effect'});
g = fit(g, dat, 'noverbose');

tc.verifyTrue(g.is_fitted);
tc.verifyEqual(g.dfe, n - 2);
tc.verifyClass(g.betas, 'statistic_image');
tc.verifyClass(g.t, 'statistic_image');
tc.verifyClass(g.contrast_estimates, 'statistic_image');
tc.verifyClass(g.contrast_t, 'statistic_image');
tc.verifyEqual(size(g.betas.dat, 2), 2);               % one beta image per regressor
tc.verifyEqual(size(g.contrast_estimates.dat, 2), 1);  % one image per contrast
tc.verifyNotEmpty(g.diagnostics.Variance_inflation_factors);
tc.verifyClass(g.df, 'fmri_data');                     % per-voxel df kept as fmri_data
tc.verifyClass(g.sigma, 'fmri_data');                  % per-voxel sigma kept as fmri_data
tc.verifyTrue(isstruct(g.input_parameters));           % nested option struct populated
end


function test_threshold_returns_glm_map(tc)
dat = canlab_get_sample_fmri_data();
n = size(dat.dat, 2);
g = glm_map('X', [ones(n, 1) zscore((1:n)')], 'level', 2);
g = fit(g, dat, 'noverbose');

w = warning('off', 'all');
c = onCleanup(@() warning(w));
g = threshold(g, .01, 'unc', 'which_map', 't');
tc.verifyClass(g, 'glm_map');
tc.verifyTrue(g.is_fitted);
end


% =====================================================================
% fmri_data.regress now returns a glm_map; constructor re-casts a struct
% =====================================================================
function test_regress_returns_glm_map(tc)
dat = canlab_get_sample_fmri_data();
n = size(dat.dat, 2);
dat.X = [zscore((1:n)') ones(n, 1)];
out = regress(dat, 0.05, 'unc', 'noverbose', 'nodisplay');

tc.verifyClass(out, 'glm_map');
% Historical struct-style field access still works via Dependent aliases
tc.verifyClass(out.b, 'statistic_image');              % alias for betas
tc.verifyClass(out.t, 'statistic_image');
tc.verifyClass(out.df, 'fmri_data');
tc.verifyClass(out.sigma, 'fmri_data');
tc.verifyTrue(isstruct(out.input_parameters));
tc.verifyTrue(isstruct(out.diagnostics));
tc.verifyNotEmpty(out.diagnostics.Variance_inflation_factors);
end


function test_construct_from_struct_and_aliases(tc)
% A struct with regress-style field names re-casts into a glm_map, and the
% out-struct aliases read the canonical properties.
S = struct();
S.b = statistic_image; S.b.dat = randn(50, 2);
S.t = statistic_image; S.t.dat = randn(50, 2);
S.variable_names = {'slope', 'intercept'};
S.C = [1 0]';
S.contrast_names = {'slope'};
S.analysis_name = 'from_struct_test';

g = glm_map(S);
tc.verifyClass(g, 'glm_map');
tc.verifyEqual(g.analysis_name, 'from_struct_test');
tc.verifyEqual(g.betas.dat, S.b.dat);                  % b -> betas
tc.verifyEqual(g.regressor_names, {'slope', 'intercept'});  % variable_names -> regressor_names
tc.verifyEqual(g.contrasts, [1 0]');                   % C -> contrasts
tc.verifyEqual(g.contrast_names, {'slope'});
% Round-trip: alias getters return the canonical values
tc.verifyEqual(g.b.dat, g.betas.dat);
tc.verifyEqual(g.variable_names, g.regressor_names);
end


% =====================================================================
% Event-mode (1st-level) design build via wrapped fmri_glm_design_matrix
% =====================================================================
function test_event_mode_build_design(tc)
TR = 2; nscan = 120;
ons = {[10 40 70 100]', [25 55 85 115]'};
d = fmri_glm_design_matrix(TR, 'nscan', nscan, 'units', 'secs', ...
    'onsets', ons, 'condition_names', {'A', 'B'});
g = glm_map(d);

tc.verifyEqual(g.level, 1);                            % wrapping a design implies 1st-level
tc.verifyEqual(g.TR, TR);                              % read-through to design.TR
tc.verifyEmpty(g.X);                                   % not built yet

w = warning('off', 'all');
c = onCleanup(@() warning(w));
g = build_design(g);

tc.verifyEqual(size(g.X, 1), nscan);                   % built design matrix
% Default canonical HRF: one column per condition + one intercept = 3
tc.verifyEqual(size(g.X, 2), 3);
tc.verifyNotEmpty(g.regressor_names);
tc.verifyNumElements(g.onsets, 2);                     % onsets read-through to design

% Timing: condition A's regressor should peak ~3 TRs after its first onset
% (canonical HRF peak ~6 s = 3 TRs); first onset 10 s -> TR 6 -> peak ~TR 9.
[~, pk] = max(g.X(:, 1));
tc.verifyLessThanOrEqual(abs(pk - 9), 1);
end


function test_event_design_matches_onsets2fmridesign(tc)
% The high-resolution pipeline should match onsets2fmridesign across TRs,
% including fractional TRs.
ons = {[10 40 70 100 150]', [25 55 85 130 200]'};
for TR = [2 1.3 2.5]
    nscan = round(260 / TR);
    d = fmri_glm_design_matrix(TR, 'nscan', nscan, 'units', 'secs', ...
        'onsets', ons, 'condition_names', {'A', 'B'});
    w = warning('off', 'all'); c = onCleanup(@() warning(w)); %#ok<NASGU>
    g = glm_map(d); g = build_design(g);
    Xgt = onsets2fmridesign(ons, TR, nscan * TR, 'hrf');   % intercept last col
    tc.verifyEqual(size(g.X, 1), size(Xgt, 1));
    tc.verifyGreaterThan(corr(g.X(:, 1), Xgt(:, 1)), 0.999);
    tc.verifyGreaterThan(corr(g.X(:, 2), Xgt(:, 2)), 0.999);
end
end


function test_replace_basis_set(tc)
TR = 2; nscan = 120;
d = fmri_glm_design_matrix(TR, 'nscan', nscan, 'units', 'secs', ...
    'onsets', {[10 40 70 100]' [25 55 85 115]'}, 'condition_names', {'A', 'B'});
w = warning('off', 'all'); c = onCleanup(@() warning(w)); %#ok<NASGU>

g = glm_map(d); g.is_timeseries = true; g = build_design(g);
rng(0); sim = fmri_data; sim.dat = randn(50, nscan);
g = add_contrasts(g, [1 -1 0], {'AmB'});
g = fit(g, sim, 'noverbose');
tc.verifyEqual(g.num_regressors, 3);
tc.verifyTrue(g.is_fitted);

[xBF_hires, ~] = fmri_spline_basis(TR, 'length', 20, 'nbasis', 4, 'order', 3);

% Without data: rebuild, clear stale fit/contrasts, refresh diagnostics
g2 = replace_basis_set(g, 1, xBF_hires);
tc.verifyEqual(g2.num_regressors, 6);          % 4 spline (A) + 1 hrf (B) + intercept
tc.verifyFalse(g2.is_fitted);                  % fitted maps cleared
tc.verifyEqual(g2.num_contrasts, 0);           % contrasts cleared (columns changed)
tc.verifyNotEmpty(g2.diagnostics.Variance_inflation_factors);  % diagnostics refreshed

% With data: re-fit on the new design
g3 = replace_basis_set(g, 1, xBF_hires, 'data', sim, 'noverbose');
tc.verifyTrue(g3.is_fitted);
tc.verifyEqual(size(g3.betas.dat, 2), 6);
end


% =====================================================================
% import_SPM (SPM12/SPM25 first-level)
% =====================================================================
function test_import_SPM_and_fit(tc)
rng(1);
SPM = local_make_synthetic_spm();
g = import_SPM(glm_map, SPM, 'noverbose');

tc.verifyEqual(g.level, 1);
tc.verifyTrue(g.is_timeseries);
tc.verifyEqual(g.TR, SPM.xY.RT);
tc.verifyEqual(size(g.X), size(SPM.xX.X));
tc.verifyEqual(g.regressor_names, SPM.xX.name);
tc.verifyNumElements(g.onsets, 2);                     % onsets read-through
tc.verifyEqual([g.condition_names{:}], {'Cue', 'Pain'});

% Fit against a matching synthetic timeseries
nscan = size(SPM.xX.X, 1);
sim = fmri_data;
sim.dat = randn(25, nscan);
g = fit(g, sim, 'noverbose');

tc.verifyTrue(g.is_fitted);
tc.verifyEqual(g.dfe, nscan - size(SPM.xX.X, 2));
end


% =====================================================================
% Input-validation error paths
% =====================================================================
function test_fit_without_design_errors(tc)
g = glm_map;                                            % no X, no design
dat = fmri_data; dat.dat = randn(10, 5);
tc.verifyError(@() fit(g, dat, 'noverbose'), 'glm_map:NoDesign');
end


function test_AR_requires_timeseries_errors(tc)
X = [ones(12, 1) randn(12, 1)];
g = glm_map('X', X, 'level', 2);                       % is_timeseries == false
dat = fmri_data; dat.dat = randn(8, 12);
tc.verifyError(@() fit(g, dat, 'AR', 1, 'noverbose'), 'glm_map:ARnotTimeseries');
end


function test_add_contrasts_size_mismatch_errors(tc)
X = [ones(15, 1) randn(15, 2)];                        % 3 regressors
g = glm_map('X', X, 'level', 2);
tc.verifyError(@() add_contrasts(g, [1 0], {'bad'}), 'glm_map:ContrastSize');
end


% =====================================================================
% Local helpers
% =====================================================================
function SPM = local_make_synthetic_spm()
% Minimal SPM12/SPM25-shared first-level structure for import tests.
nscan = 60;
TR = 1.5;
SPM = struct();
SPM.xY.RT = TR;
SPM.nscan = nscan;
SPM.xBF = struct('name', 'hrf', 'length', 32, 'order', 1, 'dt', TR / 16, ...
    'T', 16, 'T0', 1, 'UNITS', 'secs', 'Volterra', 1);

P = struct('name', 'none', 'P', [], 'h', 0);
U(1) = struct('dt', TR / 16, 'name', {{'Cue'}},  'ons', [12 30 48]', 'dur', [0 0 0]',   'P', P, 'u', [], 'pst', []);
U(2) = struct('dt', TR / 16, 'name', {{'Pain'}}, 'ons', [18 36 54]', 'dur', [6 6 6]',   'P', P, 'u', [], 'pst', []);
SPM.Sess = struct('U', U, 'C', struct('C', [], 'name', {{}}), 'row', 1:nscan, 'col', 1:3, 'Fc', []);

SPM.xX.X = [randn(nscan, 2) ones(nscan, 1)];
SPM.xX.name = {'Sn(1) Cue*bf(1)', 'Sn(1) Pain*bf(1)', 'Sn(1) constant'};
SPM.xX.iH = [1 2]; SPM.xX.iC = []; SPM.xX.iB = 3; SPM.xX.iG = [];
end
