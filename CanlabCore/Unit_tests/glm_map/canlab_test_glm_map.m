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

tc.verifyNumElements(g.vif, 3);
tc.verifyNumElements(g.contrast_vif, 1);
tc.verifyGreaterThan(g.condition_number, 0);
tc.verifyFalse(g.rank_deficient);
tc.verifyTrue(isstruct(g.collinearity_report));
end


function test_diagnostics_flags_rank_deficiency(tc)
X = [ones(20, 1) (1:20)' (1:20)'];                     % duplicate column
g = glm_map('X', X, 'level', 2);
w = warning('off', 'all');
c = onCleanup(@() warning(w));
g = diagnostics(g, 'noverbose');

tc.verifyTrue(g.rank_deficient);
tc.verifySize(g.collinearity_report.duplicate_column_pairs, [1 2]);
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
tc.verifyNotEmpty(g.vif);
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
tc.verifyGreaterThan(size(g.X, 2), 0);
tc.verifyNotEmpty(g.regressor_names);
tc.verifyNumElements(g.onsets, 2);                     % onsets read-through to design
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
