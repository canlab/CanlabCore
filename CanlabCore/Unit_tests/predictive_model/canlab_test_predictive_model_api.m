function tests = canlab_test_predictive_model_api
%CANLAB_TEST_PREDICTIVE_MODEL_API End-to-end test of the @predictive_model public API.
%
% Converted from Unit_tests/predictive_model_unit_test.m. Exercises the
% sklearn-style @predictive_model object surface on the DPSP Hot vs Warm
% dataset: construction, fit/predict/score, crossval (cv_splitter), bootstrap,
% weight_map_object, permutation_test, calibrate/predict_proba,
% select_features, grid_search, stability_selection, @pipeline,
% fmri_data.predict 'newapi' routing, and the regression/report/summary path.
% Visualisation methods (plot/confusionchart/montage) are smoke-tested and
% skipped on a headless runner. Complements the xval_* wrapper tests, which
% cover the legacy entry points rather than the object's own methods.
%
% The shared fit -> crossval -> bootstrap chain is computed once in setupOnce
% and cached; method-specific tests build their own small models. Small
% nboot/nperm/grid values keep this to a few minutes; it is still much slower
% than the rest of the suite. Skipped if DPSP_hotwarm is not on the path.

tests = functiontests(localfunctions);
end


% =====================================================================
% Fixtures
% =====================================================================

function setupOnce(tc)   %#ok<*DEFNU>
tc.assumeNotEmpty(which('load_image_set'), 'load_image_set not on path');
try
    hw_obj = load_image_set('DPSP_hotwarm', 'noverbose');
catch ME
    tc.assumeFail(['DPSP_hotwarm sample not available: ' ME.message]);
    return
end

% Restrict to a sparse gray-matter ROI so the many cross-validations below
% are fast: mask to gray matter, then deterministically thin to every 8th
% voxel (~20k features instead of ~195k). hw_obj and X stay in lockstep, so
% weight_map_object / montage still back-project consistently. Falls back to
% the full volume if the gray-matter mask is not on the path.
gm_file = which('gray_matter_mask.img');
if ~isempty(gm_file)
    hw_obj = remove_empty(apply_mask(hw_obj, gm_file));
    thin = get_wh_image(hw_obj, 1);
    thin.dat = zeros(size(hw_obj.dat, 1), 1);
    thin.dat(1:8:end) = 1;
    hw_obj = remove_empty(apply_mask(hw_obj, thin));
end

X  = double(hw_obj.dat');
Y  = hw_obj.Y;
id = grp2idx(hw_obj.metadata_table.subj_id);

tc.TestData.hw_obj = hw_obj;
tc.TestData.X  = X;
tc.TestData.Y  = Y;
tc.TestData.id = id;

% Shared expensive chain, computed once: in-sample fit, cross-validation,
% and bootstrap on the cross-validated model.
tc.TestData.pm_insample = fit( ...
    predictive_model('algorithm', 'svm', 'task', 'classification', 'random_state', 0), ...
    X, Y, 'id', id);

pm_cv = predictive_model('algorithm', 'svm', 'task', 'classification', ...
    'modeloptions', {'KernelFunction', 'linear'}, 'random_state', 0);
pm_cv = crossval(pm_cv, X, Y, ...
    'cv', cv_splitter.stratified_group_kfold(5), 'groups', id, ...
    'scoring', 'balanced_accuracy');
tc.TestData.pm_cv = pm_cv;

% Minimum-viable nboot: just enough to compute boot std / z / p and exercise
% the empirical p-floor. This is a smoke test, not real inference.
tc.TestData.pm_boot = bootstrap(pm_cv, X, Y, 'nboot', 5, 'groups', id, 'verbose', false);
end


function setup(tc)
tc.TestData.PrevFigVis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
end


function teardown(tc)
close all force
if isfield(tc.TestData, 'PrevFigVis')
    set(0, 'DefaultFigureVisible', tc.TestData.PrevFigVis);
end
end


% =====================================================================
% Construction / in-sample fit / predict / score
% =====================================================================

function test_construction_hyperparameters(tc)
pm = predictive_model('algorithm', 'svm', 'task', 'classification', 'random_state', 0);
tc.verifyEqual(pm.algorithm, 'svm');
tc.verifyEqual(pm.task, 'classification');
tc.verifyFalse(pm.is_fitted, 'should start unfitted');
tc.verifyTrue(pm.is_classifier);
end


function test_fit_insample(tc)
pm = tc.TestData.pm_insample;
tc.verifyTrue(pm.is_fitted);
tc.verifyEqual(pm.fit_type, 'insample');
tc.verifyLessThanOrEqual(numel(pm.weights.w), size(tc.TestData.X, 2));
end


function test_predict(tc)
[yhat, scores] = predict(tc.TestData.pm_insample, tc.TestData.X);
tc.verifyNumElements(yhat, numel(tc.TestData.Y));
tc.verifyEqual(size(scores, 1), numel(tc.TestData.Y));
end


function test_score_balanced_accuracy(tc)
pm = tc.TestData.pm_insample;
pm.scorer = cv_scorer.balanced_accuracy();
v = score(pm, tc.TestData.X, tc.TestData.Y);
tc.verifyGreaterThanOrEqual(v, 0);
tc.verifyLessThanOrEqual(v, 1);
end


% =====================================================================
% crossval / bootstrap / weight map
% =====================================================================

function test_crossval_stratified_group_kfold(tc)
pm = tc.TestData.pm_cv;
tc.verifyEqual(pm.fit_type, 'crossval');
tc.verifyEqual(pm.cv_partition.nfolds, 5);
tc.verifyNumElements(pm.fitted_values.yfit, numel(tc.TestData.Y));
tc.verifyFalse(isnan(pm.error_metrics.balanced_accuracy.value), 'metric populated');
end


function test_bootstrap_weight_stats(tc)
nboot = 5;   % must match setupOnce
pm = tc.TestData.pm_boot;
tc.verifyEqual(size(pm.weights.boot_w, 2), nboot, 'boot_w should have nboot columns');
tc.verifyNumElements(pm.weights.z, numel(pm.weights.w));
tc.verifyNumElements(pm.weights.p, numel(pm.weights.w));
tc.verifyTrue(all(pm.weights.p >= 0 & pm.weights.p <= 1), 'p in [0,1]');
tc.verifyTrue(all(abs(pm.weights.z) < 1e6), 'z floored, not 1e15');
p_floor = 2 / (nboot + 1);
tc.verifyGreaterThanOrEqual(min(pm.weights.p), p_floor - eps, 'empirical p floor');
end


function test_weight_map_object(tc)
pm = tc.TestData.pm_boot;
hw = tc.TestData.hw_obj;
[pm_w, si_full] = weight_map_object(pm, hw);
[~,    si_thr]  = weight_map_object(pm, hw, 'use', 'thresh_fdr');
tc.verifyClass(si_full, 'statistic_image');
tc.verifyEqual(size(si_full.dat, 1), size(hw.dat, 1));
tc.verifyEqual(size(si_thr.dat, 1), size(hw.dat, 1));
tc.verifyNotEmpty(pm_w.weights.weight_obj, 'weight_obj cached on pm');
tc.verifyClass(pm_w.weights.weight_obj, 'statistic_image');
end


% =====================================================================
% permutation / calibration / feature selection
% =====================================================================

function test_permutation_test_within_subjects(tc)
pm = predictive_model('algorithm', 'svm', 'task', 'classification', ...
    'modeloptions', {'KernelFunction', 'linear'}, 'random_state', 0);
pm.cv = cv_splitter.stratified_group_kfold(3);
nperm = 3;   % minimum-viable: enough to form a null and a defined p-value
pm = permutation_test(pm, tc.TestData.X, tc.TestData.Y, ...
    'nperm', nperm, 'groups', tc.TestData.id, 'verbose', false);
tc.verifyNumElements(pm.permutation_results.null_scores, nperm);
tc.verifyGreaterThan(pm.permutation_results.p_value, 0);
tc.verifyLessThanOrEqual(pm.permutation_results.p_value, 1);
% DPSP is paired (Hot+Warm per subject) -> auto should pick within_subjects.
tc.verifyEqual(pm.permutation_results.permutation, 'within_subjects');
end


function test_calibrate_predict_proba(tc)
pm_cal = calibrate(tc.TestData.pm_cv, tc.TestData.X, tc.TestData.Y);
tc.verifyTrue(isfield(pm_cal.fitted_values, 'calibrator'), 'calibrator present');
P = predict_proba(pm_cal, tc.TestData.X);
tc.verifyTrue(all(P >= 0 & P <= 1), 'predict_proba in [0,1]');
end


function test_select_features(tc)
pm = predictive_model('algorithm', 'svm', 'task', 'classification');
pm = select_features(pm, tc.TestData.X, tc.TestData.Y, 'k', 500, 'verbose', false);
tc.verifyEqual(pm.diagnostics.feature_selection.n_selected, 500);
end


% =====================================================================
% grid search / stability selection / pipeline / newapi
% =====================================================================

function test_grid_search(tc)
pm = predictive_model('algorithm', 'svm', 'task', 'classification', ...
    'cv', cv_splitter.stratified_group_kfold(3));
pm = grid_search(pm, tc.TestData.X, tc.TestData.Y, ...
    struct('BoxConstraint', [0.1 1]), 'groups', tc.TestData.id, 'verbose', false);
tc.verifyTrue(isfield(pm.diagnostics.grid_search, 'best_score'));
end


function test_stability_selection(tc)
pm = stability_selection( ...
    predictive_model('algorithm', 'linear_svm', 'task', 'classification'), ...
    tc.TestData.X, tc.TestData.Y, 'nboot', 5, 'k', 1000, 'threshold', 0.6, ...
    'groups', tc.TestData.id, 'verbose', false);
ss = pm.diagnostics.stability_selection;
tc.verifyNumElements(ss.selection_freq, size(tc.TestData.X, 2));
tc.verifyTrue(all(ss.selection_freq >= 0 & ss.selection_freq <= 1), 'freq in [0,1]');
end


function test_pipeline_pca_svm(tc)
est  = predictive_model('algorithm', 'svm', 'task', 'classification');
pipe = pipeline({ {'pca', 'k', 20} }, est);
pipe = crossval(pipe, tc.TestData.X, tc.TestData.Y, ...
    'groups', tc.TestData.id, 'cv', cv_splitter.stratified_group_kfold(5));
tc.verifyEqual(pipe.fit_type, 'crossval');
tc.verifyNumElements(pipe.weights.w, size(tc.TestData.X, 2), ...
    'pipeline should back-project weights to full feature space');
si_pipe = weight_map_object(pipe, tc.TestData.hw_obj);
tc.verifyClass(si_pipe, 'statistic_image');
end


function test_fmri_data_predict_newapi(tc)
[cverr, st, ~, pm_api] = predict(tc.TestData.hw_obj, 'algorithm_name', 'cv_svm', ...
    'nfolds', 5, 'newapi', 'verbose', 0);
tc.verifyClass(pm_api, 'predictive_model');
tc.verifyEqual(pm_api.fit_type, 'crossval');
tc.verifyTrue(isfield(st, 'dist_from_hyperplane_xval'));
tc.verifyTrue(isfinite(cverr));
end


% =====================================================================
% regression metrics / report_accuracy / summary
% =====================================================================

function test_regression_metrics_and_reports(tc)
rng(7);
Y  = tc.TestData.Y;
Yr = double(Y) + 0.5 * randn(size(Y));
pm_reg = predictive_model('algorithm', 'pcr', 'task', 'regression');
pm_reg = crossval(pm_reg, tc.TestData.X, Yr, ...
    'cv', cv_splitter.group_kfold(5), 'groups', tc.TestData.id);

tc.verifyTrue(isfield(pm_reg.error_metrics, 'predicted_r2'));
tc.verifyTrue(isfield(pm_reg.error_metrics, 'out_of_sample_r2'));

% predicted_r2 must equal 1 - PRESS/SST computed by hand.
yf = pm_reg.fitted_values.yfit(:);
yo = Yr(:);
pr2_manual = 1 - sum((yo - yf).^2) / sum((yo - mean(yo)).^2);
tc.verifyEqual(pm_reg.error_metrics.predicted_r2.value, pr2_manual, 'AbsTol', 1e-9);

acc_r = report_accuracy(pm_reg, 'noverbose');
tc.verifyEqual(acc_r.task, 'regression');
tc.verifyFalse(isnan(acc_r.predicted_r2));

s_reg = summary(pm_reg, 'noverbose');
tc.verifyEqual(s_reg.provenance.fit_type, 'crossval');

% Classification report_accuracy: ROC-derived fields populated.
acc_c = report_accuracy(tc.TestData.pm_cv, 'noverbose');
tc.verifyEqual(acc_c.task, 'classification');
tc.verifyFalse(isnan(acc_c.auc));
tc.verifyFalse(isnan(acc_c.npv));
end


% =====================================================================
% Visualisation (smoke; skipped on a headless runner)
% =====================================================================

function test_rocplot_metrics(tc)
% rocplot('noplot') returns metrics without a display.
ROC = rocplot(tc.TestData.pm_boot, 'noplot');
tc.verifyTrue(isfield(ROC, 'AUC'));
tc.verifyGreaterThanOrEqual(ROC.AUC, 0);
tc.verifyLessThanOrEqual(ROC.AUC, 1);
end


function test_display_methods_smoke(tc)
pm = tc.TestData.pm_boot;
hw = tc.TestData.hw_obj;
try
    h = plot(pm); close(h);             % classification dispatcher (violin + ROC)
    confusionchart(pm); close(gcf);
    montage(pm, hw); close all force;   % delegates to weight-map montage
catch ME
    if strcmp(canlab_classify_environment_error(ME), 'genuine')
        rethrow(ME);
    end
    tc.assumeFail(['display method needs a graphics environment: ' ME.message]);
end
end
