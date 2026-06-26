function tests = canlab_test_xval_regression_multisubject_featureselect
%CANLAB_TEST_XVAL_REGRESSION_MULTISUBJECT_FEATURESELECT Feature-select regression.
%
% Converted from Unit_tests/xval_regression_multisubject_featureselect_unit_test.m.
% Runs xval_regression_multisubject_featureselect ('ols') on the DPSP Hot-Warm
% contrast maps against a synthetic continuous outcome. Beyond the usual
% predictive_model checks, it pins the per-fold weight padding: weights.w_perfold
% must span the FULL feature space, not the selected-feature count (the bug this
% test originally guarded). Skipped if the DPSP sample data is not on the path.

tests = functiontests(localfunctions);
end


function setupOnce(tc)   %#ok<*DEFNU>
[hot, warm, ok] = canlab_get_dpsp_hot_warm();
tc.assumeTrue(ok, 'DPSP sample data not on path');

hvw = image_math(hot, warm, 'minus');

rng(0);
[p, n] = size(hvw.dat);
keep_vox = randsample(p, min(2000, p));
X = double(hvw.dat(keep_vox, :))';

b_true = zeros(size(X, 2), 1);
b_true(1:30) = randn(30, 1);
Y = X * b_true + 0.5 * std(X * b_true) * randn(n, 1);

tc.TestData.X = X;
tc.TestData.Y = Y;
tc.TestData.pm = xval_regression_multisubject_featureselect('ols', {Y}, {X}, ...
    'holdout_method', 'balanced4', 'noverbose');
end


function test_class_and_state(tc)
pm = tc.TestData.pm;
tc.verifyClass(pm, 'predictive_model');
tc.verifyTrue(pm.is_fitted, 'is_fitted should be true');
end


function test_canonical_path_fields_populated(tc)
% mean_vox_weights is only populated when pcsquash is enabled; this no-PCA
% test does not exercise that path, so it is not checked here.
pm = tc.TestData.pm;
tc.verifyNotEmpty(pm.fitted_values.subjfit);
tc.verifyNotEmpty(pm.weights.w_perfold);
tc.verifyNotEmpty(pm.error_metrics.pred_err.value);
tc.verifyNotEmpty(pm.error_metrics.pred_err_null.value);
tc.verifyNotEmpty(pm.error_metrics.var_reduction.value);
tc.verifyNotEmpty(pm.Y);
end


function test_w_perfold_spans_full_feature_space(tc)
% Regression guard: w_perfold must have one row per ORIGINAL feature, not per
% selected feature.
pm = tc.TestData.pm;
tc.verifyEqual(size(pm.weights.w_perfold, 1), size(tc.TestData.X, 2), ...
    'weights.w_perfold should span the full feature count');
end


function test_variance_explained(tc)
r2 = tc.TestData.pm.error_metrics.r_squared.value;
tc.verifyNotEmpty(r2);
tc.verifyGreaterThan(r2(1), 0, 'expected r_squared > 0');
end


function test_validate_object_accepts(tc)
tc.TestData.pm.validate_object('noverbose');
end
