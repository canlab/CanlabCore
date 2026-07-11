function tests = canlab_test_xval_regression_multisubject
%CANLAB_TEST_XVAL_REGRESSION_MULTISUBJECT OLS+PCA multisubject regression.
%
% Converted from Unit_tests/xval_regression_multisubject_unit_test.m. Runs
% xval_regression_multisubject ('ols' + 'pca') on the DPSP Hot-Warm contrast
% maps against a synthetic continuous outcome and checks the returned
% @predictive_model object. Model fit once in setupOnce; skipped if the DPSP
% sample data is not on the path.

tests = functiontests(localfunctions);
end


function setupOnce(tc)   %#ok<*DEFNU>
[hot, warm, ok] = canlab_get_dpsp_hot_warm();
tc.assumeTrue(ok, 'DPSP sample data not on path');

hot_vs_warm = image_math(hot, warm, 'minus');

rng(0);
[p, n] = size(hot_vs_warm.dat);
keep_vox = randsample(p, min(5000, p));
X = double(hot_vs_warm.dat(keep_vox, :))';        % n x v

b_true = zeros(size(X, 2), 1);
b_true(1:50) = randn(50, 1);
Y = X * b_true + 0.5 * std(X * b_true) * randn(n, 1);

tc.TestData.X = X;
tc.TestData.Y = Y;
tc.TestData.n = n;
tc.TestData.pm = xval_regression_multisubject('ols', {Y}, {X}, ...
    'pca', 'ndims', 10, 'holdout_method', 'balanced4', 'noverbose');
end


function test_class_and_state(tc)
pm = tc.TestData.pm;
tc.verifyClass(pm, 'predictive_model');
tc.verifyTrue(pm.is_fitted, 'is_fitted should be true');
end


function test_canonical_path_fields_populated(tc)
pm = tc.TestData.pm;
tc.verifyNotEmpty(pm.fitted_values.subjfit);
tc.verifyNotEmpty(pm.weights.subjbetas);
tc.verifyNotEmpty(pm.weights.mean_vox_weights);
tc.verifyNotEmpty(pm.error_metrics.pred_err.value);
tc.verifyNotEmpty(pm.error_metrics.pred_err_null.value);
tc.verifyNotEmpty(pm.error_metrics.var_reduction.value);
tc.verifyNotEmpty(pm.error_metrics.r_each_subject.value);
tc.verifyNotEmpty(pm.Y);
tc.verifyNotEmpty(pm.inputParameters);
end


function test_shapes_and_variance_explained(tc)
pm = tc.TestData.pm;
tc.verifyNumElements(pm.fitted_values.subjfit, 1, ...
    'subjfit cell length should equal number of datasets');
tc.verifyNumElements(pm.fitted_values.subjfit{1}, tc.TestData.n);
tc.verifyEqual(size(pm.weights.mean_vox_weights, 1), size(tc.TestData.X, 2), ...
    'mean_vox_weights should have one row per voxel');

r2 = pm.error_metrics.r_squared.value;
tc.verifyNotEmpty(r2);
tc.verifyGreaterThan(r2(1), 0, 'OLS+PCA should explain non-trivial variance');
end


function test_validate_object_accepts(tc)
tc.TestData.pm.validate_object('noverbose');
end
