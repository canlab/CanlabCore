function tests = canlab_test_xval_SVR
%CANLAB_TEST_XVAL_SVR Cross-validated SVR regression via xval_SVR.
%
% Converted from Unit_tests/xval_SVR_unit_test.m. Predicts a synthetic
% continuous outcome from the DPSP Hot-Warm contrast maps via cross-validated
% linear SVR and checks the returned @predictive_model object. Model fit once
% in setupOnce; skipped if the DPSP sample data is not on the path.

tests = functiontests(localfunctions);
end


function setupOnce(tc)   %#ok<*DEFNU>
[hot, warm, ok] = canlab_get_dpsp_hot_warm();
tc.assumeTrue(ok, 'DPSP sample data not on path');

hot_vs_warm = image_math(hot, warm, 'minus');

rng(0);
[p, n] = size(hot_vs_warm.dat);
keep_vox = randsample(p, min(3000, p));
X = double(hot_vs_warm.dat(keep_vox, :))';

% Synthetic outcome predictable from a sparse subset of features.
b_true = zeros(size(X, 2), 1);
b_true(1:30) = randn(30, 1);
Y  = X * b_true + 0.5 * std(X * b_true) * randn(n, 1);
id = (1:n)';

tc.TestData.X  = X;
tc.TestData.Y  = Y;
tc.TestData.pm = xval_SVR(X, Y, id, 'nooptimize', 'norepeats', ...
    'nobootstrap', 'noverbose', 'noplot');
end


function test_class_and_state(tc)
pm = tc.TestData.pm;
tc.verifyClass(pm, 'predictive_model');
tc.verifyTrue(pm.is_fitted, 'is_fitted should be true');
tc.verifyTrue(pm.is_regressor, 'continuous Y should give is_regressor true');
end


function test_canonical_path_fields_populated(tc)
pm = tc.TestData.pm;
tc.verifyNotEmpty(pm.Y);
tc.verifyNotEmpty(pm.fitted_values.yfit);
tc.verifyNotEmpty(pm.weights.w);
tc.verifyNotEmpty(pm.error_metrics.prediction_outcome_r.value);
tc.verifyNotEmpty(pm.ml_model);
end


function test_weight_shape(tc)
pm = tc.TestData.pm;
tc.verifyEqual(size(pm.weights.w, 1), size(tc.TestData.X, 2), ...
    'weights.w should have one row per feature');
end


function test_fit_type_and_omitted_markers(tc)
pm = tc.TestData.pm;
tc.verifyEqual(pm.fit_type, 'crossval');
tc.verifyClass(pm.omitted_cases, 'logical');
tc.verifyClass(pm.omitted_features, 'logical');
end


function test_validate_object_accepts(tc)
tc.TestData.pm.validate_object('noverbose');
end
