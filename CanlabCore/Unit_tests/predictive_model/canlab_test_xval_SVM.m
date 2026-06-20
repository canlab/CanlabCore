function tests = canlab_test_xval_SVM
%CANLAB_TEST_XVAL_SVM Cross-validated SVM classification via xval_SVM.
%
% Converted from Unit_tests/xval_SVM_unit_test.m. Drives xval_SVM on the DPSP
% Hot vs Warm single-subject maps (a paired, within-person design) and checks
% the returned @predictive_model object: class/state, canonical-path fields,
% shapes/ranges, fit_type + omitted markers, validate_object, and clone.
%
% The model is fit once in setupOnce and cached; each test reads the cached
% object. Skipped if the DPSP sample data is not on the path.

tests = functiontests(localfunctions);
end


function setupOnce(tc)   %#ok<*DEFNU>
[hot, warm, ok] = canlab_get_dpsp_hot_warm();
tc.assumeTrue(ok, 'DPSP sample data not on path');

n_hot  = size(hot.dat, 2);
n_warm = size(warm.dat, 2);

rng(0);
p = size(hot.dat, 1);
keep_vox = randsample(p, min(5000, p));
X  = double([hot.dat(keep_vox, :) warm.dat(keep_vox, :)])';
Y  = [ones(n_hot, 1); -ones(n_warm, 1)];
id = [(1:n_hot)'; (1:n_warm)'];   % each subject contributes both conditions

tc.TestData.X  = X;
tc.TestData.Y  = Y;
tc.TestData.id = id;
tc.TestData.pm = xval_SVM(X, Y, id, 'nooptimize', 'norepeats', ...
    'nobootstrap', 'noverbose', 'noplot');
end


function test_class_and_state(tc)
pm = tc.TestData.pm;
tc.verifyClass(pm, 'predictive_model');
tc.verifyTrue(pm.is_fitted, 'is_fitted should be true');
tc.verifyTrue(pm.is_classifier, 'binary Y should give is_classifier true');
end


function test_canonical_path_fields_populated(tc)
pm = tc.TestData.pm;
tc.verifyNotEmpty(pm.Y);
tc.verifyNotEmpty(pm.id);
tc.verifyNotEmpty(pm.fitted_values.yfit);
tc.verifyNotEmpty(pm.fitted_values.dist_from_hyperplane_xval);
tc.verifyNotEmpty(pm.weights.w);
tc.verifyNotEmpty(pm.error_metrics.crossval_accuracy.value);
tc.verifyNotEmpty(pm.error_metrics.d_singleinterval.value);
tc.verifyNotEmpty(pm.ml_model);
tc.verifyNotEmpty(pm.cv_partition.trIdx);
tc.verifyNotEmpty(pm.cv_partition.teIdx);
tc.verifyNotEmpty(pm.cv_partition.nfolds);
end


function test_shapes_and_ranges(tc)
pm = tc.TestData.pm;
n  = numel(tc.TestData.Y);
nfeat = size(tc.TestData.X, 2);

tc.verifyNumElements(pm.fitted_values.yfit, n);
tc.verifyNumElements(pm.fitted_values.dist_from_hyperplane_xval, n);
tc.verifyEqual(size(pm.weights.w, 1), nfeat, ...
    'weights.w should have one row per feature');

cv_acc = pm.error_metrics.crossval_accuracy.value;
tc.verifyGreaterThanOrEqual(cv_acc, 0);
tc.verifyLessThanOrEqual(cv_acc, 100);

% DPSP is paired, so within-person scoring should be populated.
tc.verifyFalse(isnan(pm.error_metrics.crossval_accuracy_within.value));
tc.verifyFalse(isnan(pm.error_metrics.d_within.value));
tc.verifyTrue(pm.diagnostics.mult_obs_within_person, ...
    'should detect multiple observations per id');
end


function test_fit_type_and_omitted_markers(tc)
pm = tc.TestData.pm;
tc.verifyEqual(pm.fit_type, 'crossval');
tc.verifyClass(pm.omitted_cases, 'logical');
tc.verifyClass(pm.omitted_features, 'logical');
tc.verifyNumElements(pm.omitted_cases, numel(tc.TestData.Y));
tc.verifyNumElements(pm.omitted_features, size(tc.TestData.X, 2));
end


function test_validate_object_accepts(tc)
% validate_object throws on an invalid object; reaching the end is a pass.
tc.TestData.pm.validate_object('noverbose');
end


function test_clone_clears_fitted_state(tc)
pm2 = clone(tc.TestData.pm);
tc.verifyFalse(pm2.is_fitted, 'clone should clear fitted state');
tc.verifyEqual(pm2.modeloptions, tc.TestData.pm.modeloptions, ...
    'clone should preserve modeloptions');
end
