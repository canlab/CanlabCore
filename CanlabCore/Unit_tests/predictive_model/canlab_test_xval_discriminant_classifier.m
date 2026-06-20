function tests = canlab_test_xval_discriminant_classifier
%CANLAB_TEST_XVAL_DISCRIMINANT_CLASSIFIER Cross-validated LDA classification.
%
% Converted from Unit_tests/xval_discriminant_classifier_unit_test.m. Drives
% cross-validated LDA (fitcdiscr) on the DPSP Hot vs Warm single-subject maps
% (100 random voxels, so n > p) and checks the returned @predictive_model
% object. Model fit once in setupOnce; skipped if DPSP data is not on the path.

tests = functiontests(localfunctions);
end


function setupOnce(tc)   %#ok<*DEFNU>
[hot, warm, ok] = canlab_get_dpsp_hot_warm();
tc.assumeTrue(ok, 'DPSP sample data not on path');

rng(0);
p = size(hot.dat, 1);
keep_vox = randsample(p, 100);   % LDA needs n > p
X = double([hot.dat(keep_vox, :) warm.dat(keep_vox, :)])';
labels = int32([ones(size(hot.dat, 2), 1); 2 * ones(size(warm.dat, 2), 1)]);

tc.TestData.pm = xval_discriminant_classifier(X, labels, 'nFolds', 5, ...
    'verbose', false, 'doplot', false);
end


function test_class_and_state(tc)
pm = tc.TestData.pm;
tc.verifyClass(pm, 'predictive_model');
tc.verifyTrue(pm.is_fitted, 'is_fitted should be true');
end


function test_canonical_path_fields_populated(tc)
pm = tc.TestData.pm;
tc.verifyNotEmpty(pm.Y);
tc.verifyNotEmpty(pm.fitted_values.yfit);
tc.verifyNotEmpty(pm.fitted_values.predictions);
tc.verifyNotEmpty(pm.fitted_values.Y_per_fold);
tc.verifyNotEmpty(pm.error_metrics.accuracy.value);
tc.verifyNotEmpty(pm.error_metrics.overallAccuracy.value);
tc.verifyNotEmpty(pm.cv_partition.trIdx);
tc.verifyNotEmpty(pm.cv_partition.teIdx);
tc.verifyNotEmpty(pm.cv_partition.nfolds);
end


function test_fold_models_match_nfolds(tc)
pm = tc.TestData.pm;
tc.verifyNumElements(pm.fold_models, pm.cv_partition.nfolds, ...
    'fold_models length should equal nfolds');
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
