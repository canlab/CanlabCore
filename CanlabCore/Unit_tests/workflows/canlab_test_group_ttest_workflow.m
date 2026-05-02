function tests = canlab_test_group_ttest_workflow
%TEST_GROUP_TTEST_WORKFLOW End-to-end smoke test of the canonical group-level pipeline.
%
% Validates that the chain
%     load_image_set -> ttest -> threshold -> region
% runs without error on the sample dataset and produces objects of the
% expected classes. This is an integration check, not a numerical-correctness
% check; per-step numerical invariants belong in the per-class test files.

tests = functiontests(localfunctions);
end


function test_full_workflow(tc)
imgs = load_image_set('emotionreg', 'noverbose');
tc.verifyClass(imgs, 'fmri_data');

t = ttest(imgs);
tc.verifyClass(t, 'statistic_image');

t = threshold(t, 0.005, 'unc');
tc.verifyClass(t, 'statistic_image');

r = region(t);
tc.verifyClass(r, 'region');
end
