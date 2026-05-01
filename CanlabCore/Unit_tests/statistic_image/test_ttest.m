function tests = test_ttest
%TEST_TTEST Voxelwise one-sample t-test on an fmri_data object.
tests = functiontests(localfunctions);
end


function test_ttest_returns_statistic_image(tc)
obj = get_sample_fmri_data();
t = ttest(obj);
tc.verifyClass(t, 'statistic_image');
tc.verifyGreaterThan(size(t.dat, 1), 0);
end


function test_threshold_returns_statistic_image(tc)
obj = get_sample_fmri_data();
t = ttest(obj);
t_thr = threshold(t, 0.005, 'unc');
tc.verifyClass(t_thr, 'statistic_image');
end


function test_threshold_marks_significant_voxels(tc)
% After thresholding, some voxels should be flagged in .sig (or none, if
% the threshold is too strict). We just verify the field exists and is
% the right shape.
obj = get_sample_fmri_data();
t = ttest(obj);
t_thr = threshold(t, 0.05, 'unc');
tc.verifyTrue(isprop(t_thr, 'sig'));
tc.verifyEqual(size(t_thr.sig, 1), size(t_thr.dat, 1));
end
