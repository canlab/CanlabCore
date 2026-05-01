function tests = test_region_from_stat
%TEST_REGION_FROM_STAT Building a region object from a thresholded map.
tests = functiontests(localfunctions);
end


function test_region_returns_region_object(tc)
% Use a permissive threshold so we get at least one cluster on the sample.
obj = get_sample_fmri_data();
t = ttest(obj);
t = threshold(t, 0.05, 'unc');
r = region(t);
tc.verifyClass(r, 'region');
end


function test_region_count_nonnegative(tc)
obj = get_sample_fmri_data();
t = ttest(obj);
t = threshold(t, 0.05, 'unc');
r = region(t);
tc.verifyGreaterThanOrEqual(numel(r), 0);
end
