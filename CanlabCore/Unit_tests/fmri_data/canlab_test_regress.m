function tests = canlab_test_regress
%CANLAB_TEST_REGRESS Voxelwise regression smoke tests.

tests = functiontests(localfunctions);
end


function test_regress_with_intercept_only(tc)
% With a constant predictor, regress should run and produce stat output.
obj = canlab_get_sample_fmri_data();
n = size(obj.dat, 2);
obj.X = ones(n, 1);
out = regress(obj, 0.05, 'unc', 'noverbose', 'nodisplay');
tc.verifyTrue(isstruct(out) || isobject(out));
end


function test_regress_with_one_random_predictor(tc)
% Set a random regressor and verify regress runs end-to-end.
rng(0);
obj = canlab_get_sample_fmri_data();
n = size(obj.dat, 2);
obj.X = randn(n, 1);
out = regress(obj, 0.05, 'unc', 'noverbose', 'nodisplay');
tc.verifyTrue(isstruct(out) || isobject(out));
end
