function tests = canlab_test_resample_space
%CANLAB_TEST_RESAMPLE_SPACE Resampling and space-comparison utilities.

tests = functiontests(localfunctions);
end


function test_compare_space_self_is_zero(tc)
% Comparing an object to itself should report identical space.
% compare_space documents return 0 for same-space; in practice it returns
% the result of any(...), which is logical false. We accept either.
obj = canlab_get_sample_fmri_data();
result = compare_space(obj, obj);
tc.verifyFalse(logical(result), 'compare_space(obj, obj) should indicate same space');
end


function test_resample_space_to_self_preserves_image_count(tc)
% Resampling onto the same space should be a no-op for image count.
obj = canlab_get_sample_fmri_data();
n_imgs = size(obj.dat, 2);
resampled = resample_space(obj, obj);
tc.verifyEqual(size(resampled.dat, 2), n_imgs);
tc.verifyClass(resampled, class(obj));
end
