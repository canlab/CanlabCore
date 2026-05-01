function tests = test_apply_mask
%TEST_APPLY_MASK Basic apply_mask behavior.
tests = functiontests(localfunctions);
end


function test_apply_mask_preserves_image_count(tc)
% apply_mask changes voxels (rows) but never images (cols).
obj = get_sample_fmri_data();
mask_file = which('brainmask_canlab.nii');
tc.assumeNotEmpty(mask_file, 'brainmask_canlab.nii not on path');

obj_masked = apply_mask(obj, mask_file);
tc.verifyClass(obj_masked, 'fmri_data');
tc.verifyEqual(size(obj_masked.dat, 2), size(obj.dat, 2));
end


function test_apply_mask_returns_same_class(tc)
% apply_mask should preserve object class (fmri_data in -> fmri_data out).
obj = get_sample_fmri_data();
mask_file = which('brainmask_canlab.nii');
tc.assumeNotEmpty(mask_file, 'brainmask_canlab.nii not on path');

obj_masked = apply_mask(obj, mask_file);
tc.verifyEqual(class(obj_masked), class(obj));
end
