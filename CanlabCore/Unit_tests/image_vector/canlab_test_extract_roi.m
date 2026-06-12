function tests = canlab_test_extract_roi
%CANLAB_TEST_EXTRACT_ROI Data extraction from regions and masks.

tests = functiontests(localfunctions);
end


function test_extract_roi_averages_with_mask_file(tc)
% Extract one mean per image using the canonical brain mask. This is the
% simplest invocation of extract_roi_averages and exercises the typical
% (data, mask_filename) path that doesn't need on-the-fly resampling.
mask_file = which('brainmask_canlab.nii');
tc.assumeNotEmpty(mask_file, 'brainmask_canlab.nii not on path');

imgs = canlab_get_sample_fmri_data();
out = extract_roi_averages(imgs, mask_file);
tc.verifyTrue(isstruct(out) || isa(out, 'region'));
tc.verifyGreaterThan(numel(out), 0);

% from image object
mask = fmri_data(mask_file);
out = extract_roi_averages(imgs, mask_file);

tc.verifyTrue(isstruct(out) || isa(out, 'region'));
tc.verifyGreaterThan(numel(out), 0);

end


function test_extract_roi_averages_dat_shape_matches_images(tc)
% Each extracted region's .dat should have one row per image (30 for emotionreg).
mask_file = which('brainmask_canlab.nii');
tc.assumeNotEmpty(mask_file, 'brainmask_canlab.nii not on path');

imgs = canlab_get_sample_fmri_data();
out = extract_roi_averages(imgs, mask_file);
tc.assumeGreaterThan(numel(out), 0, 'no regions to check');
first = out(1);
tc.verifyTrue(isfield(first, 'dat') || isprop(first, 'dat'));
tc.verifyEqual(size(first.dat, 1), size(imgs.dat, 2));
end
