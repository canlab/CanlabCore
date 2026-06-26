function tests = canlab_test_check_roi_extraction
%CANLAB_TEST_CHECK_ROI_EXTRACTION Multi-region ROI extraction + write/reload round-trip.
%
% Complements canlab_test_extract_roi (single-mask extraction) by exercising
% the two things that test does not:
%   1. Multi-region extraction with 'unique_mask_values' (one average per
%      integer label in a parcellation), using atlas_labels_combined.img.
%   2. Invariance of the extracted region averages across a write-to-disk and
%      reload cycle - a regression guard on the NIfTI I/O path.
%
% Converted from the old standalone script Unit_tests/old_to_integrate/
% check_roi_extraction.m, which printed PASS/FAIL but never asserted.

tests = functiontests(localfunctions);
end


function test_unique_mask_values_roundtrip(tc)   %#ok<*DEFNU>
mask_file = which('atlas_labels_combined.img');
tc.assumeNotEmpty(mask_file, 'atlas_labels_combined.img not on path');

mask_image = fmri_data(mask_file, 'noverbose');

% Synthetic per-region timeseries: region i carries the signal ts*sqrt(i),
% so every region has a distinct, known mean trajectory across images.
wh_region = mask_image.dat;
regions   = unique(wh_region);
nimgs     = 20;
ts        = linspace(-10, 10, nimgs);

dat = mask_image;
dat.dat = zeros(size(mask_image.dat, 1), nimgs);
for i = 1:numel(regions)
    idx = wh_region == regions(i);
    dat.dat(idx, :) = repmat(ts .* sqrt(i), sum(idx), 1);
end

cl1      = extract_roi_averages(dat, mask_image, 'unique_mask_values');
all_reg1 = cat(2, cl1(:).dat);

% One row per image, one column per non-zero region label.
tc.verifyEqual(size(all_reg1, 1), nimgs);
tc.verifyGreaterThan(size(all_reg1, 2), 1);

% Round-trip through disk in a scratch folder, then re-extract.
tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
dat.fullpath = fullfile(pwd, 'test_roi_image.nii');
write(dat, 'overwrite');
reloaded = fmri_data(dat.fullpath, 'noverbose');

cl2      = extract_roi_averages(reloaded, mask_image, 'unique_mask_values');
all_reg2 = cat(2, cl2(:).dat);

tc.verifyEqual(size(all_reg2), size(all_reg1));
% Small absolute slack absorbs single-precision NIfTI write rounding (the
% original script used a 1e-3 relative threshold for the same reason).
tc.verifyEqual(all_reg2, all_reg1, 'AbsTol', 1e-2, 'RelTol', 1e-3);
end
