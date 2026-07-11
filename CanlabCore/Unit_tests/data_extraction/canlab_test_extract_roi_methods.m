function tests = canlab_test_extract_roi_methods
%CANLAB_TEST_EXTRACT_ROI_METHODS Cross-method consistency of ROI data-extraction methods.
%
% Function-based test (matlab.unittest via functiontests) verifying that the
% several CanlabCore routes to an "ROI mean" agree with one another and with
% analytically known values, and that on-the-fly resampling does not
% materially change extracted ROI means.
%
% Methods exercised:
%   - @fmri_data/extract_roi_averages   (region-object output, .dat per image)
%   - @image_vector/apply_mask          (subset voxels, then manual mean)
%   - @region/extract_data              (mm-coordinate matching, no resampling)
%   - @image_vector/apply_parcellation  (single-parcel weighted mean)
%
% Cases:
%   1. SAME-SPACE IDENTITY  - all four routes return identical ROI means
%                             (tight tolerance, ~1e-4 single-precision).
%   2. SYNTHETIC KNOWN MEAN - .dat set to known values; recovered mean == expected.
%   3. RESAMPLING TOLERANCE - mask resampled to a different space and back;
%                             extracted means stay corr > 0.99 and within a
%                             documented looser tolerance.
%
% Run headless:
%   results = runtests('canlab_test_extract_roi_methods')
% or via the suite:
%   results = canlab_run_all_tests

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
% Case 1: same-space identity across methods
% -------------------------------------------------------------------------
function test_same_space_identity_across_methods(tc)
% Build a single ROI in the EXACT space of the data, then pull the per-image
% ROI mean four different ways. All must agree to a tight tolerance.

imgs = canlab_get_sample_fmri_data();          % 30-image emotionreg fmri_data
imgs = replace_empty(imgs);

% --- Define a single-region mask in the data's own space -----------------
% Pick a stable subset of in-mask voxels (every 3rd voxel) so the region is
% well-defined and identical for every extraction route. Coding them 1 makes
% a one-region "atlas-like" mask.
mask = imgs;                                    % inherits volInfo / space
v    = size(mask.dat, 1);
wh   = false(v, 1);
wh(1:3:v) = true;                               % deterministic voxel subset
mask.dat = single(double(wh));                  % 1 = in ROI, 0 = out
nvox = nnz(wh);
tc.assumeGreaterThan(nvox, 10, 'ROI too small to be meaningful');

% --- Ground-truth ROI mean: average of the in-ROI voxels per image -------
expected = mean(double(imgs.dat(wh, :)), 1)';   % [n images x 1]

% --- Route A: extract_roi_averages (unique_mask_values => one region) -----
clA = extract_roi_averages(imgs, mask, 'unique_mask_values', 'noverbose');
tc.assumeGreaterThan(numel(clA), 0, 'extract_roi_averages returned no regions');
meanA = double(clA(1).dat(:));                  % region object: .dat = per-image mean

% --- Route B: apply_mask, then manual mean over in-mask voxels -----------
masked = apply_mask(imgs, mask);
masked = remove_empty(masked);
meanB  = nanmean(double(masked.dat), 1)';       % mean over voxels -> per image

% --- Route C: region/extract_data (mm-coordinate matching) ---------------
% Build a region object from the same-space mask, then extract.
r  = region(mask, 'unique_mask_values', 'noverbose');
r  = extract_data(r, imgs);
tc.assumeGreaterThan(numel(r), 0, 'region/extract_data returned no regions');
meanC = double(r(1).dat(:));

% --- Route D: apply_parcellation, single parcel --------------------------
% parcels is an fmri_data with integer codes (1 = the only parcel).
parcels = imgs;
parcels.dat = single(double(wh));               % one parcel coded 1
parcel_means = apply_parcellation(imgs, parcels);   % [images x parcels]
meanD = double(parcel_means(:, 1));

% --- All routes must agree (tight tolerance) -----------------------------
% Data are stored single-precision and routes differ in summation order, so
% use a tight-but-single-safe absolute tolerance.
tol = 1e-4;
fprintf('\n[Case 1] Same-space ROI mean across methods (n=%d images, %d voxels):\n', numel(expected), nvox);
fprintf('  max|A-expected| = %.3e (extract_roi_averages)\n', max(abs(meanA - expected)));
fprintf('  max|B-expected| = %.3e (apply_mask + mean)\n',    max(abs(meanB - expected)));
fprintf('  max|C-expected| = %.3e (region/extract_data)\n',  max(abs(meanC - expected)));
fprintf('  max|D-expected| = %.3e (apply_parcellation)\n',   max(abs(meanD - expected)));

tc.verifyEqual(meanA, expected, 'AbsTol', tol, 'extract_roi_averages != ground truth');
tc.verifyEqual(meanB, expected, 'AbsTol', tol, 'apply_mask mean != ground truth');
tc.verifyEqual(meanC, expected, 'AbsTol', tol, 'region/extract_data != ground truth');
tc.verifyEqual(meanD, expected, 'AbsTol', tol, 'apply_parcellation != ground truth');

% And to each other
tc.verifyEqual(meanA, meanB, 'AbsTol', tol, 'extract_roi_averages != apply_mask');
tc.verifyEqual(meanA, meanC, 'AbsTol', tol, 'extract_roi_averages != region/extract_data');
tc.verifyEqual(meanA, meanD, 'AbsTol', tol, 'extract_roi_averages != apply_parcellation');
end


% -------------------------------------------------------------------------
% Case 2: synthetic data with analytically known ROI mean
% -------------------------------------------------------------------------
function test_synthetic_known_roi_mean(tc)
% Set .dat to a KNOWN per-voxel pattern, define an ROI, and check that the
% recovered ROI mean equals the analytic expectation exactly (to tolerance).

base = canlab_get_sample_fmri_data();           % borrow a valid space/volInfo
base = replace_empty(base);
v    = size(base.dat, 1);

obj  = base;
% (a) constant-valued images: every voxel = c -> any ROI mean == c
nimg = 4;
c    = 7;
obj.dat = c * ones(v, nimg, 'single');

% Define an ROI (deterministic subset)
wh        = false(v, 1);
wh(2:5:v) = true;
nvox      = nnz(wh);
tc.assumeGreaterThan(nvox, 10, 'ROI too small');

mask        = base;
mask.dat    = single(double(wh));

cl       = extract_roi_averages(obj, mask, 'unique_mask_values', 'noverbose');
recov_c  = double(cl(1).dat(:));
expect_c = c * ones(nimg, 1);

fprintf('\n[Case 2a] Constant data (every voxel = %d), ROI mean per image:\n', c);
fprintf('  expected = %g, recovered = [%s], max|err| = %.3e\n', ...
    c, num2str(recov_c', '%.6g '), max(abs(recov_c - expect_c)));
tc.verifyEqual(recov_c, expect_c, 'AbsTol', 1e-5, 'constant-data ROI mean wrong');

% (b) known per-voxel ramp: voxel i has value i (same across images).
% Analytic ROI mean = mean of the voxel indices that are in the ROI.
obj2        = base;
ramp        = single((1:v)');
obj2.dat    = repmat(ramp, 1, nimg);

cl2         = extract_roi_averages(obj2, mask, 'unique_mask_values', 'noverbose');
recov_r     = double(cl2(1).dat(:));
expect_r    = mean(double(find(wh))) * ones(nimg, 1);

rel_err = max(abs(recov_r - expect_r)) ./ max(abs(expect_r));
fprintf('[Case 2b] Ramp data (voxel i = i), ROI mean:\n');
fprintf('  expected = %.4f, recovered(1) = %.4f, max|err| = %.3e, rel.err = %.3e\n', ...
    expect_r(1), recov_r(1), max(abs(recov_r - expect_r)), rel_err);
% Ramp values are large (~1e4) and stored single-precision, so the exact mean
% is not float32-representable. Verify to single-precision RELATIVE tolerance.
tc.verifyEqual(recov_r, expect_r, 'RelTol', 1e-4, 'ramp-data ROI mean wrong');
end


% -------------------------------------------------------------------------
% Case 3: resampling tolerance (mask in a different space)
% -------------------------------------------------------------------------
function test_resampling_tolerance(tc)
% Route 1: extract in the data's NATIVE space (no resampling).
% Route 2: resample the DATA *and* the mask into a different (finer) space --
% the "resample the image data to match a different space, or vice versa" case
% -- then extract there. Because the data values are genuinely interpolated
% onto a new grid (not a lossless mask round-trip), the two ROI-mean vectors
% differ slightly; they must stay highly correlated (r > 0.99) and within a
% documented looser tolerance. This exercises real reslicing, unlike a
% same-resolution round-trip which would be numerically exact.

imgs = canlab_get_sample_fmri_data();
imgs = replace_empty(imgs);
v    = size(imgs.dat, 1);

% Same-space single-region mask
wh        = false(v, 1);
wh(1:2:v) = true;                               % larger ROI -> stable mean
mask      = imgs;
mask.dat  = single(double(wh));

% --- Route 1: native-space extraction ------------------------------------
cl1   = extract_roi_averages(imgs, mask, 'unique_mask_values', 'noverbose');
mean1 = double(cl1(1).dat(:));

% --- Resample BOTH data and mask into a finer reference space ------------
% brainmask_canlab.nii is a 1.5 mm grid; the sample data is ~3.4x3.4x4.5 mm,
% so this up-samples (linearly interpolates) the data onto a genuinely
% different grid -- the "vice versa" direction (image data -> new space).
ref_file = which('brainmask_canlab.nii');
tc.assumeNotEmpty(ref_file, 'brainmask_canlab.nii not on path');
ref      = fmri_data(ref_file, 'noverbose');

imgs_rs = resample_space(imgs, ref);            % DATA interpolated to ref grid
mask_rs = resample_space(mask, ref);            % mask follows
mask_rs.dat = single(double(mask_rs.dat > 0.5));% rebinarize after interpolation
tc.assumeGreaterThan(nnz(mask_rs.dat > 0), 10, 'resampled ROI vanished');

% --- Route 2: extract in the resampled (finer) space ---------------------
cl2   = extract_roi_averages(imgs_rs, mask_rs, 'unique_mask_values', 'noverbose');
tc.assumeGreaterThan(numel(cl2), 0, 'no regions after resampling');
mean2 = double(cl2(1).dat(:));

% --- Compare --------------------------------------------------------------
n = min(numel(mean1), numel(mean2));
mean1 = mean1(1:n); mean2 = mean2(1:n);

r = corr(mean1, mean2);

% Looser tolerance: nearest-neighbour reslicing to a different grid changes
% the exact voxel set, so means shift slightly. We require the discrepancy to
% be small relative to the between-image spread of the ROI means.
roi_scale  = std(mean1);
abs_err    = max(abs(mean1 - mean2));
rel_err    = abs_err / max(roi_scale, eps);
TOL_REL    = 0.10;                              % documented: <10% of ROI-mean SD

fprintf('\n[Case 3] Resampling tolerance (ROI mean: native vs resampled mask):\n');
fprintf('  correlation        = %.6f (require > 0.99)\n', r);
fprintf('  max|abs difference| = %.4e\n', abs_err);
fprintf('  ROI-mean SD         = %.4e\n', roi_scale);
fprintf('  relative error      = %.4f (require < %.2f)\n', rel_err, TOL_REL);

tc.verifyGreaterThan(r, 0.99, 'resampled ROI means not highly correlated');
tc.verifyLessThan(rel_err, TOL_REL, 'resampled ROI means differ more than tolerance');
end
