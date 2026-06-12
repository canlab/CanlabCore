function tests = canlab_test_cat_split
%CANLAB_TEST_CAT_SPLIT Basic image math: cat, split, mean, get_wh_image.
%
% Uses the emotionreg sample (30 images) as a fixture and exercises the
% combine/subset operators without going near full-volume computation.

tests = functiontests(localfunctions);
end


function test_cat_doubles_image_count(tc)
% cat([a; a]) along the image axis should give 2x the image count.
obj = canlab_get_sample_fmri_data();
n = size(obj.dat, 2);
combined = cat(obj, obj);
tc.verifyClass(combined, 'fmri_data');
tc.verifyEqual(size(combined.dat, 2), 2 * n);
end


function test_split_then_cat_round_trip(tc)
% split() splits according to obj.images_per_session, so we set that first
% to two halves and verify cat(split) preserves the original image count.
obj = canlab_get_sample_fmri_data();
n = size(obj.dat, 2);
half = floor(n / 2);
obj.images_per_session = [half, n - half];

halves = split(obj);
tc.verifyEqual(numel(halves), 2);

% halves is a cell array of fmri_data objects.
rebuilt = cat(halves{1}, halves{2});
tc.verifyEqual(size(rebuilt.dat, 2), n);
end


function test_get_wh_image_returns_subset(tc)
% Pull the first 5 images; .dat should narrow accordingly, voxels unchanged.
obj = canlab_get_sample_fmri_data();
sub = get_wh_image(obj, 1:5);
tc.verifyClass(sub, 'fmri_data');
tc.verifyEqual(size(sub.dat, 2), 5);
tc.verifyEqual(size(sub.dat, 1), size(obj.dat, 1));
end


function test_mean_collapses_to_single_image(tc)
% mean across images should produce one image with the same voxel count.
obj = canlab_get_sample_fmri_data();
m = mean(obj);
tc.verifyEqual(size(m.dat, 2), 1);
tc.verifyEqual(size(m.dat, 1), size(obj.dat, 1));
end
