function tests = test_constructor
%TEST_CONSTRUCTOR fmri_data constructor invariants.
tests = functiontests(localfunctions);
end


function test_empty_constructor_returns_fmri_data(tc)
% fmri_data with no args produces an object of the right class with the
% standard volInfo populated from the default brain mask.
obj = fmri_data();
tc.verifyClass(obj, 'fmri_data');
tc.verifyNotEmpty(obj.volInfo, '.volInfo should be populated by the default mask');
end


function test_constructor_from_nifti(tc)
% Construct from a known 4-D NIfTI in Sample_datasets.
img = which('Wager_2008_emo_reg_vs_look_neg_contrast_images.nii.gz');
tc.assumeNotEmpty(img, 'Sample image not found on path');

obj = fmri_data(img, 'noverbose');
tc.verifyClass(obj, 'fmri_data');
tc.verifyEqual(size(obj.dat, 2), 30, ...
    'emotionreg sample is expected to contain 30 contrast images');
tc.verifyGreaterThan(size(obj.dat, 1), 0, ...
    '.dat should have at least one in-mask voxel');
end
