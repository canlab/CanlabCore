function tests = canlab_test_similarity
%CANLAB_TEST_SIMILARITY Spatial-similarity / annotation methods.

tests = functiontests(localfunctions);
end


function test_jackknife_similarity_runs(tc)
% jackknife_similarity computes leave-one-out similarity across images.
% Uses canonical 'doplot'/'verbose' flags (logical) per the function's
% inputParser, not a bare 'noplot' string.
obj = canlab_get_sample_fmri_data();
result = jackknife_similarity(obj, 'doplot', false, 'verbose', false);
tc.verifyTrue(isnumeric(result) || isstruct(result));
end


function test_image_similarity_plot_runs(tc)
% image_similarity_plot needs a basis set; the bucknermaps variant pulls
% one in by name. Skip cleanly if the maps repo isn't on path.
tc.assumeNotEmpty(which('Buckner_Yeo_networks_7.nii.gz') | ...
                  which('Yeo_7networks.nii') | ...
                  which('rBucknerlab_wholebrain.nii'), ...
    'Buckner network maps not on path');
obj = canlab_get_sample_fmri_data();
m = mean(obj);
fn = @() image_similarity_plot_bucknermaps(m, 'noplot');
tc.verifyWarningFree(fn);
end
