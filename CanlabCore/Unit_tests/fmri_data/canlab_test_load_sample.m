function tests = canlab_test_load_sample
%TEST_LOAD_SAMPLE load_image_set keyword resolution for the sample dataset.
tests = functiontests(localfunctions);
end


function test_load_image_set_emotionreg(tc)
obj = load_image_set('emotionreg', 'noverbose');
tc.verifyClass(obj, 'fmri_data');
tc.verifyEqual(size(obj.dat, 2), 30);
end


function test_load_image_set_attaches_paths(tc)
% After load, fullpath / image_names should be populated for provenance.
obj = load_image_set('emotionreg', 'noverbose');
tc.verifyNotEmpty(obj.fullpath, '.fullpath should record the source image(s)');
end
