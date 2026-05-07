function tests = canlab_test_misc
%CANLAB_TEST_MISC Misc utilities: flip, enforce_variable_types, history.

tests = functiontests(localfunctions);
end


function test_flip_preserves_class_and_image_count(tc)
% L-R flip should not change object class or image count. flip() requires
% a single-image (3-D) object, so reduce emotionreg to its mean first.
obj = mean(canlab_get_sample_fmri_data());
flipped = flip(obj);
tc.verifyClass(flipped, class(obj));
tc.verifyEqual(size(flipped.dat, 2), size(obj.dat, 2));
end


function test_flip_twice_round_trips(tc)
% Flipping twice should return an object with the same shape as the input.
obj = mean(canlab_get_sample_fmri_data());
flipped_twice = flip(flip(obj));
tc.verifyEqual(size(flipped_twice.dat), size(obj.dat));
end


function test_flip_errors_on_multi_image_object(tc)
% flip() must reject objects with more than one image rather than
% silently producing nonsense via the 3-D slice loop.
obj = canlab_get_sample_fmri_data();   % 30 images
tc.verifyError(@() flip(obj), 'image_vector:flip:multiImage');
end


function test_enforce_variable_types_returns_same_class(tc)
% enforce_variable_types should not change the object's class.
obj = canlab_get_sample_fmri_data();
casted = enforce_variable_types(obj);
tc.verifyClass(casted, class(obj));
tc.verifyEqual(size(casted.dat), size(obj.dat));
end
