function tests = test_replace_remove_empty
%TEST_REPLACE_REMOVE_EMPTY Round-trip and idempotence of empty-voxel state.
%
% These are the most important invariants in the toolbox: many subtle bugs
% have come from methods that assume the wrong state. We pin down:
%   1. replace_empty -> remove_empty -> replace_empty restores .dat row count
%   2. running the same transform twice is a no-op (idempotence).

tests = functiontests(localfunctions);
end


function test_round_trip_returns_to_padded_size(tc)
obj = get_sample_fmri_data();
obj_padded = replace_empty(obj);
n_padded = size(obj_padded.dat, 1);

obj_compressed = remove_empty(obj_padded);
obj_restored = replace_empty(obj_compressed);

tc.verifyEqual(size(obj_restored.dat, 1), n_padded, ...
    'replace_empty(remove_empty(replace_empty(obj))) should preserve voxel count');
end


function test_remove_empty_does_not_grow_dat(tc)
obj = get_sample_fmri_data();
obj_padded = replace_empty(obj);
obj_compressed = remove_empty(obj_padded);
tc.verifyLessThanOrEqual(size(obj_compressed.dat, 1), size(obj_padded.dat, 1));
end


function test_replace_empty_idempotent(tc)
obj = replace_empty(get_sample_fmri_data());
obj2 = replace_empty(obj);
tc.verifyEqual(size(obj.dat), size(obj2.dat));
end


function test_remove_empty_idempotent(tc)
obj = remove_empty(get_sample_fmri_data());
obj2 = remove_empty(obj);
tc.verifyEqual(size(obj.dat), size(obj2.dat));
end
