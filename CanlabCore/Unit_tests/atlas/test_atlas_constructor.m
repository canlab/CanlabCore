function tests = test_atlas_constructor
%TEST_ATLAS_CONSTRUCTOR Smoke test for atlas() with no external files.
tests = functiontests(localfunctions);
end


function test_empty_atlas_constructor(tc)
atl = atlas();
tc.verifyClass(atl, 'atlas');
end
