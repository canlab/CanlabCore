function tests = canlab_test_walkthrough_4_masking_and_writing
%CANLAB_TEST_WALKTHROUGH_4_MASKING_AND_WRITING Smoke test of canlab_help_4_masking_and_writing_nifti_files.
%
% Tier B integration test. See canlab_test_walkthrough_1_installing_tools
% for shared rationale and conventions. The WorkingFolderFixture in setup
% provides a clean cwd, which matters here because the walkthrough writes
% NIfTI output files.

tests = functiontests(localfunctions);
end


function setup(tc)             %#ok<*DEFNU>
tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
close all force
end


function teardown(tc)
close all force
end


function test_runs_to_completion(tc)
script_name = 'canlab_help_4_masking_and_writing_nifti_files';
tc.assumeNotEmpty(which(script_name), ...
    'CANlab_help_examples (example_help_files/) must be on the MATLAB path');
evalc(script_name);
tc.verifyTrue(true);
end
