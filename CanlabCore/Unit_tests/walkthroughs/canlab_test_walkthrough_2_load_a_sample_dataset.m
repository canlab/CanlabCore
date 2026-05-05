function tests = canlab_test_walkthrough_2_load_a_sample_dataset
%CANLAB_TEST_WALKTHROUGH_2_LOAD_A_SAMPLE_DATASET Smoke test of canlab_help_2_load_a_sample_dataset.
%
% Tier B integration test: runs the walkthrough end-to-end and counts a
% return without error as a pass. See canlab_test_walkthrough_1_installing_tools
% for shared rationale.

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
script_name = 'canlab_help_2_load_a_sample_dataset';
tc.assumeNotEmpty(which(script_name), ...
    'CANlab_help_examples (example_help_files/) must be on the MATLAB path');
evalc(script_name);
tc.verifyTrue(true);
end
