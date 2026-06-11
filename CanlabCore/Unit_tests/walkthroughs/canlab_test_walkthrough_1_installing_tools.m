function tests = canlab_test_walkthrough_1_installing_tools
%CANLAB_TEST_WALKTHROUGH_1_INSTALLING_TOOLS End-to-end smoke test of canlab_help_1_installing_tools.
%
% Wrapper around CANlab_help_examples/example_help_files/canlab_help_1_installing_tools.m.
% Acts as a Tier B integration test: confirms the walkthrough runs end-to-end
% on the current toolbox without erroring. Does not check pixel content,
% statistics, or output values; reaching the end of the script counts as a pass.
%
% Skipped (Incomplete) if CANlab_help_examples is not on the MATLAB path.

tests = functiontests(localfunctions);
end


function setup(tc)             %#ok<*DEFNU>
% Run each walkthrough in its own temp working directory so output files
% (NIfTI writes, saved figures) don't collide with prior runs.
tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
close all force
end


function teardown(tc)
close all force
end


function test_runs_to_completion(tc)
script_name = 'canlab_help_1_installing_tools';
tc.assumeNotEmpty(which(script_name), ...
    'CANlab_help_examples (example_help_files/) must be on the MATLAB path');
evalc(script_name);
tc.verifyTrue(true);
end
