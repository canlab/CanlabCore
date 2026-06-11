function tests = canlab_test_walkthrough_7_multivariate_prediction
%CANLAB_TEST_WALKTHROUGH_7_MULTIVARIATE_PREDICTION Smoke test of canlab_help_7_multivariate_prediction_basics.
%
% Tier B integration test. See canlab_test_walkthrough_1_installing_tools
% for shared rationale and conventions.

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
script_name = 'canlab_help_7_multivariate_prediction_basics';
tc.assumeNotEmpty(which(script_name), ...
    'CANlab_help_examples (example_help_files/) must be on the MATLAB path');
evalc(script_name);
tc.verifyTrue(true);
end
