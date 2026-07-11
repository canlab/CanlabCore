function tests = canlab_test_walkthrough_2b_basic_image_visualization
%CANLAB_TEST_WALKTHROUGH_2B_BASIC_IMAGE_VISUALIZATION Hardened smoke test: Basic image visualization walkthrough.
%
% Tier B integration test. Runs a frozen snapshot of the CANlab_help_examples
% walkthrough "canlab_help_2b_basic_image_visualization" through the headless, fault-tolerant
% harness canlab_run_walkthrough_snapshot. The snapshot lives under
% Unit_tests/walkthroughs/private/ (a verbatim copy; refresh it from
% example_help_files/ when the tutorial changes). Reaching the end of the
% script without a genuine error counts as a pass; sections that need a
% display, interactive input, or optional data files are skipped, not failed.
%
% These are SLOW and are skipped by the per-push suite. Run them with
%   canlab_run_all_tests('Walkthroughs', 'only')   % nightly tier
% See canlab_run_walkthrough_snapshot for the hardening details.

tests = functiontests(localfunctions);
end


function setupOnce(tc)          %#ok<*DEFNU>
here = fileparts(mfilename('fullpath'));
tc.TestData.snapshot = fullfile(here, 'private', 'canlab_help_2b_basic_image_visualization.m');
end


function setup(tc)
% Clean cwd per run (snapshots write NIfTI/figure files) and headless graphics.
tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
tc.TestData.PrevFigVis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
close all force
end


function teardown(tc)
close all force
if isfield(tc.TestData, 'PrevFigVis')
    set(0, 'DefaultFigureVisible', tc.TestData.PrevFigVis);
end
end


function test_runs_to_completion(tc)
canlab_run_walkthrough_snapshot(tc, tc.TestData.snapshot);
end
