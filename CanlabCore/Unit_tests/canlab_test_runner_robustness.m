function tests = canlab_test_runner_robustness
%CANLAB_TEST_RUNNER_ROBUSTNESS Discovery must not abort on a stray non-test file.
%
% canlab_run_all_tests globs every canlab_test_*.m under Unit_tests and
% concatenates the suites. A file that matches the name pattern but is not a
% valid matlab.unittest test (an old plain-assert script, or a function whose
% internal name does not match the filename) makes TestSuite.fromFile throw
% MATLAB:unittest:TestSuite:NonTestFile. Without guarding, that single error
% aborts discovery of the entire suite. canlab_safe_suite_from_file warn-skips
% such files instead; these tests pin that behavior.

tests = functiontests(localfunctions);
end


function test_nontest_file_is_warn_skipped(tc)   %#ok<*DEFNU>
% A canlab_test_*.m whose function name != filename and that never calls
% functiontests is not a valid test file: skip it, warn, return empty.
folder = tc.applyFixture( ...
    matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
fpath = local_write(folder, 'canlab_test_bogus_nontest.m', { ...
    'function some_other_name()'
    'assert(true);'
    'end'});

% Use lastwarn rather than verifyWarning so the check does not depend on the
% ambient warning display state (a caller may have warnings off; lastwarn
% still records the id even when the warning event is not displayed).
lastwarn('', '');
[suite, ok] = canlab_safe_suite_from_file(fpath);
[~, warn_id] = lastwarn;

tc.verifyFalse(ok, 'expected ok=false for a non-test file');
tc.verifyEmpty(suite, 'expected an empty suite for a non-test file');
tc.verifyEqual(warn_id, 'canlab_run_all_tests:skippedNonTestFile', ...
    'expected a skippedNonTestFile warning');
end


function test_valid_test_file_is_loaded(tc)
% A well-formed functiontests file loads normally (ok=true, non-empty suite).
folder = tc.applyFixture( ...
    matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
fpath = local_write(folder, 'canlab_test_valid_dummy.m', { ...
    'function tests = canlab_test_valid_dummy'
    'tests = functiontests(localfunctions);'
    'end'
    ''
    'function test_trivial(tc)'
    'tc.verifyTrue(true);'
    'end'});

[suite, ok] = canlab_safe_suite_from_file(fpath);
tc.verifyTrue(ok, 'expected ok=true for a valid test file');
tc.verifyNotEmpty(suite);
tc.verifyEqual(numel(suite), 1);
end


function fpath = local_write(folder, name, lines)
% Write a cellstr of lines to folder/name and return the full path. Uses
% fopen/fprintf (not writelines) so the test runs on older MATLAB too.
fpath = fullfile(folder, name);
fid = fopen(fpath, 'w');
tc_assert_open(fid, fpath);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', lines{:});
end


function tc_assert_open(fid, fpath)
if fid < 0
    error('canlab_test_runner_robustness:cannotWrite', ...
        'Could not open %s for writing.', fpath);
end
end
