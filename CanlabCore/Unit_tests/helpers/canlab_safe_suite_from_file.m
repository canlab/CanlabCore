function [suite, ok] = canlab_safe_suite_from_file(fpath)
%CANLAB_SAFE_SUITE_FROM_FILE Build a test suite from one file without aborting on bad files.
%
%   [suite, ok] = canlab_safe_suite_from_file(fpath)
%
%   Wraps matlab.unittest.TestSuite.fromFile so that a file which is not a
%   valid test (e.g. a stray canlab_test_*.m that is an old plain-assert
%   script or whose internal function name does not match the filename)
%   does not throw and abort discovery of the whole suite. Such a file is
%   warn-skipped and an empty Test array is returned instead.
%
%   This matters because canlab_run_all_tests globs every canlab_test_*.m
%   under Unit_tests and concatenates the results; a single NonTestFile
%   error from fromFile would otherwise take down the entire run.
%
% :Inputs:
%   **fpath:**  absolute path to a candidate .m test file.
%
% :Outputs:
%   **suite:**  a matlab.unittest.Test array (empty if the file is not a
%               valid test file).
%   **ok:**     logical, true if fromFile succeeded, false if the file was
%               skipped.
%
% :See also: canlab_run_all_tests, matlab.unittest.TestSuite

ok = true;
suite = matlab.unittest.Test.empty;

try
    suite = matlab.unittest.TestSuite.fromFile(fpath);
catch ME
    ok = false;
    if strcmp(ME.identifier, 'MATLAB:unittest:TestSuite:NonTestFile')
        warning('canlab_run_all_tests:skippedNonTestFile', ...
            'Skipping %s: not a valid matlab.unittest test file.', fpath);
    else
        warning('canlab_run_all_tests:fromFileError', ...
            'Skipping %s: %s (%s)', fpath, ME.message, ME.identifier);
    end
end
end
