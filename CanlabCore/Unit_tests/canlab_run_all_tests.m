function results = canlab_run_all_tests(varargin)
%CANLAB_RUN_ALL_TESTS Run the CanlabCore unit test suite.
%
%   results = canlab_run_all_tests
%   results = canlab_run_all_tests('JUnit', 'test-results.xml')
%   results = canlab_run_all_tests('Tag', 'Core')              % run only tagged tests
%   results = canlab_run_all_tests('Tag', '~RequiresMasks')    % skip a tag
%
%   Discovers test files matching canlab_test_*.m anywhere under this
%   folder (excluding old_to_integrate/) and runs them via matlab.unittest.
%   Adds the parent CanlabCore tree to the path before running. Returns a
%   matlab.unittest.TestResult array suitable for both interactive use and
%   CI; in CI, throw on failure with assertSuccess(results).
%
%   Note on filename pattern: test files are namespaced as canlab_test_*
%   to avoid colliding with other toolboxes a user might have on their
%   path. matlab.unittest's TestSuite.fromFolder only auto-discovers
%   files whose names start or end with "test", so this runner does its
%   own glob and uses TestSuite.fromFile on each match.

p = inputParser;
p.addParameter('JUnit', '', @(x) ischar(x) || isstring(x));
p.addParameter('Tag', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.XMLPlugin
import matlab.unittest.selectors.HasTag
import matlab.unittest.Verbosity

this_dir = fileparts(mfilename('fullpath'));
canlabcore_root = fileparts(this_dir);
addpath(genpath(canlabcore_root));

% Glob for canlab_test_*.m files; skip old_to_integrate/.
% Note: TestSuite.fromFile returns matlab.unittest.Test arrays (Test is a
% concrete subclass of TestSuite), so we initialize with Test.empty rather
% than TestSuite.empty — mixing the abstract empty with concrete Tests
% errors during horzcat.
matches = dir(fullfile(this_dir, '**', 'canlab_test_*.m'));
suite = matlab.unittest.Test.empty;
for k = 1:numel(matches)
    fpath = fullfile(matches(k).folder, matches(k).name);
    if contains(fpath, [filesep 'old_to_integrate' filesep])
        continue
    end
    suite = [suite, TestSuite.fromFile(fpath)]; %#ok<AGROW>
end

tag = char(p.Results.Tag);
if ~isempty(tag)
    if startsWith(tag, '~')
        suite = suite.selectIf(~HasTag(tag(2:end)));
    else
        suite = suite.selectIf(HasTag(tag));
    end
end

runner = TestRunner.withTextOutput('OutputDetail', Verbosity.Concise);

junit_path = char(p.Results.JUnit);
if ~isempty(junit_path)
    junit_dir = fileparts(junit_path);
    if ~isempty(junit_dir) && ~exist(junit_dir, 'dir')
        mkdir(junit_dir);
    end
    runner.addPlugin(XMLPlugin.producingJUnitFormat(junit_path));
end

results = runner.run(suite);

n_passed = sum([results.Passed]);
n_failed = sum([results.Failed]);
n_incomplete = sum([results.Incomplete]);
fprintf('\n=== %d passed, %d failed, %d incomplete ===\n', ...
    n_passed, n_failed, n_incomplete);

if nargout == 0
    clear results
end
end
