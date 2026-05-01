function results = run_all_tests(varargin)
%RUN_ALL_TESTS Run the CanlabCore unit test suite.
%
%   results = run_all_tests
%   results = run_all_tests('JUnit', 'test-results.xml')
%   results = run_all_tests('Tag', 'Core')              % run only tagged tests
%   results = run_all_tests('Tag', '~RequiresMasks')    % skip a tag (prefix with ~)
%
%   Discovers tests in this folder and its subfolders, excluding
%   old_to_integrate/. Adds the parent CanlabCore tree to the path before
%   running. Returns a matlab.unittest.TestResult array. Suitable for both
%   interactive use and CI; in CI, throw on failure with assertSuccess(results).

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

suite = TestSuite.fromFolder(this_dir, 'IncludingSubfolders', true);

% Drop anything under old_to_integrate/ — those are legacy scripts pending
% rewrite, not matlab.unittest tests.
if ~isempty(suite)
    base_folders = {suite.BaseFolder};
    keep = ~contains(base_folders, [filesep 'old_to_integrate']);
    suite = suite(keep);
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
