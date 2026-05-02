# CanlabCore unit tests

This is the unit-test suite for CanlabCore. It uses MATLAB's built-in
`matlab.unittest` framework. There is no external runner — you can run any
single test by `runtests('canlab_test_constructor')`, or the whole suite
with `canlab_run_all_tests`.

Test files (and the runner and shared helpers) are namespaced with a
`canlab_` prefix to avoid colliding with similarly-named files from other
toolboxes a user may have on their MATLAB path.

## Running locally

```matlab
cd CanlabCore/Unit_tests
results = canlab_run_all_tests;                          % full suite
results = canlab_run_all_tests('Tag', 'Core');           % only tagged 'Core'
results = canlab_run_all_tests('JUnit', 'results.xml');  % also write JUnit XML
```

`canlab_run_all_tests` adds the parent `CanlabCore/` tree to the MATLAB
path with subfolders, so you do not need to set up paths first. SPM and
Neuroimaging_Pattern_Masks must already be on your path if a test
requires them; tests that need them call `assumeNotEmpty(which('...'))`
and skip gracefully when absent.

## Layout

```
Unit_tests/
├── canlab_run_all_tests.m            entry point used locally and in CI
├── helpers/                          shared fixtures
│   └── canlab_get_sample_fmri_data.m
├── fmri_data/                        tests of the @fmri_data class
├── image_vector/                     tests of @image_vector behavior shared by subclasses
├── statistic_image/                  tests of @statistic_image (ttest, threshold)
├── region/                           tests of @region construction and methods
├── atlas/                            tests of @atlas
├── workflows/                        end-to-end pipeline smoke tests
└── old_to_integrate/                 legacy ad-hoc scripts pending rewrite as proper tests
```

`old_to_integrate/` is excluded from discovery by `canlab_run_all_tests`.
Files there should be ported into one of the other folders as proper
`matlab.unittest` tests (and renamed to `canlab_test_*.m`) over time.

## Writing a new test

Function-based tests are the default style. A test file is a function
named `canlab_test_<thing>.m` that declares its sub-tests:

```matlab
function tests = canlab_test_thing
tests = functiontests(localfunctions);
end

function test_some_specific_behavior(tc)
    obj = canlab_get_sample_fmri_data();
    result = method_under_test(obj);
    tc.verifyEqual(size(result.dat), size(obj.dat));
end
```

Conventions:

- One file per method or invariant. File name `canlab_test_<thing>.m`,
  outer function name matching the file. Inner sub-test function names
  follow the matlab.unittest convention `test_<specific_behavior>` —
  they are local to the file and don't need the `canlab_` prefix.
- Use `tc.verify*` (records failure, keeps going) for normal assertions.
  Reserve `tc.assert*` for setup preconditions where the rest of the
  test cannot run.
- Use `tc.assumeNotEmpty(which('foo.nii'))` when a test requires a file
  that may not be on every contributor's path. The test will be marked
  *Incomplete* rather than failed.
- Prefer `canlab_get_sample_fmri_data()` (in `helpers/`) over re-loading
  the sample directly so the sample-load path is centralized.
- A test should not write outside a temp directory. If you need to write,
  use `tempname` and clean up with `tc.addTeardown(@delete, ...)`.
- Tests should not depend on each other or on execution order.

### Why `canlab_test_*` files instead of `test_*`?

`matlab.unittest`'s `TestSuite.fromFolder` only auto-discovers files
whose names start or end with `test` (case-insensitive), so the runner
glob-discovers `canlab_test_*.m` and uses `TestSuite.fromFile` directly.
The result is the same; the difference is that `canlab_test_constructor`
on the path won't shadow some other toolbox's `test_constructor`.

### Tags

Apply tags via class-based tests when a group of tests needs to be
selectable as a unit (e.g., `'RequiresMasks'`, `'Slow'`, `'GUI'`).
Tags filter via `canlab_run_all_tests('Tag', 'Slow')` or
`canlab_run_all_tests('Tag', '~GUI')`.

## CI

Every push and PR runs the suite on the latest supported MATLAB release
on Ubuntu. SPM and Neuroimaging_Pattern_Masks are checked out as
siblings. See `.github/workflows/test.yml` for the full configuration.
