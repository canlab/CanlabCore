# CanlabCore unit tests

This is the unit-test suite for CanlabCore. It uses MATLAB's built-in
`matlab.unittest` framework. There is no external runner — you can run any
single test by `runtests('test_constructor')`, or the whole suite with
`run_all_tests`.

## Running locally

```matlab
cd CanlabCore/Unit_tests
results = run_all_tests;                          % full suite
results = run_all_tests('Tag', 'Core');           % only tagged 'Core' tests
results = run_all_tests('JUnit', 'results.xml');  % also write JUnit XML
```

`run_all_tests` adds the parent `CanlabCore/` tree to the MATLAB path with
subfolders, so you do not need to set up paths first. SPM and
Neuroimaging_Pattern_Masks must already be on your path if a test requires
them; tests that need them call `assumeNotEmpty(which('...'))` and skip
gracefully when absent.

## Layout

```
Unit_tests/
├── run_all_tests.m            entry point used locally and in CI
├── helpers/                   shared fixture loaders (e.g. get_sample_fmri_data)
├── fmri_data/                 tests of the @fmri_data class
├── image_vector/              tests of @image_vector behavior shared by subclasses
├── statistic_image/           tests of @statistic_image (ttest, threshold)
├── region/                    tests of @region construction and methods
├── atlas/                     tests of @atlas
├── workflows/                 end-to-end pipeline smoke tests
└── old_to_integrate/          legacy ad-hoc scripts pending rewrite as proper tests
```

`old_to_integrate/` is excluded from discovery by `run_all_tests`. Files
there should be ported into one of the other folders as proper
`matlab.unittest` tests over time.

## Writing a new test

Function-based tests are the default style. A test file is a function that
declares its sub-tests:

```matlab
function tests = test_thing
tests = functiontests(localfunctions);
end

function test_some_specific_behavior(tc)
    obj = get_sample_fmri_data();
    result = method_under_test(obj);
    tc.verifyEqual(size(result.dat), size(obj.dat));
end
```

Conventions:

- One file per method or invariant. File name `test_<thing>.m`, sub-test
  function name `test_<specific_behavior>`.
- Use `tc.verify*` (records failure, keeps going) for normal assertions.
  Reserve `tc.assert*` for setup preconditions where the rest of the test
  cannot run.
- Use `tc.assumeNotEmpty(which('foo.nii'))` when a test requires a file
  that may not be on every contributor's path. The test will be marked
  *Incomplete* rather than failed.
- Prefer `get_sample_fmri_data()` (in `helpers/`) over re-loading directly
  so the sample-load path is centralized.
- A test should not write outside a temp directory. If you need to write,
  use `tempname` and clean up with `tc.addTeardown(@delete, ...)`.
- Tests should not depend on each other or on execution order.

### Tags

Apply tags via class-based tests when a group of tests needs to be
selectable as a unit (e.g., `'RequiresMasks'`, `'Slow'`, `'GUI'`). Tags
filter via `run_all_tests('Tag', 'Slow')` or `run_all_tests('Tag', '~GUI')`.

## CI

Every push and PR runs the suite on the latest supported MATLAB release on
Ubuntu. SPM and Neuroimaging_Pattern_Masks are checked out as siblings.
See `.github/workflows/test.yml` for the full configuration.
