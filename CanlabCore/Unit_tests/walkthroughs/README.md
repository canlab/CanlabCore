# Walkthrough integration tests (Tier B)

These tests smoke-test the end-to-end CANlab tutorials ("walkthroughs") from
the [CANlab_help_examples](https://github.com/canlab/CANlab_help_examples)
repository, to catch the case where a tutorial silently bit-rots against the
current toolbox.

They are **slow** and depend on more sibling repos and optional data than the
fast per-push unit suite, so they are a separate tier, run on the nightly
`tests-walkthroughs` GitHub Actions workflow, **not** on every push. They are
skipped by `canlab_run_all_tests` unless you ask for them:

```matlab
canlab_run_all_tests('Walkthroughs', 'only')      % just the walkthroughs
canlab_run_all_tests('Walkthroughs', 'include')   % unit suite + walkthroughs
canlab_run_all_tests                              % default: walkthroughs skipped
```

## Design: snapshots + a hardening harness

The tutorials in CANlab_help_examples are teaching scripts. They should stay
clean tutorials — full of `orthviews`, `surface`, optional signature/atlas
data, and the occasional interactive prompt — and **not** be retrofitted with
test scaffolding. CanlabCore CI should also not depend on an external repo's
un-hardened scripts at run time. So instead of executing the tutorials in
place, this folder works from frozen copies:

```
walkthroughs/
  canlab_test_walkthrough_*.m     <- thin matlab.unittest wrappers (discovered)
  private/
    canlab_help_*.m               <- VERBATIM snapshots of the tutorials
  README.md
```

- `private/` is excluded from `genpath`, so the snapshots never go on the
  MATLAB path and never collide with the real tutorials when both repos are
  checked out on CI. The harness reads them by absolute path.
- Each `canlab_test_walkthrough_*.m` wrapper just points at its snapshot and
  calls the shared harness. The hardening lives in one place
  (`Unit_tests/helpers/canlab_run_walkthrough_snapshot.m`), not duplicated into
  each copy — so refreshing a snapshot is a clean file overwrite.

### What the harness does

`canlab_run_walkthrough_snapshot(tc, snapshot_path)`:

1. Forces `DefaultFigureVisible 'off'` (headless).
2. Splits the snapshot at its `%%` cell boundaries and runs each cell in its
   own `try/catch`, in one shared workspace, so a graphics-only section that
   fails on a headless runner does not abort the compute sections.
3. Classifies any caught error via
   `Unit_tests/helpers/canlab_classify_environment_error.m`:
   - **graphics / input / data / cascade** → section **skipped** (recorded),
     execution continues.
   - **genuine** → the harness stops and the test **fails** with an informative
     report (which section, the offending line, error id and message).
4. Maps outcomes to `matlab.unittest`:
   - any genuine failure → **Failed**
   - every section skipped for environment reasons → **Incomplete**
   - at least one section ran → **Passed** (with a logged skip summary)

So the nightly goes red only on a *genuine* regression, and tells you exactly
where; missing-display / missing-data / prompt conditions become skips.

## Refreshing a snapshot

When a tutorial changes upstream, re-copy it verbatim:

```bash
cp ../../../../CANlab_help_examples/example_help_files/canlab_help_5_regression_walkthrough.m \
   private/canlab_help_5_regression_walkthrough.m
```

No re-hardening needed — the harness handles it. Keep the copies verbatim so
they can be diffed against upstream.

## Why a separate nightly tier (not folded into the per-push suite)

Measured wall-time for the full set is ~5 min on a fast workstation and an
estimated ~15–20 min on the 2-core GitHub Actions Linux runner. The per-push
`tests` suite is the fast gate that every PR blocks on; adding 15–20 min plus
graphics/data-dependent flakiness to it would slow iteration and make the
required check unreliable. End-to-end tutorial runs are exactly the kind of
slow, environment-sensitive integration test that belongs in a nightly tier.
