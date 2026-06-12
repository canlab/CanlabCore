# `predictive_model` — development notes & continuation plan

**Last updated:** 2026-06-09
**Branch:** `predictive_model_obj_dev`  ·  **PR:** [#76 → master](https://github.com/canlab/CanlabCore/pull/76) (open)
**Companion:** the original design plan is `predictive_model_plan_claude.txt` (root). This file is the *current-state + where-to-go-next* handoff.

---

## 0. TL;DR — where things stand

The scikit-learn–style `@predictive_model` API is **built and working** end-to-end, with a 5-part tutorial series, a methods reference, and a 17-section unit test. PR #76 collects the whole branch for merge into `master`. The original plan (Phases 0–6 in `predictive_model_plan_claude.txt`) is essentially **complete**. What remains is review/merge, optional polish, and the future-work backlog in §7.

If you are resuming cold: read §3 (architecture), then §6 (how to verify nothing is broken), then pick from §7.

---

## 1. What the branch delivers (for reviewers)

A **value-class** `@predictive_model` object + verb API for multivariate classification/regression on brain data:

- **Object:** consolidated property layout — hyperparameters (`algorithm`, `task`, `modeloptions`, `cv`, `scorer`, `random_state`, …) kept disjoint from fitted state in categorised sub-structs (`fitted_values`, `weights`, `error_metrics`, `cv_partition`, `ml_model`, `fold_models`, `bootstrap_results`, `permutation_results`, `diagnostics`). Legacy flat `Dependent` aliases were **removed**.
- **Three creation/training paths:** standalone sklearn API; `fmri_data.predict(..., 'newapi')` (4th output is the `pm`); and the refactored `xval_*` wrappers (all return a `predictive_model`).
- **Verbs:** `fit`/`predict`/`score`, `crossval`, `bootstrap`, `permutation_test`, `grid_search`, `calibrate`/`predict_proba`, `select_features`, `stability_selection`.
- **CV infrastructure:** `cv_splitter`, `cv_scorer`, and `@pipeline` (leakage-free per-fold transform refitting).
- **Metrics/reporting:** predicted R² & out-of-sample R² (Ch. 39.4), `report_accuracy`, `summary`.
- **Viz:** `plot` (auto-adds weights + permutation panels), `plot_permutation`, `rocplot`, `confusionchart`/`confusion_matrix`, `weight_map_object` (+ `montage`/`surface` delegates).

---

## 2. What *this* session added (most recent work)

1. **`weight_map_object(pm, reference)`** — builds the weight `statistic_image` and **caches it on `pm.weights.weight_obj`** so `montage(pm)`/`surface(pm)` need no source. Consolidated the old `weight_image` builder into it (single method, `[pm, si]` return). Wired into `fmri_data.predict 'newapi'`, `@pipeline`, and `xval_*` help.
2. **Predicted R² / out-of-sample R²** (Wager & Lindquist, *Principles of fMRI* Ch. 39.4) computed in `crossval` via static `predictive_model.regression_metrics_from_cv`. Plus **`report_accuracy`** (metric struct + table) and **`summary`** (provenance + metrics).
3. **`plot_permutation`** method + **`plot(pm)`** enhanced to auto-draw the weight montage and permutation histogram when present.
4. **Effect-size investigation (resolved):** within-person Cohen's *dz* being *below* the single-interval *d* while forced-choice accuracy is *higher* is **correct, not a bug** — `dz/d = 1/√(2(1−r))` with within-subject r≈0.39 (<0.5); accuracy compares Φ(dz) vs Φ(d/2). Documented in Part 1.
5. **Docs:** Part 3 (full algorithm registry incl. `pcr`/`lassopcr`, rewritten Scoring section, Repeated-CV section, thresholded montage/surface + permutation figures, completed shrinkage-sweep code, `summary` demo); Part 1 (weight-map + summary sections, twochoice ROC output); methods.md (creation-paths section + tutorial links table); constructor/static-method help. Repointed help links from `CANlab_help_examples`/readthedocs → `canlab.github.io`. Synced Part 1 & 3 `.m` Live Scripts.

---

## 3. Architecture / where things live

```
CanlabCore/@predictive_model/         <- the object + methods (one .m per method)
  predictive_model.m                  <- classdef: properties, constructor, ALL static methods
  fit.m predict.m score.m crossval.m
  bootstrap.m permutation_test.m grid_search.m
  calibrate.m predict_proba.m select_features.m stability_selection.m
  weight_map_object.m montage.m surface.m
  plot.m plot_permutation.m plot_predicted_vs_observed.m rocplot.m
  confusion_matrix.m confusionchart.m report_accuracy.m summary.m

CanlabCore/@pipeline/                  <- transform-steps + estimator composer
  pipeline.m crossval.m

CanlabCore/Statistics_tools/Cross_validation/
  @cv_splitter/cv_splitter.m           <- fold factories (kfold, group_kfold, ...)
  @cv_scorer/cv_scorer.m               <- metrics (balanced_accuracy, r2, roc_auc, ...)

CanlabCore/@fmri_data/predict.m        <- 'newapi' route -> @predictive_model (4th output)
CanlabCore/Statistics_tools/Cross_validated_Regression/
  xval_SVM.m xval_SVR.m ...            <- thin wrappers returning a predictive_model
CanlabCore/Statistics_tools/xval_discriminant_classifier.m

CanlabCore/Unit_tests/
  predictive_model_unit_test.m         <- 17-section end-to-end test
  compare_predict_legacy_vs_newapi.m   <- legacy vs newapi parity diagnostic

docs/predictive_model_methods.md       <- methods reference + creation paths + tutorial table
docs/markdown_tutorials/multivariate_classification_with_SVM/
  multivariate_decoding_part1..5{.md,.m}   <- 5-part tutorial (md + Live Script .m)
  pngs/                                 <- tutorial figures
  make_predictive_model_codemap.m       <- regenerates the API code-map PNG
predictive_model_plan_claude.txt        <- original design plan (Phases 0-6)
```

**Key idea to remember:** the numeric core (`fit`/`predict`/`crossval`) knows nothing about brains; neuroimaging awareness lives only in the *adapter* methods (`weight_map_object`, `montage`, `surface`) that use a source image's `volInfo` to map weights back to voxels.

---

## 4. Invariants & design decisions (don't break these)

- **Value semantics.** Every mutating method returns a NEW object: `pm = method(pm, ...)`. Never rely on in-place mutation.
- **Hyperparameters vs fitted state are disjoint.** `clone(pm)` = fresh unfitted copy with same hyperparameters (used for per-fold refit). `is_fitted(pm)` tests for any fitted state.
- **One canonical path per datum.** No flat aliases. Read `pm.fitted_values.yfit`, `pm.weights.w`, `pm.error_metrics.<name>.value` (metrics are `(value, descrip)` tuples — use `predictive_model.metric_value(...)` to unwrap).
- **`pcr`/`lassopcr` are special-cased in `fit.m`** (not in `algorithm_registry`), and are verified to reproduce the legacy `cv_pcr`/`cv_lassopcr` predictions. Don't "simplify" them into the registry without re-checking parity against the legacy path.
- **Bootstrap on strongly-regularised `linear_svm`** gives near-identical weights across resamples → empirical p floors at `2/(nboot+1)`, FDR can be empty. This is expected; the recommended inference there is `stability_selection` (documented in Part 3 §7, §13).
- **Group labels** may be numeric or cell-string; the splitter normalises them. Don't add `==`-on-cellstring logic.

---

## 5. Two subtle results worth remembering

- **Effect sizes vs accuracy (Part 1):** within-person *dz* can be *below* single-interval *d* while forced-choice accuracy is *above* single-interval accuracy. Not a bug — different denominators (`dz/d = 1/√(2(1−r))`) and different accuracy maps (Φ(dz) vs Φ(d/2)). Three conventions float around: object `d_singleinterval`, object `d_within` (paired dz), and `roc_plot 'twochoice'` d_a (= √2·z(acc)). They legitimately differ.
- **Predicted R² ≠ per-fold-averaged `r2` scorer.** `predicted_r2` pools all held-out predictions against the grand mean (PRESS/SST); the `r2` `cv_scorer` averages per-fold R² (each fold's own mean), akin to scikit-learn's default. Both are exposed; keep them distinct.

---

## 6. How to verify nothing is broken (resume checklist)

```matlab
% from a directory above the repos:
canlab_toolbox_setup            % path setup

cd CanlabCore/Unit_tests
predictive_model_unit_test      % should print "PASS" (17 sections)

% parity check (legacy vs newapi), folds held constant:
compare_predict_legacy_vs_newapi
```

**Environment gotchas (MCP MATLAB / interactive):**
- After editing a class method file, run `clear classes; rehash toolboxcache;` before re-testing — MATLAB caches the old method.
- Lingering `onCleanup` objects can block `clear classes` (warning is harmless).
- `run_matlab_file` in a long-lived MATLAB session has served a **stale cached parse** of a just-edited script; prefer `evaluate_matlab_code` (or a fresh session) and `check_matlab_code` for static validation.
- Tutorial `.m` files intentionally leave some lines without a trailing `;` (to echo output in the Live Editor) — the "add semicolon" code-analyzer hints on those are expected, not errors.

---

## 7. Future development backlog (pick from here)

**A. Merge & release**
- Get PR #76 reviewed and merged to `master`. After merge, the root scratch files (`predictive_model_plan_claude.txt`, `predictive_model_update1_claude.txt.rtf`, `test_pm_code.m`, `bootstrap_weightmap_output.txt`) should be triaged — keep the plan, drop the scratch.
- Add `@predictive_model` and the cv classes to `docs/Object_methods.md` so they appear in the object index.

**B. API completeness / quality**
- **Multiclass beyond ECOC:** `report_accuracy`/`rocplot`/`confusionchart` assume binary for the ROC-derived metrics. Add graceful multiclass handling (per-class one-vs-rest AUC, macro-averaged sens/spec).
- **`pm.plot` for regression:** currently scatter + (conditional) weights/perm. Consider adding a residual plot and a calibration plot.
- **Nested CV as a first-class method** (`nested_crossval`) rather than the hand-rolled loop in Part 3 §11 — returns an unbiased score plus the per-fold chosen hyperparameters.
- **`save`/`load` ergonomics:** confirm a fitted `pm` (with cached `weight_obj`) round-trips through `save`; document recommended persistence.
- **Calibration metrics:** Brier score / reliability curve in `report_accuracy` when `calibrate` has run.

**C. Performance / scale**
- Bootstrap and permutation loops are the slow paths on wide data. Consider optional `parfor` (there's a `use_parallel` hyperparameter — wire it through `bootstrap`/`permutation_test`/`grid_search`).
- `gp` (Gaussian process) and kernel methods need a PCA step for ~200k voxels; the `@pipeline` path handles it but isn't the default. Consider an auto-reduce guard with a warning.

**D. Inference**
- `stability_selection`: add a calibrated selection threshold (e.g., per Meinshausen–Bühlmann error control) rather than a fixed frequency cutoff.
- Permutation: add cluster-level / max-statistic correction option for the *weight map* (currently permutation is on the overall CV score only).

**E. Docs / outreach**
- Sync any remaining tutorial `.m` files if Parts 2/4/5 `.md` get further edits (Parts 1 & 3 are synced).
- The X-post draft lives at `docs/X_post_draft_pm.rtf` (intentionally **not committed**).
- Consider a short "migration guide" for users with scripts that used the removed flat aliases (`pm.yfit` → `pm.fitted_values.yfit`, etc.).

---

## 8. Concrete first step when resuming

1. `git checkout predictive_model_obj_dev && git pull`
2. Run the resume checklist (§6). If the unit test passes, the branch is healthy.
3. Check PR #76 status; if merged, branch from `master`. If not, address review comments.
4. Pick a backlog item from §7 (highest leverage: **B-nested CV method** and **B-multiclass metrics**, since those are the two places the current API is thinnest).
