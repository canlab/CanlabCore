# Multivariate classification with SVM — Part 2

> All of Part 1's workflow now in the new sklearn-style
> `@predictive_model` API: construct → fit → cross-validate →
> bootstrap → threshold → visualise, plus calibration, permutation
> testing, hyperparameter search, and feature selection. Same DPSP
> Hot-vs-Warm dataset.

## What's new vs Part 1

Part 1 used `xval_SVM(X, Y, id, ...)` — a single-function wrapper that
returns a `predictive_model` already populated with cv predictions,
weights, and optionally bootstrap/permutation results. Part 2 takes
the same workflow apart into composable methods on the
`@predictive_model` class, so each step is independently
inspectable and re-runnable.

```
% Part 1                                  % Part 2
S = xval_SVM(X, Y, id, ...);              pm = predictive_model('algorithm','svm','task','classification');
                                          pm = crossval(pm, X, Y, 'cv', cv_splitter.stratified_group_kfold(5), 'groups', id);
                                          pm = bootstrap(pm, X, Y, 'nboot', 1000, 'groups', id);
                                          pm = permutation_test(pm, X, Y, 'nperm', 1000, 'groups', id);
                                          si = weight_image(pm, hw_obj);
                                          si = threshold(si, .05, 'fdr');
                                          montage(region(si));
```

Both return the same `@predictive_model` object; the new API just
gives you more control over each step.

## 1. Load DPSP (one line each)

```matlab
hw_obj = load_image_set('DPSP_hotwarm');         % Hot (+1) vs Warm  (-1)
rf_obj = load_image_set('DPSP_rejectorfriend');  % Rejector (+1) vs Friend (-1)
```

Both keywords (a) apply the gray-matter mask, (b) `cat()` + `remove_empty`,
(c) set `±1` effects-coded labels in `.Y`, and (d) attach
`subj_id` to `.metadata_table`. The objects can be passed directly
into `fmri_data.predict` or, with a `grp2idx`, into the new
predictive_model API.

```matlab
X  = double(hw_obj.dat');
Y  = hw_obj.Y;
id = grp2idx(hw_obj.metadata_table.subj_id);
```

## 2. Construct a predictive_model

`predictive_model` is now a thin **value-class** wrapper around the
MATLAB Statistics & ML Toolbox model objects. Construction sets
hyperparameters only — no data has been touched yet.

```matlab
pm = predictive_model( ...
    'algorithm',     'svm', ...                       % registry name
    'task',          'classification', ...
    'modeloptions',  {'KernelFunction', 'linear'}, ...
    'random_state',  2026);                            % rng seed for reproducibility
```

Registry of algorithms (see `predictive_model.algorithm_registry()`):

| Name             | Underlying fitter | Task           |
|------------------|-------------------|----------------|
| `svm`            | `fitcsvm`         | classification |
| `linear_svm`     | `fitclinear`      | classification |
| `logistic`       | `fitclinear`      | classification |
| `lda`            | `fitcdiscr`       | classification |
| `svr`            | `fitrsvm`         | regression     |
| `linear_svr`     | `fitrlinear`      | regression     |

## 3. fit / predict / score (the sklearn triad)

```matlab
pm = fit(pm, X, Y, 'id', id);

[yhat, scores] = predict(pm, X);             % yhat ∈ {-1,+1}, scores = decision values
acc = score(pm, X, Y);                       % balanced_accuracy by default
```

After `fit`:
- `pm.fit_type` is `'insample'`
- `pm.ml_model` is the underlying `ClassificationSVM`
- `pm.weights.w` is the linear-kernel coefficient vector

A pre-fit data-quality pass is automatic — `predictive_model.detect_bad_data`
flags NaN/Inf cases and all-NaN / zero-variance / Inf features.
Bad rows/cols are dropped before training and recorded in
`pm.omitted_cases` and `pm.omitted_features` so `predict()` can
re-apply the same mask to new data.

## 4. Cross-validation — the sklearn-style splitter API

```matlab
pm = crossval(pm, X, Y, ...
    'cv',      cv_splitter.stratified_group_kfold(5), ...
    'groups',  id, ...
    'scoring', 'balanced_accuracy');
```

Splitter factories (`cv_splitter.X(...)`):

| Name                       | Use case                                            |
|----------------------------|-----------------------------------------------------|
| `kfold(k)`                 | plain k-fold                                        |
| `stratified_kfold(k)`      | preserves class proportions per fold                |
| `group_kfold(k)`           | no group spans train/test                           |
| `stratified_group_kfold(k)`| both (recommended for paired classifier designs)    |
| `leave_one_group_out()`    | one fold per group                                  |
| `repeated_kfold(k, n)`     | k folds × n random repeats                          |
| `shuffle_split(n, size)`   | random train/test splits                            |
| `holdout(size)`            | single random split                                 |
| `custom_partition(ids)`    | user-supplied integer fold assignments              |

After `crossval`:
- `pm.fit_type` is `'crossval'`
- `pm.fitted_values.yfit` is the held-out predictions (full length)
- `pm.error_metrics.balanced_accuracy.value` is the mean across folds
- `pm.error_metrics.balanced_accuracy_per_fold.value` is the per-fold vector
- `pm.fold_models{f}` is the trained model from fold f
- `pm.ml_model` is the model re-trained on the FULL dataset (the
  "ship-it" version)

## 5. Visualise the weight map (pre-thresholding)

```matlab
si_raw = weight_image(pm, hw_obj);
create_figure('SVM weights — raw');
montage(si_raw);
```

`weight_image(pm, source)` maps `pm.weights.w` back into the source
fmri_data's voxel space (using its `volInfo` + `removed_voxels`)
and returns a `@statistic_image`. Voxels the model dropped at fit
time (in `omitted_features`) become zeros.

## 6. Bootstrap inference — and a word about its limit

```matlab
pm = bootstrap(pm, X, Y, 'nboot', 1000, 'groups', id);
```

Subjects (not observations) are resampled with replacement when
`groups` is provided. After bootstrap:

| Field                  | What                                               |
|------------------------|----------------------------------------------------|
| `pm.weights.boot_w`    | `[p × nboot]` bootstrap weight samples             |
| `pm.weights.boot_w_mean`/`_ste` | mean and SD of those samples              |
| `pm.weights.z`         | mean / SE_floored                                  |
| `pm.weights.p`         | **continuity-corrected empirical two-tailed p**    |
|                        | bounded in `[2/(nboot+1), 1]`                      |
| `pm.weights.p_wald`    | Wald-style `2*(1-normcdf(|z|))` (legacy formula)   |
| `pm.weights.fdr_thr`   | FDR(p, .05) on the empirical p                     |
| `pm.weights.fdr_sig`   | logical mask, p ≤ fdr_thr                          |
| `pm.weights.thresh_fdr`| `pm.weights.w` masked by `fdr_sig`                 |

**Caveat — read this if your bootstrap results look "everything significant" or "nothing significant":**
L2-regularised linear SVM (the `linear_svm` algorithm and the
`'highdimensional', true` path in xval_SVM) gives bootstrap
weights that are *nearly numerically identical across resamples* —
the regularisation drives the optimiser to almost the same solution
regardless of which 118 subjects (with replacement) you fit on. As
a consequence:

- The Wald z statistic blows up (numerator ≈ small, denominator ≈
  numerical zero), and the legacy `2*(1-normcdf(|z|))` p collapses to
  `eps`, so naïve FDR flags ~every voxel as "significant".
- The empirical bootstrap p we now use floors at `2/(nboot+1)` —
  honest, but it means for `linear_svm` with default Lambda, *every
  voxel hits the floor* and the FDR threshold collapses.

This is a real limit of bootstrap inference on a strongly
regularised model. To get finer-grained voxel-wise inference:

1. **Reduce regularisation** — pass `'Lambda', 1e-6` in
   `modeloptions` so the model has room to vary across bootstraps.
2. **Drop to a kernel SVM** — `'algorithm', 'svm'` (the default,
   non-`highdimensional` path) uses `fitcsvm`'s dual formulation,
   which produces more variable weights at the cost of slower
   fitting on wide data.
3. **Use a permutation test on a univariate statistic** —
   `permutation_test(pm, X, Y)` on the cross-validated overall score
   is a clean alternative for asking "is the *model* significantly
   better than chance?" rather than "which voxels are significantly
   non-zero?".
4. **Stability selection** — feature selection by counting how often
   a voxel appears in the top-k across bootstraps. Not yet a
   first-class method on the class but easy to compute from
   `pm.weights.boot_w`.

## 7. Visualise the thresholded weight map

After bootstrap, the canonical "what does the classifier rely on"
map is `pm.weights.w` masked by `pm.weights.fdr_sig`:

```matlab
si_thr = weight_image(pm, hw_obj, 'use', 'thresh_fdr');
create_figure('SVM weights — FDR thresholded');
montage(si_thr);
```

Or equivalently, `weight_image` attaches `si.p` and `si.sig` so
you can re-threshold with the statistic_image methods:

```matlab
si = weight_image(pm, hw_obj);
si = threshold(si, .05, 'fdr');
montage(region(si));
```

## 8. Permutation test — is the model itself significant?

```matlab
pm = permutation_test(pm, X, Y, 'nperm', 1000, 'groups', id);
disp(pm.permutation_results.p_value);
```

For paired-within-subject designs (DPSP Hot+Warm per subject), the
permutation is done WITHIN each subject — Y is shuffled for each
subject's observations independently, preserving the subject
structure while breaking the class-vs-brain mapping. For
between-subjects designs (each subject has one class), the
permutation happens at the subject level.

## 9. Calibrated probabilities

SVM decision values are not probabilities. To get well-calibrated
probabilities, fit a Platt sigmoid on the cross-validated decision
values:

```matlab
pm = calibrate(pm, X, Y);                    % Platt scaling (default)
P  = predict_proba(pm, X_new);               % calibrated class-1 probs in [0,1]
```

`calibrate(..., 'method', 'isotonic')` uses isotonic regression
(pool-adjacent-violators) instead of a parametric sigmoid.

## 10. Hyperparameter search

```matlab
pg = struct();
pg.BoxConstraint = [0.1 1 10 100];
pg.KernelScale   = [0.5 1 2 5];

pm = grid_search(pm, X, Y, pg, 'groups', id);

pm.diagnostics.grid_search.best_params
pm.diagnostics.grid_search.best_score
pm.diagnostics.grid_search.scores            % full grid
```

`grid_search` argmaxes (or argmins, for less-is-better scorers like
`mse`) the mean cv score across the Cartesian product of `pg`
values, then refits at the winner.

## 11. Univariate feature selection

```matlab
pm = select_features(pm, X, Y, 'k', 1000);   % keep top-1000 by univariate p

pm.diagnostics.feature_selection.n_selected
pm.diagnostics.feature_selection.selected     % logical mask
```

The selection is unioned into `pm.omitted_features`, so subsequent
`predict()` calls on the full-width X automatically drop the
non-selected columns. **Important caveat:** when used in conjunction
with `crossval`, the selection is computed on the *whole* dataset
— this leaks held-out information into feature choice. For honest
CV-respecting feature selection, run `select_features` *per fold*
(future feature: a Pipeline composer).

## 12. End-to-end on the rejection task

The same pipeline on `DPSP_rejectorfriend`:

```matlab
rf_obj = load_image_set('DPSP_rejectorfriend');
X  = double(rf_obj.dat');
Y  = rf_obj.Y;
id = grp2idx(rf_obj.metadata_table.subj_id);

pm = predictive_model('algorithm','svm','task','classification', ...
                       'modeloptions', {'KernelFunction','linear'}, ...
                       'random_state', 2026);

pm = crossval(pm, X, Y, 'cv', cv_splitter.stratified_group_kfold(5), ...
              'groups', id, 'scoring', 'balanced_accuracy');

fprintf('Rejector-vs-Friend cv balanced accuracy: %.1f%%\n', ...
    100 * pm.error_metrics.balanced_accuracy.value);
```

## 13. Reproduce as a unit test

The full flow is also `CanlabCore/Unit_tests/predictive_model_unit_test.m`,
which runs the entire pipeline with small `nboot=25` / `nperm=10`
so the test completes in a few minutes. Use it as a sanity check
when extending the API.

```matlab
cd CanlabCore/Unit_tests
predictive_model_unit_test
```
