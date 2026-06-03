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
                                          plot(pm);                 % scores-by-class + ROC
                                          montage(pm, hw_obj, 'use', 'thresh_fdr', 'regions');
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
into `fmri_data.predict` or into the new predictive_model API.

```matlab
X  = double(hw_obj.dat');
Y  = hw_obj.Y;
id = hw_obj.metadata_table.subj_id;   % cell of subject-id strings — used as-is
```

The grouping vector may be numeric or a cell array of id strings; the
splitter normalises it for you (a `grp2idx(...)` is no longer required).

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

### Scoring options

`'scoring'` accepts any registered `cv_scorer` name (or a `cv_scorer`
object). The scorer determines both what `crossval` optimises/reports
and whether continuous scores are needed:

| Name | Task | Needs scores? | Higher better? |
|---|---|---|---|
| `accuracy`, `balanced_accuracy`, `f1` | classification | no | yes |
| `roc_auc`, `log_loss` | classification | yes | yes / no |
| `r2`, `pearson_r` | regression | no | yes |
| `mse`, `rmse`, `mae` | regression | no | no |

`balanced_accuracy` is the safe default for classification (robust to
class imbalance); `roc_auc` is threshold-free and a good choice when
you care about ranking rather than a fixed decision boundary.

> **Tip:** `'groups'` accepts a numeric vector OR a cell array of
> subject-id strings (e.g. `hw_obj.metadata_table.subj_id`) directly —
> the splitter normalises non-numeric labels for you, so the `grp2idx`
> call above is optional.

## 5. Visualise cross-validated performance

The new top-level visualisation methods read straight from the
cross-validated `pm` — no manual extraction:

```matlab
plot(pm);                 % classification: scores-by-class violin + ROC panel
rocplot(pm);              % ROC curve alone; returns AUC / sens / spec struct
confusionchart(pm);       % row-normalized confusion chart
```

`plot(pm)` dispatches on `pm.task` (scatter for regression, the
two-panel violin+ROC for classification). `rocplot` wraps CANlab's
`roc_plot` on the cross-validated decision values; `confusionchart`
wraps MATLAB's confusion chart on `pm.fitted_values.yfit`.

## 6. Visualise the weight map (pre-thresholding)

```matlab
montage(pm, hw_obj);                 % one-line delegate
% equivalently:
si_raw = weight_image(pm, hw_obj);
montage(si_raw);
surface(pm, hw_obj);                 % cortical-surface rendering
```

`weight_image(pm, source)` maps `pm.weights.w` back into the source
fmri_data's voxel space (using its `volInfo` + `removed_voxels`)
and returns a `@statistic_image`. `montage(pm, source)` /
`surface(pm, source)` are thin delegates that build that image and
render it. Voxels the model dropped at fit time (in
`omitted_features`) become zeros.

## 7. Bootstrap inference — and a word about its limit

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
4. **Stability selection** — count how often each voxel lands in the
   top-k by `|weight|` across bootstrap resamples. This is now a
   first-class method, `stability_selection(pm, X, Y, ...)` (see
   §12), and is the recommended inference for high-dimensional
   regularised models where the bootstrap z/p collapses.

## 8. Visualise the thresholded weight map

After bootstrap, the canonical "what does the classifier rely on"
map is `pm.weights.w` masked by `pm.weights.fdr_sig`:

```matlab
montage(pm, hw_obj, 'use', 'thresh_fdr');            % delegate
% or build the image yourself:
si_thr = weight_image(pm, hw_obj, 'use', 'thresh_fdr');
montage(si_thr);
```

Or, since `weight_image` attaches `si.p` and `si.sig`, re-threshold
with the statistic_image methods and outline contiguous clusters:

```matlab
montage(pm, hw_obj, 'use', 'thresh_fdr', 'regions');  % delegate -> region montage
% equivalently:
si = weight_image(pm, hw_obj);
si = threshold(si, .05, 'fdr');
montage(region(si));
```

## 9. Permutation test — is the model itself significant?

```matlab
pm = permutation_test(pm, X, Y, 'nperm', 1000, 'groups', id);
disp(pm.permutation_results.p_value);
disp(pm.permutation_results.permutation);          % which scheme was used
disp(pm.permutation_results.permutation_descrip);  % one-line explanation
```

### Permutation schemes

| `'permutation'` value | What it does | When to use |
|---|---|---|
| `'auto'` *(default)* | Detect from data: no groups → `free`; groups + Y constant per group → `between_subjects`; groups + Y varies per group → `within_subjects` | Safe default; the resolved scheme is stored on the output |
| `'within_subjects'` *(paired-design gold standard)* | Permute Y INDEPENDENTLY within each subject's observations | Each subject contributes both classes (DPSP Hot+Warm, drug A vs B within-subject, etc.). Preserves subject-level pattern; breaks only the class-vs-brain mapping |
| `'between_subjects'` | Reassign each subject's Y to another subject's at random | Each subject contributes ONE class (patient vs control). Errors with a warning if Y varies within a subject |
| `'free'` | Observation-level shuffle of Y, ignoring groups | Truly independent observations only. **Inflates false positives** for grouped data |

The output stores which scheme was actually used in
`pm.permutation_results.permutation` and a one-line explanation in
`pm.permutation_results.permutation_descrip`, so you always know
what null you tested against.

**Why within-subject is the gold standard for paired designs:**
the null is "there's no Hot-vs-Warm signal in the brain for *each*
subject." Within-subject shuffling preserves the subject-level
brain pattern (no subject confound), the class balance (still 59
Hot + 59 Warm globally), and the within-subject correlation
structure; the only thing it breaks is the actual Hot/Warm
assignment for each subject's two observations. That's exactly
the thing the classifier is supposed to be exploiting, so it's
the right thing to randomize.

## 10. Calibrated probabilities

SVM decision values are not probabilities. To get well-calibrated
probabilities, fit a Platt sigmoid on the cross-validated decision
values:

```matlab
pm = calibrate(pm, X, Y);                    % Platt scaling (default)
P  = predict_proba(pm, X_new);               % calibrated class-1 probs in [0,1]
```

`calibrate(..., 'method', 'isotonic')` uses isotonic regression
(pool-adjacent-violators) instead of a parametric sigmoid.

## 11. Hyperparameter search

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

## 12. Univariate feature selection

```matlab
pm = select_features(pm, X, Y, 'k', 1000);   % keep top-1000 by univariate p

pm.diagnostics.feature_selection.n_selected
pm.diagnostics.feature_selection.selected     % logical mask
```

The selection is unioned into `pm.omitted_features`, so subsequent
`predict()` calls on the full-width X automatically drop the
non-selected columns. **Important caveat:** when used in conjunction
with `crossval`, the selection is computed on the *whole* dataset
— this leaks held-out information into feature choice. For honest,
CV-respecting feature selection, wrap the screen as a custom step in
an `@pipeline` (whose `crossval` refits every step on each fold's
training rows only):

```matlab
fs = struct( ...
  'fit_transform', @(Xtr,Ytr) deal(Xtr, top_k_mask(Xtr, Ytr, 2000)), ...
  'transform',     @(mask, Xin) Xin(:, mask));            % apply learned mask
est  = predictive_model('algorithm','svm','task','classification');
pipe = pipeline({fs}, est);
pipe = crossval(pipe, X, Y, 'cv', cv_splitter.stratified_group_kfold(5), 'groups', id);
```

where `top_k_mask` returns a logical column mask of the 2000 most
Y-correlated voxels computed **on the training rows only**.

## 13. Stability selection — voxel inference when bootstrap z/p collapses

As §7 warned, a strongly regularised model gives near-identical
bootstrap weights, so the bootstrap z/p (and the FDR mask) become
useless. **Stability selection** asks a different, more robust
question: *on each bootstrap, which voxels land in the top-k by
`|weight|`, and how often does each voxel make that cut?* Voxels that
are reliably top-ranked across resamples are "stable" — the classifier
leans on them regardless of which subjects it sees.

```matlab
pm = stability_selection(pm, X, Y, ...
    'nboot',     200, ...      % bootstrap resamples
    'k',         2000, ...     % top-k voxels by |w| per bootstrap
    'threshold', 0.6, ...      % "stable" if selected in >= 60% of boots
    'groups',    id);          % resample whole subjects

ss = pm.diagnostics.stability_selection;
fprintf('%d stable voxels (of %d)\n', ss.n_stable, numel(ss.selection_freq));
```

After it runs, `ss` holds:

| Field | What |
|---|---|
| `ss.selection_count` | `[p × 1]` times each voxel was in the top-k |
| `ss.selection_freq`  | `[p × 1]` that count / `valid_boots`, in `[0,1]` |
| `ss.stable`          | logical mask, `selection_freq >= threshold` |
| `ss.n_stable`        | number of stable voxels |
| `ss.valid_boots`     | bootstraps that produced a usable weight vector |

The selection frequency is itself a brain map — push it through
`weight_image` to montage where the classifier is *stably* reading
signal:

```matlab
freq_obj = hw_obj;                       % borrow volInfo + removed_voxels
freq_obj.dat = ss.selection_freq;
montage(freq_obj);                       % selection-frequency map in [0,1]
```

Reference: Meinshausen & Bühlmann, *Stability Selection*, J. R. Stat.
Soc. B (2010).

## 14. End-to-end on the rejection task

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

## 15. Reproduce as a unit test

The full flow is also `CanlabCore/Unit_tests/predictive_model_unit_test.m`,
which runs the entire pipeline with small `nboot=25` / `nperm=10`
so the test completes in a few minutes. Use it as a sanity check
when extending the API.

```matlab
cd CanlabCore/Unit_tests
predictive_model_unit_test
```
