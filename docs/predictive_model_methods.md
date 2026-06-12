# `predictive_model` methods, organized by area

`predictive_model` is a scikit-learn-style object for multivariate
predictive modelling. It holds the hyperparameters that *specify* a model
(algorithm, task, cross-validation scheme, scorer, …) and, after fitting,
the *artifacts* it produces (predictions, weight maps, bootstrap
statistics, error metrics, the underlying trained MATLAB model). It is the
canonical output type of `xval_SVM`, `xval_SVR`, `fmri_data.predict`
(`'newapi'`), and friends; the constructor also accepts a plain struct
(such as the legacy wrapper outputs) and routes each field to its
categorised home, with translations for legacy names like
`SVMModel`/`SVRModel`, `vox_weights`, `wZ`/`wP`, etc.

Type `methods(my_obj)` in MATLAB for the live list on any instance.

## Three ways to create one

You never have to build a `predictive_model` field by field — there are three
entry points, and they all return the **same** object type, so everything in
this reference applies regardless of how you made it.

1. **Standalone, scikit-learn style** — construct from hyperparameters and
   drive it with the method verbs (`fit` → `predict` → `crossval` →
   `bootstrap`/`permutation_test` → `weight_map_object`/`plot`). Most control;
   operates on numeric `(X, Y, groups)`.

   ```matlab
   pm = predictive_model('algorithm','svm','task','classification');
   pm = crossval(pm, X, Y, 'groups', id);
   pm = bootstrap(pm, X, Y, 'nboot', 1000, 'groups', id);
   pm = weight_map_object(pm, dat);     % attach a brain map for montage(pm)
   summary(pm);
   ```

2. **From an `fmri_data` object** — `fmri_data.predict(..., 'newapi')` takes a
   brain-data object directly, runs the cross-validated prediction, and returns
   a `predictive_model` as its 4th output (already carrying a weight map):

   ```matlab
   [cverr, stats, optout, pm] = predict(dat, 'algorithm_name','cv_svm', ...
                                        'nfolds', 5, 'newapi');
   ```

3. **Via an `xval_*` wrapper** — single-call functions that take numeric
   matrices and return a populated `predictive_model`: `xval_SVM`, `xval_SVR`,
   `xval_discriminant_classifier`, `xval_lasso_brain`, `xval_cross_classify`,
   …. These never see an image, so attach a weight map afterward with
   `weight_map_object(pm, reference_image)`:

   ```matlab
   pm = xval_SVM(X, Y, id);
   pm = weight_map_object(pm, dat);     % then montage(pm) / surface(pm)
   ```

The constructor also accepts a **plain struct** (e.g. a legacy wrapper's STATS
output): `predictive_model(stats_struct)` maps each field to its consolidated
home.

## Tutorials

Hands-on walkthroughs in the multivariate-decoding series — each has a
rendered Markdown page and a runnable plain-text Live Script (`.m`, opens in
the Live Editor on R2025a+):

| Part | Topic | Markdown | Live Script (.m) |
|---|---|---|---|
| 1 | Classification basics with SVM (ROC, confusion, effect sizes, weight maps, apply-to-test) | [part 1 .md](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part1_classification_with_SVM.md) | [.m](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part1_classification_with_SVM.m) |
| 2 | Classification **and** regression; one-line loaders; `xval_*` family; predicted R² | [part 2 .md](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part2_classification_and_regression.md) | [.m](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part2_classification_and_regression.m) |
| 3 | The sklearn-style `predictive_model` API (fit/predict/crossval/bootstrap/permutation, tuning, calibration, stability selection) | [part 3 .md](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part3_predictive_model_api.md) | [.m](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part3_predictive_model_api.m) |
| 4 | Cross-classification (Woo et al., 2014: does a pain pattern decode rejection?) | [part 4 .md](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part4_cross_classification.md) | [.m](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part4_cross_classification.m) |
| 5 | Algorithms, multiclass ECOC, tuning, and inference | [part 5 .md](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part5_algorithms_and_tuning.md) | [.m](markdown_tutorials/multivariate_classification_with_SVM/multivariate_decoding_part5_algorithms_and_tuning.m) |

## Design

- **Value semantics.** Every mutating method returns a NEW object — write
  `pm = fit(pm, X, Y)`, not `fit(pm, X, Y)`.
- **Hyperparameters vs fitted state are disjoint.** `clone(pm)` returns a
  fresh, unfitted copy with identical hyperparameters (used for per-fold
  refitting). `is_fitted(pm)` tests for any populated fitted state.
- **Numeric core, image adapters on top.** `fit`/`predict`/`crossval`
  operate on numeric `(X, Y, groups)`. Neuroimaging awareness lives in the
  adapter methods (`weight_map_object`, `montage`, `surface`) that map the
  coefficient vector back into voxel space using a source image's
  `volInfo`.

## Properties (consolidated layout)

### Hyperparameters (set by user / constructor)

| Property | Description |
|---|---|
| `algorithm` | Registry name: `svm`, `linear_svm`, `logistic`, `lda`, `ecoc`, `svr`, `lasso`, `ridge`, `gp`, `pcr` (PCA+OLS), `lassopcr` (PCA+LASSO+relaxed-OLS), … |
| `task` | `'classification'` or `'regression'` |
| `modeloptions` | Cell of fitter options, e.g. `{'KernelFunction','linear'}` |
| `cv` | A `cv_splitter` (kfold, stratified_group_kfold, …) |
| `scorer` | A `cv_scorer` (accuracy, roc_auc, r2, rmse, …) |
| `random_state`, `standardize`, `use_parallel`, `nboot`, `nperm`, `do_calibrate` | reproducibility / behaviour flags |
| `class_labels`, `Y_name`, `X_name` | labels for plotting |

### Data + fit metadata

| Property | Description |
|---|---|
| `Y` | Outcome actually fit (after bad-case removal) |
| `id` | Grouping/subject vector |
| `omitted_cases` / `omitted_features` | Logical masks (over the original cases / features) of what was dropped |
| `inputParameters` | Struct of stored fit-time parameters (standardization mean/std, …) |
| `fit_type` | `'crossval'` \| `'insample'` \| `'test'` — provenance of `yfit` AND `weights` |

### Fitted state (categorised sub-structs)

| Property | Key fields |
|---|---|
| `fitted_values` | `.yfit`, `.scores`, `.score_type`, `.dist_from_hyperplane_xval`, within-person arrays, `.calibrator` |
| `weights` | `.w`, `.intercept`, `.w_perfold`, `.boot_w*`, `.z`, `.p`, `.fdr_thr`, `.fdr_sig`, `.thresh_fdr`, `.weight_obj` |
| `error_metrics` | `(value, descrip)` tuples: `.crossval_accuracy`, `.r2`, `.rmse`, `.prediction_outcome_r`, `.d_within`, … |
| `cv_partition` | `.trIdx`, `.teIdx`, `.nfolds` |
| `ml_model` / `fold_models` | trained MATLAB model(s) |
| `bootstrap_results` / `permutation_results` | full bootstrap / permutation output |
| `diagnostics` | `.mult_obs_within_person`, `.grid_search`, `.feature_selection`, `.stability_selection`, … |
| `cross_classify` | cross-classification output (from `xval_cross_classify`) |
| `note` / `history` | cell arrays of commentary / provenance lines |

`error_metrics` entries are `(value, descrip)` tuples; read e.g.
`pm.error_metrics.crossval_accuracy.value`.

## Fit / predict / score

| Method | One-liner |
|---|---|
| `predictive_model` | Constructor: name/value hyperparameters, or a legacy struct |
| `fit` | Train on the full sample (in-sample fit) |
| `predict` | Apply the fitted model: returns `[yhat, scores]` |
| `score` | Evaluate `obj.scorer` on `(Y, predict(obj, X))` |
| `crossval` | Cross-validate: refit per fold, predict held-out, score; auto within-person stats when `groups` given |
| `clone` / `is_fitted` / `is_classifier` / `is_regressor` / `validate_object` | model bookkeeping |

## Inference / tuning

| Method | One-liner |
|---|---|
| `bootstrap` | Bootstrap weights → SE, Z, empirical p, FDR threshold + significant mask |
| `permutation_test` | Null distribution by shuffling `Y` (free / between- / within-subjects schemes) |
| `grid_search` | Exhaustive hyperparameter search via CV |
| `stability_selection` | Per-bootstrap top-k selection frequency (high-dim feature inference) |
| `select_features` | Univariate / top-k-correlation feature pre-screen |
| `calibrate` / `predict_proba` | Platt / isotonic probability calibration and calibrated `P(class)` |

## Reporting

| Method | One-liner |
|---|---|
| `report_accuracy` | Model-type-relevant performance metrics as a struct + printed table. Classification: accuracy / balanced acc / AUC / sensitivity / specificity / PPV / NPV / *d* / forced-choice acc (from the cross-validated ROC). Regression: prediction–outcome *r*, predicted R², out-of-sample R², RMSE / MAE / MSE |
| `summary` | Human-readable overview: identity (algorithm, task, fit provenance), data (n obs / features / groups), CV scheme, which inference is available (bootstrap / permutation / calibration / cross-classification), and the `report_accuracy` performance block |

**Regression R² variants** (Wager & Lindquist, *Principles of fMRI* Ch. 39.4), computed in `crossval` from the pooled cross-validated predictions and exposed in `error_metrics`:

- `predicted_r2` = 1 − PRESS / SST, where PRESS = Σ(yᵢ − ŷᵢ)² over the held-out (cross-validated) predictions and SST = Σ(yᵢ − ȳ)² uses the grand mean. This is §39.4's *predicted R²*.
- `out_of_sample_r2` = 1 − PRESS / Σ(yᵢ − ȳ_train(fold))², using each fold's *training* mean as the baseline.
- Both can be **negative** when the model predicts worse than the mean. Distinct from the per-fold-averaged `r2` cv_scorer (each fold's own test mean, then averaged — akin to scikit-learn's default).

## Visualization

| Method | One-liner |
|---|---|
| `plot` | Task-dispatched core figure (predicted-vs-observed scatter for regression; scores-by-class + ROC for classification), **plus** a weight-map montage and permutation-null histogram when those artifacts are present (`'noweights'` / `'noperm'` / `'surface'` to control) |
| `plot_predicted_vs_observed` | Scatter (regression) or violin (binary) of predictions vs `Y` |
| `plot_permutation` | Histogram of the permutation null with observed-score and α-level critical-value lines (after `permutation_test`) |
| `rocplot` | ROC curve from cross-validated scores (wraps `roc_plot`) |
| `confusionchart` / `confusion_matrix` | Confusion chart / raw + normalized confusion matrices |
| `weight_map_object` | Map the weight vector into a `statistic_image` (with bootstrap `.p` / `.sig`) and cache it on `pm.weights.weight_obj`; second output is the image (`[~, si] = ...`) |
| `montage` / `surface` | Render the weight map on slices / cortical surface (thin delegates) |

## Cross-validation infrastructure

These standalone classes back `crossval` and are shared with `@pipeline`:

| Class | One-liner |
|---|---|
| `cv_splitter` | Fold generators: `kfold`, `stratified_kfold`, `group_kfold`, `stratified_group_kfold`, `logo`, `holdout`, `shuffle_split`, `repeated_kfold`, `custom_partition` |
| `cv_scorer` | Metrics: `accuracy`, `balanced_accuracy`, `f1`, `roc_auc`, `log_loss`, `mse`, `rmse`, `mae`, `r2`, `pearson_r` |

## Pipelines

`@pipeline` composes zero or more transform steps (`center`, `zscore`,
`pca`, or a custom `fit_transform`/`transform` struct) feeding a final
`predictive_model` estimator. Its `crossval` refits **every step per fold**
(leakage-free PCR, standardize-then-SVM, …), and `weight_map_object(pipe,
source)` back-projects the estimator weights through the steps into voxel
space.

## Worked example

```matlab
dat = load_image_set('DPSP_hotwarm');          % 118 imgs (59 subj x hot/warm), Y = +/-1
X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;

pm = predictive_model('algorithm','svm','task','classification');
pm = crossval(pm, X, Y, 'groups', id, ...
              'cv', cv_splitter.stratified_group_kfold(5));
pm.error_metrics.crossval_accuracy.value        % cross-validated accuracy (%)

pm = bootstrap(pm, X, Y, 'nboot', 1000, 'groups', id);
plot(pm);                                       % scores-by-class + ROC
montage(pm, dat, 'use', 'thresh_fdr', 'regions');  % FDR-thresholded clusters
```
