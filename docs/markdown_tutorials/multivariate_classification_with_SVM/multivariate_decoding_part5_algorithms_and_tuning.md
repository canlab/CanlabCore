# Multivariate decoding — Part 5: algorithms, tuning, and inference

> **Multivariate decoding tutorial series**
> 1. [Classification basics with SVM](multivariate_decoding_part1_classification_with_SVM.md) — train and cross-validate a linear SVM (Hot vs Warm); ROC, confusion matrix, effect sizes; apply to a held-out test set.
> 2. [Classification and regression](multivariate_decoding_part2_classification_and_regression.md) — the difference between the two, the one-line dataset loaders, the `xval_*` wrapper family, and `fmri_data.predict` end-to-end for both.
> 3. [The sklearn-style `predictive_model` API](multivariate_decoding_part3_predictive_model_api.md) — fit / predict / crossval / bootstrap / permutation, nested-CV tuning, calibration, stability selection.
> 4. [Cross-classification](multivariate_decoding_part4_cross_classification.md) — does a pain pattern decode social rejection? (Woo et al., 2014).
> 5. **Algorithms, tuning, and inference** *(this part)* — compare SVM / SVR / lasso / ridge / GP, ECOC multiclass, grid search, stability selection.

> One dataset, many estimators. This part uses the `@predictive_model`
> registry to run **binary** and **multiclass (ECOC)** classification and
> **regression** (SVR, lasso, ridge, Gaussian process) on the DPSP data,
> compares cross-validated performance across algorithms, tunes
> hyperparameters with `grid_search`, and does high-dimensional feature
> inference with `stability_selection` — then closes with a decision
> guide for *when to use which*.

## 1. The algorithm registry

`fit` dispatches `obj.algorithm` through a one-row-per-algorithm
registry (`predictive_model.algorithm_registry()`). Switching estimators
is a one-word change; everything else (cv, scoring, bootstrap, plotting)
is identical.

| Name | Fitter | Task | Notes |
|---|---|---|---|
| `svm` | `fitcsvm` | classification | kernel SVM (linear default); dual form |
| `linear_svm` | `fitclinear` | classification | high-dim L2-SVM; fast on wide data |
| `logistic` | `fitclinear` | classification | logistic loss; gives probabilities |
| `lda` | `fitcdiscr` | classification | linear discriminant |
| `ecoc` | `fitcecoc` | classification | **multiclass** via error-correcting codes |
| `knn`, `naive_bayes`, `tree_classifier`, `rf_classifier`, `nnet_classifier` | various | classification | non-linear baselines |
| `svr` | `fitrsvm` | regression | kernel SVR |
| `linear_svr` | `fitrlinear` | regression | high-dim linear SVR |
| `lasso` | `fitrlinear` (L1) | regression | sparse |
| `ridge` | `fitrlinear` (L2) | regression | dense, stable |
| `pcr` | PCA + OLS | regression | principal-components regression; reproduces legacy `cv_pcr` / default `cv_lassopcr` |
| `lassopcr` | PCA + LASSO + relaxed-OLS | regression | shrinkage by `lasso_num` (path step) or `estimateparam` (nested-CV lambda), then OLS-refit on the non-zero components; reproduces legacy `cv_lassopcr` shrinkage |
| `gp` | `fitrgp` | regression | Gaussian process; needs few features |
| `tree_regressor`, `rf_regressor`, `nnet_regressor` | various | regression | non-linear baselines |

## 1b. What each algorithm is (and the data it suits)

A short field guide. For whole-brain fMRI (*p ≫ n* — far more voxels than
images), the **linear** models are the default; the nonlinear ones need many
observations relative to features, so reach for them only on ROI/parcel
summaries or after PCA reduction.

**Classification**

- **`svm` (`fitcsvm`)** — the max-margin classifier: finds the hyperplane that
  separates the two classes with the widest gap. Solves the *dual* problem and
  forms an *n × n* kernel matrix, so it supports **nonlinear kernels** (RBF,
  polynomial) but scales with the **number of observations**. Best for
  small/moderate *n* and modest feature counts, or when you want a nonlinear
  boundary. The workhorse for ROI-scale data.
- **`linear_svm` (`fitclinear`)** — the *same* max-margin idea, but a primal
  large-scale solver built for **high-dimensional** data; **linear only**.
  Scales with the **number of features**, so it stays fast at hundreds of
  thousands of voxels where `fitcsvm` bogs down. *svm vs linear_svm*: identical
  objective on linear, different solver and kernel support — use `linear_svm`
  once you're past a few thousand features (whole-brain), `svm` for smaller
  feature counts or a nonlinear kernel.
- **`logistic` (`fitclinear`, logistic loss)** — like `linear_svm` but fits a
  probability model; its scores are (near-)calibrated **probabilities** rather
  than signed distances. Use when you want `P(class)` directly.
- **`lda` (`fitcdiscr`)** — linear discriminant analysis: models each class as a
  Gaussian with a **shared covariance** and draws the optimal linear boundary.
  Fast and well-behaved when that assumption roughly holds; **natively
  multiclass** (no binary reduction needed). Needs regularization
  (`Gamma`/`Delta`) when *p > n*.
- **`knn` (`fitcknn`)** — labels a point by majority vote of its *k* nearest
  neighbours. No training, purely *local*, no weight map; suffers in
  high-dimensional space ("curse of dimensionality"). For low-dimensional /
  reduced data.
- **`naive_bayes` (`fitcnb`)** — assumes features are conditionally independent
  given the class. Cheap and surprisingly strong as a baseline, but the
  independence assumption is poor for correlated voxels.
- **`tree_classifier` / `rf_classifier` / `nnet_classifier`** — a single decision
  tree, a bagged-tree ensemble (random forest), and a feed-forward neural net.
  These capture **nonlinear** structure and interactions, at the cost of needing
  more data and giving a less interpretable map. Nonlinear baselines; rarely the
  first choice on raw voxels.
- **`ecoc` (`fitcecoc`)** — the **multiclass** wrapper (>2 classes); see §4.

**Regression**

- **`svr` (`fitrsvm`)** — support-vector regression: the SVM idea for a
  continuous target (fit within an ε-insensitive tube). Kernel-capable; scales
  with *n*. Often the best **ranker** on brain data (tracks the outcome well).
- **`linear_svr` (`fitrlinear`)** — high-dimensional linear SVR; the wide-data
  counterpart of `svr` (linear only, scales with features).
- **`lasso` (`fitrlinear`, L1)** — linear regression with an **L1 penalty** that
  drives many coefficients to exactly zero → a **sparse**, few-voxel model. Use
  when you expect a compact predictive set and want built-in feature selection.
- **`ridge` (`fitrlinear`, L2)** — linear regression with an **L2 penalty** that
  shrinks coefficients smoothly → a **dense, stable** map. The robust default
  when predictive signal is spread across many correlated voxels (typical for
  fMRI).
- **`pcr`** — principal-components regression: PCA-reduce, then OLS on the
  components. A classic, well-conditioned linear decoder for *p ≫ n*; reproduces
  CANlab's legacy `cv_pcr`.
- **`lassopcr`** — PCR with a LASSO selection of components and a relaxed-OLS
  refit; adds sparsity at the component level (tune via `lasso_num` or nested-CV
  `estimateparam`). The CANlab standard for many published signatures.
- **`gp` (`fitrgp`)** — Gaussian-process regression: a flexible nonparametric
  model with built-in uncertainty, but it forms an *n × n* kernel and **does not
  scale to ~200k voxels** — use inside a `pca` `@pipeline` (§5).
- **`tree_regressor` / `rf_regressor` / `nnet_regressor`** — nonlinear regression
  baselines, analogous to their classifier cousins.

*lasso vs ridge vs pcr/lassopcr*: all four are penalized/projected **linear**
regressions; they differ in how they cope with correlated, high-dimensional
predictors. `ridge` shrinks everything (dense), `lasso` selects a sparse subset,
`pcr` projects onto top variance directions, `lassopcr` selects among those
projections. On smooth, distributed fMRI signal, the dense ones (`ridge`, `pcr`)
are often the most stable; `lasso` wins when the truth really is sparse.

## 2. Setup

```matlab
hw = load_image_set('DPSP_hotwarm');         % binary: Hot (+1) vs Warm (-1)
rf = load_image_set('DPSP_rejectorfriend');

X  = double(hw.dat');  Y  = hw.Y;  id = hw.metadata_table.subj_id;
```

## 3. Binary classification — compare algorithms on identical folds

Cross-validate several classifiers over the **same** splitter so the
comparison is fold-matched (`predict_test_suite` does exactly this for
you, but here it's spelled out):

```matlab
cv   = cv_splitter.stratified_group_kfold(5);
algs = {'svm','linear_svm','logistic'};   % see note on lda below

for a = algs
    pm = crossval(predictive_model('algorithm',a{1},'task','classification'), ...
                  X, Y, 'groups', id, 'cv', cv, 'scoring','balanced_accuracy');
    fprintf('%-12s cv bal-acc = %.3f\n', a{1}, pm.error_metrics.balanced_accuracy.value);
end
```

> **Why no `lda` here:** `fitcdiscr` estimates a *p × p* class covariance, which
> is singular and intractable at ~195k raw voxels. LDA is a fine linear
> classifier, but on whole-brain data it needs heavy regularization
> (`Gamma`/`Delta`) or PCA reduction first (e.g. an `lda` estimator inside a
> `pca` `@pipeline`).

Or in one call:

```matlab
[cverr, yhat, pm_all, results] = predict_test_suite(hw, 'nfolds', id);
disp(results);    % table: algorithm, cv_score, cv_error, n_fdr_sig
```

## 4. Multiclass classification with ECOC

SVMs are intrinsically **binary** — they find one hyperplane between two
classes. To handle *K > 2* classes you either use a learner that is *natively*
multiclass (`lda`, `naive_bayes`, trees/forests, kNN, neural nets all handle
>2 classes directly) or **reduce** the K-class problem to a set of binary
problems. There are three standard reductions:

- **One-vs-all (one-vs-rest, OvA/OvR).** Train *K* binary classifiers, each
  "class *k* vs all the others." At test time, pick the class whose classifier
  is most confident. Cheap (*K* models) but the binary problems are imbalanced
  (1 class vs K−1) and the scores aren't always comparable across classifiers.
- **One-vs-one (OvO, all-pairs).** Train one classifier per **pair** of classes
  — *K(K−1)/2* of them — each trained only on that pair's data; classify by
  majority vote. More models, but each is a small, balanced, easier problem.
  This is **MATLAB's `fitcecoc` default**.
- **Error-correcting output codes (ECOC).** The general framework that contains
  both. Assign each class a **codeword** — a row in a coding matrix of
  `+1 / −1 / 0` ("ignore") — and train one binary learner **per column** (each
  column defines which classes are positive, negative, or left out). At test
  time the *K* learners produce a predicted codeword; assign the class whose row
  is **closest** (Hamming / loss-weighted distance). If the codewords are
  designed with redundancy (well-separated rows), a few wrong binary learners
  still decode to the right class — *error-correcting*, exactly like codes in
  communications. OvA and OvO are just two particular coding matrices; richer
  designs (dense/sparse random, BCH) trade more learners for more robustness.

`ecoc` → `fitcecoc` uses one-vs-one with SVM learners by default; pass
`{'Coding','onevsall'}` (or a custom design / different binary `Learners`) in
`modeloptions` to change it. **For linear SVM specifically you must pick a
reduction** (OvA/OvO/ECOC); if you'd rather avoid that, `lda` gives a native
multiclass linear decoder in one shot.

Build a 4-class problem by stacking the two DPSP tasks — Hot, Warm, Rejecter,
Friend:

```matlab
Xall  = [double(hw.dat'); double(rf.dat')];
Ycls  = [ (hw.Y==1)*1 + (hw.Y==-1)*2 ; (rf.Y==1)*3 + (rf.Y==-1)*4 ];  % 1..4
idall = [hw.metadata_table.subj_id; rf.metadata_table.subj_id];

pm = predictive_model('algorithm','ecoc','task','classification');
pm = crossval(pm, Xall, Ycls, 'groups', idall, ...
              'cv', cv_splitter.group_kfold(5), 'scoring','accuracy');

fprintf('ECOC 4-class accuracy = %.1f%% (chance = 25%%)\n', ...
        100 * pm.error_metrics.accuracy.value);     % ~70%
confusionchart(pm);                                  % see which pairs confuse
```

On this dataset the 4-way accuracy is ~**70%** (chance 25%). The
confusion chart is the interesting part: Hot/Warm and Rejecter/Friend
each confuse *within* their own task far more than across, echoing the
Part 4 finding that the two domains are partly separable.

![ECOC 4-class confusion matrix](pngs/ecoc_confusion.png)

> **Note on scoring:** for multiclass use `accuracy` or
> `balanced_accuracy`. `roc_auc` and the within-person forced-choice
> stats are binary-only and are skipped automatically.

## 5. Regression — compare SVR / lasso / ridge / GP

DPSP ships **binary** labels, so here we treat the signed label as a
continuous target purely to *compare regression algorithms on the same
signal* — in a real study you would swap in a continuous outcome
(temperature, pain rating, …): `Y = hw.metadata_table.rating;`. The
mechanics are identical.

```matlab
cv = cv_splitter.group_kfold(5);
for a = {'svr','linear_svr','lasso','ridge'}
    pm = crossval(predictive_model('algorithm',a{1},'task','regression'), ...
                  X, Y, 'groups', id, 'cv', cv);
    fprintf('%-12s r = %.3f   R^2 = %.3f\n', a{1}, ...
        predictive_model.metric_value(pm.error_metrics,'prediction_outcome_r'), ...
        pm.error_metrics.r2.value);
end
```

Typical output (prediction–outcome correlation `r`, cross-validated
`R^2`):

| algorithm    | r     | R²     |
|--------------|-------|--------|
| `svr`        | 0.44  | 0.09   |
| `linear_svr` | 0.22  | −0.38  |
| `lasso`      | 0.18  | −0.67  |
| `ridge`      | 0.26  | −0.22  |

![Regression algorithm comparison](pngs/algorithm_comparison.png)

`r` (does the prediction *track* the outcome) and `R²` (does it track on
the *right scale*) can disagree — `svr` ranks well **and** is scaled
sensibly here; the penalized linear models rank weakly and a negative R²
means "worse than predicting the mean." This is exactly why you compare
several estimators rather than trusting one.

### Gaussian process needs dimensionality reduction

`fitrgp` builds an n×n kernel and does not scale to ~200k voxels.
Reduce first with a PCA step in an `@pipeline` (which refits the PCA
**per fold** — no leakage):

```matlab
est  = predictive_model('algorithm','gp','task','regression');
pipe = pipeline({ {'pca','k',30} }, est);          % PCA(30) -> GP
pipe = crossval(pipe, X, Y, 'groups', id, 'cv', cv);
fprintf('gp(PCA30)  R^2 = %.3f\n', pipe.error_metrics.r2.value);   % ~0.11
```

The same PCA-then-estimator pattern turns *any* sample-hungry or
feature-hungry learner (GP, kNN, kernel SVR) into something tractable on
wide neuroimaging data, and `weight_map_object(pipe, hw)` still back-projects
the weights to voxel space for linear estimators.

## 6. Hyperparameter tuning with `grid_search`

`grid_search` cross-validates the Cartesian product of a parameter grid
and refits at the winner. Because a tuned `predictive_model` is *still a
predictive_model*, you can wrap it in an outer `crossval` for proper
nested CV.

```matlab
g.BoxConstraint = [0.01 1 100];
pm = predictive_model('algorithm','svm','task','classification', ...
                      'cv', cv_splitter.stratified_group_kfold(5));
pm = grid_search(pm, X, Y, g, 'groups', id);

pm.diagnostics.grid_search.best_params       % e.g. {'BoxConstraint', 1}
pm.diagnostics.grid_search.best_score        % e.g. 0.805
pm.diagnostics.grid_search.scores            % full grid
```

On DPSP Hot/Warm the SVM is robust across two orders of magnitude of
`BoxConstraint`, peaking at the default `1` (~0.80 balanced accuracy) —
a useful sanity check that you're not leaving performance on the table,
and that the model isn't pathologically sensitive to C.

## 7. High-dimensional inference with `stability_selection`

For wide, regularised models the bootstrap z/p collapses (Part 3 §7).
`stability_selection` is the robust alternative: how often does each
voxel land in the top-k by `|weight|` across resamples?

```matlab
pm = stability_selection( ...
        predictive_model('algorithm','linear_svm','task','classification'), ...
        X, Y, 'nboot', 200, 'k', 2000, 'threshold', 0.6, 'groups', id);

ss = pm.diagnostics.stability_selection;
fprintf('%d stable voxels (selected in >= 60%% of bootstraps)\n', ss.n_stable);

% selection-frequency brain map
freq = hw; freq.dat = ss.selection_freq;
montage(freq);
```

Stable voxels are the ones the classifier relies on *regardless* of
which subjects it sees — a far more honest "where is the signal" map than
a bootstrap-z threshold on an L2 model.

## 8. When to use which — a decision guide

| Situation | Reach for |
|---|---|
| 2 classes, want a weight map | `svm` (kernel, default) or `linear_svm` for very wide data |
| 2 classes, want probabilities | `logistic`, or `svm` + `calibrate` |
| > 2 classes | `ecoc` |
| Continuous outcome, linear, interpretable map | `svr` (often the best ranker here), or `ridge` for a dense stable map |
| Continuous outcome, want sparsity | `lasso` |
| Small n, non-linear, few features (or PCA-reduced) | `gp` in a `pca` `@pipeline` |
| Hyperparameters matter | `grid_search` (wrap in outer `crossval` for nested CV) |
| Voxel-level inference on a wide regularised model | `stability_selection` (not bootstrap z/p) |
| Honest CV with any preprocessing (PCA, scaling, feature select) | `@pipeline` (refits every step per fold) |

### Rules of thumb

- **Start simple.** Linear `svm`/`svr` with default settings is a strong
  baseline and gives an interpretable weight map; only add capacity (RBF
  kernels, GP, ensembles) if a fold-matched comparison says it helps.
- **Compare on identical folds.** Always cross-validate competing
  algorithms over the *same* `cv_splitter` (or use `predict_test_suite`),
  or fold-to-fold noise will swamp the algorithm difference.
- **Match the metric to the question.** Ranking → `roc_auc` / `r` /
  forced-choice; calibrated scale → `R²` / accuracy at threshold.
- **Group by subject.** With multiple images per subject, always pass
  `'groups', subject_id` and a group-aware splitter, or held-out folds
  leak subject identity and inflate performance.

That completes the five-part multivariate-decoding walkthrough:
**Part 1** SVM classification basics, **Part 2** classification vs
regression, **Part 3** the composable sklearn-style `predictive_model`
API, **Part 4** cross-classification, and **Part 5** algorithms, tuning,
and inference.
