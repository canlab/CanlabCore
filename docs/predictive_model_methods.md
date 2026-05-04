# `predictive_model` methods, organized by area

`predictive_model` is a container for the artifacts of a fitted
multivariate predictive model: the outcome `Y`, predictions `yfit`,
weight maps `w`, bootstrap statistics, cross-validation indices and
errors, the underlying MATLAB classification/regression model object,
and effect-size summaries. It is intended as an output type for
`xval_SVM`, `xval_SVR`, and similar wrappers; the constructor accepts a
plain struct (such as the output of those wrappers) and copies matching
fields, with legacy translations for `SVMModel`/`SVRModel` and the
classification/regression `_d_singleinterval` and `_d_within` fields.

Many methods on the class are placeholders ("forthcoming") that error if
called; see the table below for which ones are implemented today. Type
`methods(my_obj)` in MATLAB for the live list on any instance.

## Properties

| Property | Description |
|---|---|
| `Y` | Observed outcomes (vector); ±1 for SVM |
| `id` | Grouping/subject vector for within-participant observations |
| `class_labels` | Cell of class names, e.g. `{'No pain' 'Pain'}` |
| `Y_name` | Name of outcome variable |
| `X_name` | Name of predictor variable |
| `modeloptions` | Cell of fitter options, e.g. `{'KernelFunction','linear'}` |
| `accfun` | Function handle computing accuracy from `(Y, yfit)` |
| `accfun_descrip` | Text description of `accfun` |
| `trIdx` | Cell of training-fold logical indices |
| `teIdx` | Cell of test-fold logical indices |
| `nfolds` | Number of cross-validation folds |
| `dist_from_hyperplane_xval` | Cross-validated distances from the SVM hyperplane |
| `yfit` | Cross-validated predicted outcomes |
| `class_probability_xval` | Cross-validated Platt-scaled class-1 probabilities |
| `crossval_accuracy` | Cross-validated accuracy scalar |
| `crossval_accuracy_descrip` | Text description of the accuracy metric |
| `classification_d_singleinterval` | Single-interval classification effect size |
| `mult_obs_within_person` | Flag for multiple observations per subject |
| `ClassificationModel` | Underlying full-data MATLAB model object (e.g. `ClassificationSVM`) |
| `w` | Model weight vector |
| `boot_w` | Bootstrap weight samples |
| `boot_w_ste` | Bootstrap standard errors of weights |
| `boot_w_mean` | Bootstrap mean of weights |
| `wZ` | Bootstrap Z-scores for weights |
| `wP` | Bootstrap p-values for weights |
| `wP_fdr_thr` | FDR threshold (q < 0.05) for `wP` |
| `boot_w_fdrsig` | Logical vector of FDR-significant weights |
| `w_thresh_fdr` | Weights thresholded at FDR q < 0.05 |
| `cverrfun` | Function handle for cross-validation error |
| `cverr` | Cross-validation error values |
| `prediction_outcome_r` | Pearson r between `yfit` and `Y` |
| `prediction_outcome_r_descrip` | Description of the prediction-outcome correlation |
| `d_singleinterval` | Single-interval effect size (regression or classification) |
| `d_within` | Within-person effect size |

## Basic operations

| Method | From | One-liner |
|---|---|---|
| `predictive_model` | `@predictive_model` | Constructor; copies a struct's matching fields and runs `validate_object` |
| `validate_object` | `@predictive_model` | Type/shape checks for all properties |

## Display and visualization

| Method | From | One-liner |
|---|---|---|
| `plot_predicted_vs_observed` | `@predictive_model` | Scatter (regression) or violin (binary classification) of `yfit` vs. `Y` |

## Tables

| Method | From | One-liner |
|---|---|---|
| `confusion_matrix` | `@predictive_model` | Raw and row-normalized confusion matrices, with optional plotting |

## Statistics (forthcoming stubs)

These methods exist as stubs and currently throw an error when called.

| Method | From | One-liner |
|---|---|---|
| `train` | `@predictive_model` | Train the model (forthcoming) |
| `test` | `@predictive_model` | Test the model on held-out data (forthcoming) |
| `crossval` | `@predictive_model` | Cross-validate the model (forthcoming) |
| `bootstrap` | `@predictive_model` | Bootstrap analysis of weights (forthcoming) |
| `permutation_test` | `@predictive_model` | Permutation test for significance (forthcoming) |
| `select_features` | `@predictive_model` | Feature selection (forthcoming) |
| `error_analysis` | `@predictive_model` | Analyze misclassified images and error patterns (forthcoming) |
| `report` | `@predictive_model` | Generate a textual performance report (forthcoming) |
| `plot` | `@predictive_model` | Default predictions-vs-outcomes plot (forthcoming) |
| `montage` | `@predictive_model` | Display a montage of relevant images (forthcoming) |
| `confusionchart` | `@predictive_model` | Confusion-chart visualization (forthcoming) |
| `rocplot` | `@predictive_model` | ROC curves (forthcoming) |
