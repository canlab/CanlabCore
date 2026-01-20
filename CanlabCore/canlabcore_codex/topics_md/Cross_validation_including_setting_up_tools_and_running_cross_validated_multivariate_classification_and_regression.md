# Cross_validation_including_setting_up_tools_and_running_cross_validated_multivariate_classification_and_regression

This topic covers the main cross-validation entry points for multivariate prediction in CANlabCore, plus supporting utilities for cross-validated regression, bootstrap weight maps, and connectivity-focused prediction workflows.

## Object methods

### fmri_data.predict
Primary entry point for cross-validated prediction from brain data, supporting multiple algorithms and error metrics.
Major options: 'nfolds' (or custom fold vector), 'error_type' (mcr/mse), 'algorithm_name' (e.g., cv_pcr, cv_lassopcr, cv_svm), 'useparallel', 'bootweights', 'bootsamples', 'numcomponents', 'lasso_num', 'hvblock', 'rolling', 'platt_scaling', 'verbose'.

### fmri_data.predict_test_suite
Runs a panel of cross-validated algorithms and plots predicted vs. observed outcomes.
Major options: 'quick' (skip extended output), 'nfolds' (fold count or custom holdout vector), passes options through to fmri_data.predict.

### fmri_data.model_mpathi
Cross-validated modeling of latent pathway connectivity with default k-fold CV and optional custom indices.
Major options: custom fold indices, output plotting of cross-validated correlations, and pathway definitions (see function header).

### fmri_data.model_brain_pathway
Cross-validated pathway modeling with default 10-fold CV for testing connection strength and latent pathway models.
Major options: custom CV indices, pathway definitions, and output controls (see function header).

## Stand-alone functions

### xval_SVM (Statistics_tools/Cross_validated_Regression/xval_SVM.m)
Cross-validated SVM classification with optional nested optimization, repeats, and bootstrapping.
Major options: 'nooptimize', 'norepeats'/'dorepeats', 'nboot'/'nobootstrap', 'noverbose', 'noplot', 'highdimensional', 'modeloptions', and custom subject IDs for leave-whole-subject-out CV.

### xval_SVR (Statistics_tools/Cross_validated_Regression/xval_SVR.m)
Cross-validated support vector regression with optional nested optimization, repeats, and bootstrapping.
Major options: 'nooptimize', 'norepeats'/'dorepeats', 'nboot'/'nobootstrap', 'noverbose', 'noplot', and subject IDs for leave-whole-subject-out CV.

### xval_classify (Statistics_tools/xval_classify.m)
Cross-validated discriminant classification using fitcdiscr with balanced folds.
Major options: holdout strategy (leave-whole-subject-out via xval_stratified_holdout_leave_whole_subject_out or balanced folds via stratified_holdout_set), nfolds, and plotting controls (see function header).

### cv_mlpcr, cv_mlpcr_bt, cv_mlpcr_wi (mlpcr/)
Cross-validated multilevel PCR wrappers; bootstrapped and within-subject variants are provided.
Major options: passed through to mlpcr2, including 'subjIDs', 'numcomponents' ([between, within]), and 'cpca'.

### xval_regression_multisubject (Statistics_tools/Cross_validated_Regression)
Cross-validated regression for single- or multi-subject data with multiple holdout schemes.
Major options: 'pca', 'ndims', 'variable', 'cov'/'covs', 'lassopath'/'ridgek'/'regparams', 'holdout_method' (loo, balanced4, categorical_covs, l4o_covbalanced, or custom), 'optimize_regularization', 'nested_choose_ndims', verbosity flags.

### xval_regression_multisubject_featureselect
Extension of xval_regression_multisubject that includes feature selection within folds.
Major options: inherits xval_regression_multisubject options plus feature selection controls (see function header).

### xval_regression_multisubject_bootstrapweightmap
Runs cross-validated regression with bootstrap weight map estimation.
Major options: inherits xval_regression_multisubject options plus bootstrap parameters (see function header).

### xval_regression_bootstrapweightmap_brainplots
Utility for visualizing bootstrap weight maps from cross-validated regression outputs.
Major options: plotting and threshold controls (see function header).

### canlab_connectivity_predict
Cross-validated prediction using connectivity features with selectable algorithms and fold counts.
Major options: 'algo' (e.g., cv_svm, cv_svr, cv_lassopcr), 'folds', 'error_type', and outcome/covariate inputs.

### stratified_holdout_set (Statistics_tools/Support_functions/stratified_holdout_set.m)
Generates stratified holdout fold assignments for classification or continuous outcomes.
Major options: nfolds or custom integer fold vector; 'hvblock' and 'rolling' for time-series CV.

### xval_stratified_holdout_leave_whole_subject_out (Statistics_tools/Cross_validated_Regression/xval_stratified_holdout_leave_whole_subject_out.m)
Creates stratified k-fold splits while keeping all observations from a subject/group in the same fold.
Major options: nfolds, 'doverbose', 'doplot', and grouping IDs.

### xval_select_holdout_set (Statistics_tools/Cross_validated_Regression/xval_select_holdout_set.m)
Selects balanced holdout sets by optimizing independence between outcomes/covariates and holdout assignment.
Major options: nfolds, holdout_proportion, covariates, 'doplot'.

### xval_select_holdout_set_categoricalcovs (Statistics_tools/Cross_validated_Regression/xval_select_holdout_set_categoricalcovs.m)
Builds holdout sets balanced across categorical covariate combinations.
Major options: categorical covariate matrix; returns one logical mask per fold.
