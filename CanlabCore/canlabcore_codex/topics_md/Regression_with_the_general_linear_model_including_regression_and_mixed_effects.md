# Regression_with_the_general_linear_model_including_regression_and_mixed_effects

This topic covers voxelwise GLM regression and mixed-effects modeling for neuroimaging data, along with regression utilities for model selection, nonlinear regression, contrast construction, and collinearity diagnostics.

## Object methods

### fmri_data

#### regress
Multiple regression at each voxel using obj.X as the design matrix, returning beta/t maps and optional contrasts.
Major options: 'robust', 'AR', 'grandmeanscale', 'C' (contrast matrix), 'nointercept', 'display'/'nodisplay', 'brainony', 'residual', 'variable_names'/'names', 'contrast_names', 'analysis_name', 'covdat', and threshold type (e.g., .001, 'unc').

#### fitlme_voxelwise
Voxelwise mixed-effects modeling using fitlme; accepts a formula string or a role table and returns fixed-effect maps and contrasts.
Major options: p-threshold and type ('unc'/'fdr'), 'C' (contrast matrix), 'contrast_names', 'analysis_name', 'doparallel', 'fitmethod' (REML/ML), 'residual', 'noverbose'.

#### rsa_regression
Representational similarity regression that fits an RDM to a design matrix and returns bootstrapped inference on design regressors.
Major options: distance metrics ('correlation' default, 'average_Euclidean', 'euclidean', 'cosine'); optional bootstrap control (see function notes).

#### dual_regression
Convenience wrapper for dual regression on fmri_data objects; calls dual_regression_fsl after resampling.
Major options: passed through to dual_regression_fsl (zscore_data, zscore_maps, n_iter, add_intercept, verbose, doplot).

## Stand-alone functions

### Mixed-effects and multilevel regression

#### igls_multicond (Statistics_tools/Iterative_Generalized_Least_Squares/igls_multicond.m)
Iterative generalized least squares for multilevel variance component estimation, with IGLS and RIGLS modes.
Major options: 'covariate', 'noverbose', 'iter'/'iterations', 'type' ('i' for igls, 'r' for rigls), 'eps'/'epsilon', 'ar'/'arorder', 'within_var'/'within', 'plot' ('all'/'slopes'/'design'), 'names'.

#### glmfit_multilevel (Statistics_tools/glmfit_multilevel.m)
Fast two-level mixed-effects GLM with random intercepts and slopes across second-level units (e.g., subjects).
Major options: 'names', 'analysisname', 'beta_names', 'robust', 'weighted'/'weight'/'var'/'s2', 'nocenter', 'noint', inference ('bootstrap'/'signperm'/'ttest'), 'nresample', 'pvals', 'permsign', plotting ('plots'/'noplots'), saving ('save'/'saveplots', 'savefile'/'savefilename'), 'verbose'/'noverbose'.

### Contrast and design utilities

#### create_orthogonal_contrast_set (Statistics_tools/create_orthogonal_contrast_set.m)
Creates an orthogonal contrast set across conditions, with rows that sum to zero and are orthogonal.
Major options: none beyond nconditions input.

### Collinearity diagnostics

#### getvif (OptimizeDesign11/core_functions/getvif.m)
Computes variance inflation factors for a design matrix, optionally plotting and selecting subsets of columns.
Major options: no_add_intercept (flag), 'wh' (column indices), 'plot'.

#### cVIF (Statistics_tools/cVIF.m)
Computes contrast-based variance inflation factors using standardized design columns and provided contrasts.
Major options: none beyond X and contrast matrix inputs; intercept columns are removed automatically.

### Model selection and nonlinear regression

#### regress_best_subsets_ga
Genetic-algorithm search for best-subset regression using AIC as the objective.
Major options: none beyond inputs X and Y (algorithm settings are internal).

#### regression_cubic_bspline
Cubic B-spline regression with optional column expansion for nonlinear effects.
Major options: 'expand_cols', 'names', 'verbose', 'doplot'.

#### monotonic_regression
Isotonic regression (nondecreasing fit) with optional plotting.
Major options: doplot (0/1).

### Dual regression (stand-alone)

#### dual_regression_fsl
FSL-style dual regression for ICA spatial maps and subject data matrices.
Major options: 'zscore_data', 'zscore_maps', 'n_iter', 'add_intercept', 'verbose', 'doplot'.
