%% Mass-univariate regression with the glm_map object
%
% This walkthrough shows two equivalent ways to run a voxelwise multiple
% regression in CanlabCore and work with the result as a |glm_map| object:
%
%   (1) the quick path  - call fmri_data.regress directly; it now *returns*
%                         a glm_map, which you then access and query.
%   (2) the estimator path - build a glm_map, attach data, and call fit(),
%                         add_contrasts(), diagnostics(), threshold(), etc.
%                         (a scikit-learn-style API).
%
% It also covers first-level (event) designs, and how the object mirrors the
% fields of the regression results structure (the variable "out" inside
% fmri_data.regress), with related options grouped into nested structs
% (.input_parameters, .input_image_metadata, .diagnostics).
%
% Requirements: CanlabCore + Neuroimaging_Pattern_Masks + SPM12 on the path.
% Every block below is copy-pasteable.
%
% See also: fmri_data.regress, glm_map, statistic_image, fmri_glm_design_matrix


%% 1. Load sample data and define a design
% load_image_set('emotionreg') returns an fmri_data object with 30 first-level
% contrast images (one per participant). For a 2nd-level (group) regression we
% put a group design in a [images x regressors] matrix X: an intercept plus one
% covariate of interest (here, a stand-in continuous predictor).

dat = load_image_set('emotionreg', 'noverbose');   % 30 images

n = size(dat.dat, 2);                               % number of images
covariate = zscore((1:n)');                         % stand-in continuous predictor
X = [covariate ones(n, 1)];                         % regressor + intercept


%% 2. Quick path: fmri_data.regress returns a glm_map
% Put the design in dat.X and call regress. The result "g" is a glm_map object.
% (Historically regress returned a plain struct; it now re-casts that struct as
% a glm_map. All the old struct-style field names still work - see step 4.)

dat.X = X;
g = regress(dat, 0.001, 'unc', 'names', {'covariate'}, 'noverbose', 'nodisplay');

g                                                   % typing the name lists ALL properties
% disp(g) shows: analysis_name, design, contrasts, input_parameters,
% input_image_metadata, betas, t, contrast_estimates, contrast_t, df, sigma,
% residuals, dfe, diagnostics, warnings, ... and expands the nested structs.


%% 3. Access and query the fitted maps
% Result maps are statistic_image / fmri_data objects living inside the glm_map.

g.betas            % statistic_image: one beta image per regressor [voxels x 2]
g.t                % statistic_image: thresholded t map per regressor
g.df               % fmri_data: per-voxel error degrees of freedom
g.sigma            % fmri_data: per-voxel residual standard deviation
g.dfe              % scalar error df (median of g.df), convenience summary

% Design diagnostics are collected in the nested .diagnostics struct. The VIF
% and leverage fields keep the names used by fmri_data.regress out.diagnostics.
g.diagnostics.Variance_inflation_factors            % VIF per regressor
g.diagnostics.Contrast_variance_inflation_factors   % contrast VIFs (cVIF), if contrasts present
g.diagnostics.condition_number                      % conditioning of X
g.diagnostics.Leverages                             % per-observation leverage
g.diagnostics.Cooks_distance                        % per-observation influence (when residuals available)

% Options the fit actually used are recorded in .input_parameters, and the
% provenance of the input images in .input_image_metadata
g.input_parameters
g.input_image_metadata

% Visualize / tabulate a chosen map (delegates to statistic_image methods)
montage(g.t);                       % montage of the thresholded t map
% table(g.t);                       % atlas-labeled table of significant regions
% r = region(g.t);                  % region object for further ROI work


%% 4. The object mirrors the regress out-struct (with back-compatible aliases)
% glm_map property names match the field names of "out" in fmri_data.regress.
% Where the historical field name differs from the canonical property, an alias
% reads/writes the same data, so legacy struct-style access is unchanged:
%
%   out.b               -> g.betas
%   out.contrast_images -> g.contrast_estimates
%   out.con_t           -> g.contrast_t
%   out.resid           -> g.residuals
%   out.variable_names  -> g.regressor_names
%   out.C               -> g.contrasts

isequal(g.b.dat, g.betas.dat)            % true: .b is an alias for .betas
isequal(g.variable_names, g.regressor_names)   % true

% You can also re-cast any regress-style struct into a glm_map yourself:
%   g2 = glm_map(out_struct);


%% 5. Estimator path: build, attach data, fit
% The same analysis as a reusable estimator. Construct the design first (no
% data), screen it, add contrasts, then fit. This separates design choices from
% the data and lets you inspect collinearity *before* spending compute.

g = glm_map('X', X, 'level', 2, 'regressor_names', {'covariate', 'intercept'});

g = add_contrasts(g, [1 0], {'covariate_effect'});   % one row per contrast

g = diagnostics(g);                                  % VIF / leverage / conditioning (no fit needed)

g = fit(g, dat);                                     % runs fmri_data.regress under the hood

g.is_fitted                                          % true
summary(g);                                          % narrative report: model, diagnostics, results
montage(g, 'contrast_t');                            % thresholded contrast t map
% table(g, 'contrast');                              % table for the contrast map
%
% summary(g) prints the analysis name, the model (level + input variables),
% design diagnostics (max VIF, contrast VIF, condition number, leverage, max
% Cook's distance), and -- once fitted -- how many maps, the statistical
% threshold, and the number of significant voxels per regressor and contrast.


%% 6. Re-threshold without refitting
% threshold() re-masks the stored t / contrast_t maps; the underlying statistic
% values are preserved, so you can sweep thresholds cheaply.

g = threshold(g, .005, 'unc', 'k', 10);              % both t and contrast_t
g = threshold(g, .05, 'fdr', 'which_map', 'contrast');  % contrast map only


%% 7. First-level (event) designs
% For within-run BOLD timeseries, wrap an fmri_glm_design_matrix in .design.
% glm_map builds X by HRF-convolving the onsets, and marking is_timeseries lets
% fit() use autoregressive (AR) error models.

  TR = 2; nscan = 200;
  onsets = {[10 40 70 100]', [25 55 85 115]'};     % seconds, two conditions
  d = fmri_glm_design_matrix(TR, 'nscan', nscan, 'units', 'secs', ...
                             'onsets', onsets, 'condition_names', {'A','B'});
  g = glm_map(d);                                  % level 1, event mode
  g.is_timeseries = true;
  g = build_design(g);                             % onsets -> X via convolution
  g = fit(g, bold_timeseries_fmri_data, 'AR', 1);  % AR(1) error model

% Or import a full SPM12/SPM25 first-level model:

  g = import_SPM(glm_map, '/path/to/SPM.mat');


%% 8. Summary
% - fmri_data.regress returns a glm_map; query it via .betas/.t/.df/.sigma,
%   the nested .diagnostics / .input_parameters / .input_image_metadata structs,
%   or the historical out-struct aliases (.b/.con_t/.contrast_images/...).
% - Or build a glm_map as an estimator: glm_map(...) -> add_contrasts ->
%   diagnostics -> fit -> threshold/table/montage.
% - typing the object name (or summary(g)) lists all properties; methods(glm_map)
%   lists all operations.
