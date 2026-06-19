%% Mass-univariate regression with the glm_map object
%
% This walkthrough shows two equivalent ways to run a voxelwise multiple
% regression in CanlabCore and work with the result as a |glm_map| object:
%
%   (1) the quick path  - call fmri_data.regress directly; it now *returns*
%                         a glm_map, which you then access, display, and query.
%   (2) the estimator path - build a glm_map, attach data, and call fit(),
%                         add_contrasts(), diagnostics(), threshold(), etc.
%                         (a scikit-learn-style API).
%
% It mirrors the worked examples in the fmri_data.regress help (naming,
% nuisance covariates, gray-matter masking, outlier handling, robust fits,
% residuals, brain-predicts-behavior, writing maps, diagnostics) and shows
% how the object groups related outputs into nested structs
% (.input_parameters, .input_image_metadata, .diagnostics).
%
% Requirements: CanlabCore + Neuroimaging_Pattern_Masks + SPM12 on the path.
% Every block below is copy-pasteable.
%
% See also: fmri_data.regress, glm_map, statistic_image, fmri_glm_design_matrix


%% 1. Load sample data and define a real predictor from the metadata table
% load_image_set('emotionreg') returns an fmri_data object with 30 first-level
% contrast images (one per participant; Wager et al. 2008, Neuron) plus a
% metadata_table. We predict each participant's brain map from their
% behavioral reappraisal success score (a 2nd-level / group regression).

obj = load_image_set('emotionreg', 'noverbose');   % 30 images
n   = size(obj.dat, 2);                             % number of images

% Predictor matrix from metadata_table. A single column here; the intercept is
% added automatically if absent. Mean-centering the regressor(s) makes the
% intercept map interpretable as the activation for the average participant.
obj.X = obj.metadata_table.Reappraisal_Success;
obj.X = obj.X - mean(obj.X);                        % mean-center
obj.X(:, end + 1) = 1;                              % regressor + intercept

X = obj.X;                                          % keep for later sections


%% 2. Quick path: fmri_data.regress returns a glm_map
% Run the regression at a liberal threshold. The result "g" is a glm_map
% object. (Historically regress returned a plain struct; it now re-casts that
% struct as a glm_map. All the old struct-style field names still work - see
% section 6.)

g = regress(obj, .05, 'unc', 'noverbose', 'nodisplay');

g                                                   % typing the name lists ALL properties

% The thresholded t map (one statistic_image per regressor, incl. intercept):
g.t


%% 3. Name the regressors and the analysis
% Naming labels the output maps (.image_labels), which propagates to montages
% and tables; analysis_name is carried on the object and printed by summary.

g = regress(obj, .05, 'unc', 'names', {'Reapp_Success' 'Intercept'}, ...
            'analysis_name', 'Emotion Regulation', 'noverbose', 'nodisplay');

g.t.image_labels                                    % {'Reapp_Success','Intercept'}


%% 4. Display the result maps
% There is one image per regressor in X. Select, visualize, and re-threshold
% them with the statistic_image display methods.

montage(g.t);                                       % all regressor t maps

t_reapp = get_wh_image(g.t, 1);                     % just the Reapp_Success map
create_figure('surface'); surface(t_reapp);

t_reapp = threshold(t_reapp, .005, 'unc');          % re-threshold a single map
create_figure('montage P<.005'); axis off; montage(t_reapp);

orthviews(g.t);                                     % orthviews of the thresholded maps

% glm_map also has object-level display wrappers that pick a map by name:
montage(g, 't');
% table(g, 't');                                    % atlas-labeled results table


%% 5. Query the object and plot design diagnostics
% Result maps are statistic_image / fmri_data objects inside the glm_map.

g.betas            % statistic_image: one beta image per regressor
g.df               % fmri_data: per-voxel error degrees of freedom
g.sigma            % fmri_data: per-voxel residual standard deviation
g.dfe              % scalar error df (median of g.df)

% Design diagnostics are collected in the nested .diagnostics struct, using
% the same field names as fmri_data.regress out.diagnostics.
g = diagnostics(g);                                 % Run diagnostics and return output in object
g.diagnostics.Variance_inflation_factors            % VIF per regressor
g.diagnostics.Contrast_variance_inflation_factors   % contrast VIFs (cVIF), if contrasts present
g.diagnostics.Leverages                             % per-observation leverage
g.diagnostics.condition_number                      % conditioning of X
g.diagnostics.Cooks_distance                        % per-observation influence

% Options the fit used live in .input_parameters; input-image provenance in
% .input_image_metadata.
g.input_parameters
g.input_image_metadata

% Plot the design diagnostics
figure;
subplot(1, 2, 1); plot(g.diagnostics.Variance_inflation_factors, 'o-'); title('VIFs'); xlabel('Regressor');
subplot(1, 2, 2); plot(g.diagnostics.Leverages, 'o-'); title('Leverage'); xlabel('Observation');


%% 6. The object mirrors the regress out-struct (with back-compatible aliases)
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

isequal(g.b.dat, g.betas.dat)                   % true: .b is an alias for .betas
isequal(g.variable_names, g.regressor_names)    % true

% You can also re-cast any regress-style struct into a glm_map yourself:
%   g2 = glm_map(out_struct);


%% 7. Add a 2nd-level nuisance covariate (mean CSF)
% Covariates for effects of no interest stabilize results. Mean gray, white,
% and CSF signals are highly correlated; CSF likely carries image-wide
% confounds rather than task signal, so it is a reasonable nuisance covariate.

gwcsf = extract_gray_white_csf(obj);                % [n x 3]: gray, white, CSF means

obj.X = [obj.metadata_table.Reappraisal_Success, gwcsf(:, 3)];
obj.X = obj.X - mean(obj.X);                        % mean-center predictors
obj.X(:, end + 1) = 1;                              % intercept

g = regress(obj, .05, 'unc', 'names', {'Reapp_Success' 'CSF_mean' 'Intercept'}, ...
            'analysis_name', 'Emotion Regulation', 'noverbose', 'nodisplay');
montage(g.t);


%% 8. Apply a gray-matter mask before thresholding (FDR)
% Correcting within hypothesized gray-matter voxels reduces the comparison
% count and can make results more significant.

gm = fmri_data(which('gray_matter_mask.nii'));
g.t = apply_mask(g.t, gm);
g.t = threshold(g.t, 0.05, 'fdr');
montage(g.t);


%% 9. Find and exclude outliers, then re-run
% 'notimeseries' omits time-series-specific calculations (this is 2nd-level
% data). We flag images that do not correlate with the others (mahal_corr).

[~, wh_outliers] = outliers(obj, 'notimeseries');
find(wh_outliers)                                   % indices flagged as outliers

obj_clean = get_wh_image(obj, ~wh_outliers);        % drop flagged images
obj_clean.X = obj.X(~wh_outliers, :);               % and the matching design rows

g_clean = regress(obj_clean, .05, 'unc', ...
            'names', {'Reapp_Success' 'CSF_mean' 'Intercept'}, 'noverbose', 'nodisplay');
g_clean.t = threshold(apply_mask(g_clean.t, gm), 0.01, 'unc');
montage(g_clean.t);


%% 10. More regress options
% Reset to the single-predictor design for these.
obj.X = X;

% Robust regression (iteratively down-weights outliers; slower than OLS):
g_robust = regress(obj, .05, 'fdr', 'robust', 'names', {'Reapp_Success' 'Intercept'}, ...
            'noverbose', 'nodisplay');

% Save voxelwise residuals as an fmri_data object (useful for denoising / QC):
g_resid = regress(obj, .001, 'unc', 'residual', 'noverbose', 'nodisplay');
g_resid.residuals                                   % fmri_data [voxels x images]

% Brain-predicts-behavior ('brainony'): a univariate map of where brain
% activity predicts obj.Y. Needs obj.Y and is slow (loops over voxels):
%   obj.Y = obj.metadata_table.Reappraisal_Success;
%   g_bxy = regress(obj, .05, 'unc', 'brainony');

% Re-threshold and re-display an existing map without refitting:
g.t = threshold(g.t, .001, 'unc');
orthviews(g.t);

% Write a beta image to disk (uncomment to run; writes to the current folder):
%   g.betas.fullpath = fullfile(pwd, 'beta_reapp_success.nii');
%   write(g.betas);


%% 11. Estimator path: build, screen, fit
% The same analysis as a reusable estimator. Construct the design first (no
% data), screen its collinearity, add contrasts, then fit. fit() calls
% fmri_data.regress under the hood and stores the result maps on the object.

g = glm_map('X', X, 'level', 2, 'regressor_names', {'Reapp_Success' 'Intercept'}, ...
            'analysis_name', 'Emotion Regulation (estimator)');

g = add_contrasts(g, [1 0], {'reapp_effect'});      % one row per contrast

g = diagnostics(g);                                 % VIF / cVIF / leverage / conditioning (no fit)

g = fit(g, obj);                                    % runs the regression

g.is_fitted                                         % true
summary(g);                                         % narrative report: model, diagnostics, results
montage(g, 'contrast_t');                           % thresholded contrast t map
% table(g, 'contrast');                             % table for the contrast map
%
% summary(g) prints the analysis name, the model (level + input variables),
% design diagnostics (max VIF, contrast VIF, condition number, leverage, max
% Cook's distance), and -- once fitted -- how many maps, the statistical
% threshold, and the number of significant voxels per regressor and contrast.


%% 12. Re-threshold without refitting
% threshold() re-masks the stored t / contrast_t maps; the underlying
% statistic values are preserved, so you can sweep thresholds cheaply.

g = threshold(g, .005, 'unc', 'k', 10);             % both t and contrast_t
g = threshold(g, .05, 'fdr', 'which_map', 'contrast');  % contrast map only


%% 13. First-level (event) designs and AR models
% For within-run BOLD timeseries, wrap an fmri_glm_design_matrix in .design;
% glm_map builds X by HRF-convolving the onsets. Marking is_timeseries lets
% fit() use autoregressive (AR) error models. build_design runs here; the fit
% and import_SPM lines need real timeseries / an SPM.mat, so they are shown
% for reference.

TR = 2; nscan = 200;
onsets = {[10 40 70 100]', [25 55 85 115]'};        % seconds, two conditions
d = fmri_glm_design_matrix(TR, 'nscan', nscan, 'units', 'secs', ...
                           'onsets', onsets, 'condition_names', {'A' 'B'});
g_evt = glm_map(d);                                 % level 1, event mode
g_evt.is_timeseries = true;
g_evt = build_design(g_evt);                        % onsets -> X via convolution
plot_design(g_evt);                                 % design matrix + VIFs

%   g_evt = fit(g_evt, bold_timeseries_fmri_data, 'AR', 1);  % AR(1) error model
%   g_evt = import_SPM(glm_map, '/path/to/SPM.mat');         % import a full 1st-level model


%% 14. Summary
% - fmri_data.regress returns a glm_map; access/display it via .betas/.t/
%   montage/surface/orthviews, query the nested .diagnostics /
%   .input_parameters / .input_image_metadata structs, or use the historical
%   out-struct aliases (.b/.con_t/.contrast_images/...).
% - Or build a glm_map as an estimator: glm_map(...) -> add_contrasts ->
%   diagnostics -> fit -> summary/threshold/table/montage.
% - typing the object name lists all properties; summary(g) prints a narrative
%   report; methods(glm_map) lists all operations.
