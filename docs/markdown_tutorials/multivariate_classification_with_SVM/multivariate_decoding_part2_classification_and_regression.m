%[text] # Multivariate decoding — Part 2: classification and regression
%[text] Part 2 of the multivariate-decoding series (Part 1 classification basics with SVM; Part 3 the predictive\_model API; Part 4 cross-classification; Part 5 algorithms & tuning). See the matching .md for the full narrative and figures.
%[text] Both classification and regression are supervised decoding: learn a map from a brain image **X** to an outcome **Y**, cross-validate, read out a weight map and an out-of-sample score. They differ only in **Y**: categorical (classification → accuracy / ROC / confusion) vs continuous (regression → prediction–outcome *r* / R²). Everything else is shared.
%%
%[text] ## 1. One-line loaders
%[text] Two keyword datasets give a ready-to-decode |fmri_data| object with |.Y| and a |metadata_table| already populated — no manual masking, concatenation, or label-building.
hw_obj = load_image_set('DPSP_hotwarm', 'noverbose');   % classification: Y = +1 (Hot) / -1 (Warm)
bmrk3  = load_image_set('bmrk3', 'noverbose');          % regression: Y = continuous pain rating
%%
%[text] ## 2. Classification end-to-end with fmri_data.predict
%[text] Pass the image object, choose an algorithm by keyword, and with |'newapi'| get back a |predictive_model| (4th output) you can visualize directly. The model carries a |weights.weight\_obj|, so |montage(pm)| / |rocplot(pm)| / |confusionchart(pm)| work with no extra arguments.
rng(2026)
[cverr, stats, optout, pm] = predict(hw_obj, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'newapi');
classification_accuracy = stats.acc
create_figure('weights'); axis off; montage(pm);
create_figure('roc'); ROC = rocplot(pm);
confusionchart(pm);
%%
%[text] ## 3. Regression end-to-end with fmri_data.predict
%[text] The same call with a regression algorithm and continuous Y. We predict pain rating from the heat-evoked maps, using whole-subject folds. For a regression algorithm on a binary outcome the |'newapi'| path auto-uses MSE; here Y is already continuous.
[~, ~, sid] = unique(bmrk3.metadata_table.subject_id, 'stable');
folds = mod(sid, 5) + 1;                                 % whole-subject 5-fold
[cverr_r, stats_r, ~, pm_r] = predict(bmrk3, 'algorithm_name', 'cv_pcr', 'nfolds', folds, 'newapi');
prediction_outcome_r = corr(stats_r.yfit, bmrk3.Y)
%[text] Beyond the correlation, report **predicted R²** (1 - PRESS/SST; Wager & Lindquist Ch. 39.4) — the variance explained by the held-out predictions. |report_accuracy| / |summary| print the model-type-relevant metric block (and, for |summary|, the CV scheme and which inference is available).
predicted_r2 = pm_r.error_metrics.predicted_r2.value
report_accuracy(pm_r);
summary(pm_r);
%%
%[text] ### Prediction–outcome correlation (the headline regression read-out)
create_figure('predicted vs observed');
plot_correlation_samefig(stats_r.yfit, bmrk3.Y);
xlabel('Cross-validated prediction'); ylabel('Observed pain rating'); set(gca, 'FontSize', 14);
%%
%[text] ### Weight maps — unthresholded and bootstrap-thresholded
%[text] The unthresholded map shows the full pattern; bootstrapping (resampling whole subjects) and thresholding shows which voxels are reliable. Note the empirical bootstrap p floors at 2/(nboot+1), so an FDR threshold over ~220k voxels can be empty — threshold uncorrected, use more bootstraps, or use stability selection (Parts 3 & 5).
create_figure('weights raw'); axis off; montage(pm_r, bmrk3);
X = double(bmrk3.dat'); Y = bmrk3.Y;
pm_r = bootstrap(pm_r, X, Y, 'nboot', 1000, 'groups', bmrk3.metadata_table.subject_id);
[~, si] = weight_map_object(pm_r, bmrk3);
si = threshold(si, .01, 'unc');                          % bootstrap p < .01, uncorrected
create_figure('weights thresholded'); axis off; montage(si);
%%
%[text] ## 4. What's next
%[text] Part 3 opens the full sklearn-style |predictive_model| API (fit/predict/crossval/bootstrap/permutation, nested-CV tuning, calibration, stability selection); Part 4 runs cross-classification; Part 5 compares algorithms and multiclass ECOC. The |xval_*| wrapper functions (e.g. |xval_SVM|, |xval_SVR|, |xval_lasso_brain|) are single-call alternatives.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
