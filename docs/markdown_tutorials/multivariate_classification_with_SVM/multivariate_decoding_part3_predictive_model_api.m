%[text] # Multivariate decoding — Part 3: the sklearn-style predictive\_model API
%[text] Part 3 of the multivariate-decoding series (Part 1 classification basics; Part 2 classification & regression; Part 4 cross-classification; Part 5 algorithms & tuning). See the matching .md for the full narrative, nested-CV section, and code map.
%[text] All of Part 1's workflow now in the new sklearn-style @predictive\_model API: construct -\> fit -\> cross-validate -\> bootstrap -\> threshold -\> visualise, plus calibration, permutation testing, hyperparameter search, and feature selection. Same DPSP Hot-vs-Warm dataset.
%[text] Open this file in MATLAB and "Save As" .mlx to get a Live Script.
%[text] ## 1. Load DPSP
hw_obj = load_image_set('DPSP_hotwarm', 'noverbose');
rf_obj = load_image_set('DPSP_rejectorfriend', 'noverbose');

X  = double(hw_obj.dat');
Y  = hw_obj.Y;
id = grp2idx(hw_obj.metadata_table.subj_id);
%%
%[text] ## Use fmri\_data.predict wrapper to create two demo predictive models
%[text] - One classification model (pm\_svm)
%[text] - One regression model (pm\_pcr) \
[cverr, stats, optout, pm_svm] = predict(hw_obj, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'newapi');

create_figure('confusionmat')
confusion_matrix(pm_svm);

create_figure('montage')
montage(pm_svm);

create_figure('surface')
surface(pm_svm);
%%
%[text] The regression model is technically not the correct thing to do with binary outcomes, but we do it here to illustrate the methods:
[cverr, stats, optout, pm_pcr] = predict(hw_obj, 'algorithm_name', 'cv_pcr', 'nfolds', 5, 'newapi');
%%
%[text] ## 2. Construct a predictive\_model object using the API
%[text] Construction sets **hyperparameters only** — no data is touched yet; you're declaring *what kind of model* you want. The two options that matter most are **algorithm** (which learner; see the registry in the .md) and **task** (|'classification'| vs |'regression'|, which sets the default scorer and how Y is read). If you omit |task| it is **inferred from Y** (≤2 unique values → classification); if you omit |scorer| it defaults to balanced\_accuracy (classification) / r2 (regression); if you omit |cv| it defaults to |stratified\_kfold(5)| — which is NOT group-aware, so for within-subject data pass a grouped splitter (below). Setting |random\_state| makes folds/bootstraps reproducible.
pm = predictive_model( ...
    'algorithm',     'svm', ...                        % which learner
    'task',          'classification', ...             % classification | regression
    'modeloptions',  {'KernelFunction', 'linear'}, ... % passed through to the fitter
    'random_state',  2026);
%%
%[text] ## 3. fit / predict / score (the sklearn triad)
%[text] These three verbs are the foundation; |crossval|/|bootstrap|/… build on them.
%[text] - **fit(pm, X, Y)** — *train* on all the data. Populates |ml\_model|, |weights.w|, and the intercept. Use for a final "ship-it" model or before applying to a separate test set.
%[text] - **predict(pm, Xnew)** — *apply* to new rows; returns labels/values and continuous scores. It applies the **full decision function w·x + b** (weights AND intercept), plus the omitted-features mask and any standardization recorded at fit time. The intercept is applied in-sample, on test data, and in every CV fold.
%[text] - **score(pm, X, Y)** — *evaluate*: predict on X and compare to Y with the model's scorer (balanced accuracy by default).
pm = fit(pm, X, Y, 'id', id);
[yhat, scores] = predict(pm, X);
in_sample_acc = score(pm, X, Y);
fprintf('In-sample balanced accuracy: %.3f\n', in_sample_acc);
%[text] **Bias-term gotcha:** the score from predict equals |X\*pm.weights.w + pm.ml\_model.Bias| (SVM/fitclinear) or |+ pm.ml\_model.intercept| (PCR/lassoPCR). But |pm.weights.w| (and the weight map from |weight\_map\_object|) is the **slope only, no intercept** — correct for a brain map, but if you apply a pattern by hand (|Xnew\*w|, dot-product / apply\_mask / image\_similarity) you must add the intercept yourself. Easiest: just call |predict(pm, Xnew)|.
%%
%[text] ## 4. Cross-validation — the splitter API
%[text] An in-sample fit only shows how well the model *memorised* the data. Cross-validation estimates **out-of-sample** performance: split into k folds; train on k−1, predict the held-out fold; every observation is predicted once by a model that never saw it. A **cv\_splitter** only decides which rows go in which fold (it knows nothing about the model): |stratified| keeps class proportions per fold; |group| keeps all rows of a subject in the same fold (essential — otherwise a subject in both train and test leaks); |stratified\_group\_kfold| does both and is the right default for paired brain designs like DPSP.
pm = crossval(pm, X, Y, ...
    'cv',      cv_splitter.stratified_group_kfold(5), ...   % the splitter
    'groups',  id, ...                                      % subject id per row
    'scoring', 'balanced_accuracy');

fprintf('CV balanced accuracy: %.1f%% (per-fold: %s)\n', ...
    100 * pm.error_metrics.balanced_accuracy.value, ...
    mat2str(100 * pm.error_metrics.balanced_accuracy_per_fold.value', 3));

% summary(pm) prints provenance + the task-relevant performance block
% (balanced accuracy, AUC, sensitivity, specificity, PPV, NPV, d); the
% Inference line fills in after bootstrap/permutation_test/calibrate.
summary(pm);
% report_accuracy(pm) prints just the metric block and returns it as a struct:
acc = report_accuracy(pm, 'noverbose');
%%
%[text] ## 4b. Repeated cross-validation (more stable estimate)
%[text] A single k-fold split wobbles with the luck of the partition; repeated k-fold runs the whole CV n times and averages (cost scales linearly).
pm_rep = crossval(clone(pm), X, Y, ...
    'cv',      cv_splitter.repeated_kfold(5, 5), ...   % 5 folds x 5 repeats
    'groups',  id, ...
    'scoring', 'balanced_accuracy');
fprintf('Repeated-CV balanced accuracy: %.1f%% (sd across folds %.3f)\n', ...
    100 * pm_rep.error_metrics.balanced_accuracy.value, ...
    std(pm_rep.error_metrics.balanced_accuracy_per_fold.value));
%%
%[text] ## 5. Visualise the raw (pre-bootstrap) weight map
%[text] weight\_map\_object caches the @statistic\_image on pm (so montage(pm) / surface(pm) work later with no source); the second output is the image.
[pm, si_raw] = weight_map_object(pm, hw_obj);
create_figure('SVM weights -- raw'); axis off
montage(pm);            % uses the cached map
snapnow
create_figure('SVM weights -- surface'); axis off
surface(pm);
snapnow
%%
%[text] ## 6. Bootstrap inference
%[text] Use a smaller nboot here for tutorial speed; bump to 1000+ for real use.
pm = bootstrap(pm, X, Y, 'nboot', 200, 'groups', id);

fprintf('After bootstrap (nboot=200):\n');
fprintf('  z range = [%g .. %g]\n', min(pm.weights.z), max(pm.weights.z));
fprintf('  empirical p range = [%g .. %g]\n', min(pm.weights.p), max(pm.weights.p));
fprintf('  FDR threshold = %.4g\n', pm.weights.fdr_thr);
fprintf('  FDR-significant voxels = %d / %d\n', sum(pm.weights.fdr_sig), numel(pm.weights.w));
%%
%[text] ## 7. Visualise the FDR-thresholded weight map
[~, si_thr] = weight_map_object(pm, hw_obj, 'use', 'thresh_fdr');
create_figure('SVM weights -- FDR thresholded'); axis off
montage(si_thr);
snapnow

% Or via statistic_image's own threshold(). For a strongly regularised model
% the FDR threshold can be empty (see the bootstrap caveat in the .md), so we
% threshold at uncorrected p < .01 here and also show the surface.
[~, si] = weight_map_object(pm, hw_obj);
si = threshold(si, .01, 'unc');

create_figure('SVM weights -- thresholded montage'); axis off
montage(si);
snapnow

create_figure('SVM weights -- thresholded surface'); axis off
surface(si, 'foursurfaces_hcp');
snapnow
%%
%[text] ## 8. Permutation test (is the cv score significantly better than chance?)
%[text] Small nperm for tutorial; bump to 1000+ for real use. DPSP is paired (Hot + Warm per subject) — 'auto' detects this and uses 'within\_subjects' (the gold standard for paired designs).
pm = permutation_test(pm, X, Y, 'nperm', 50, 'groups', id);
fprintf('Permutation test (%s): observed = %.3f, null mean = %.3f, p = %.4f\n', ...
    pm.permutation_results.permutation, ...
    pm.permutation_results.observed, ...
    mean(pm.permutation_results.null_scores, 'omitnan'), ...
    pm.permutation_results.p_value);
fprintf('  scheme: %s\n', pm.permutation_results.permutation_descrip);

% Visualise the null: histogram with a dashed alpha=.05 critical line and a
% solid observed-score line. (plot(pm) also draws this automatically whenever
% permutation results are present.) For publication use nperm >= 5000.
create_figure('permutation null');
plot_permutation(pm);
snapnow

% Force a specific scheme — useful for reporting or to compare
% paired vs free-permutation p-values:
%   pm = permutation_test(pm, X, Y, 'nperm', 50, 'groups', id, ...
%                          'permutation', 'within_subjects');
%%
%[text] ## 9. Calibrated probabilities
pm_cal = calibrate(pm, X, Y);
P = predict_proba(pm_cal, X);
fprintf('Calibrated probability: mean = %.3f, range [%.3f .. %.3f]\n', ...
    mean(P), min(P), max(P));
%%
%[text] ## 10. Hyperparameter search (tiny grid for speed)
pg = struct('BoxConstraint', [0.1 1 10], 'KernelScale', [1 2]);

% Re-fit a clean pm before searching (so grid_search doesn't inherit
% the modeloptions that bootstrap/calibrate appended).
pm_gs = predictive_model('algorithm','svm','task','classification', ...
                          'modeloptions', {'KernelFunction','linear'}, ...
                          'random_state', 2026);
pm_gs.cv = cv_splitter.stratified_group_kfold(3);
pm_gs = grid_search(pm_gs, X, Y, pg, 'groups', id, 'verbose', false);
fprintf('Best params: %s, best score = %.3f\n', ...
    strjoin(string(pm_gs.diagnostics.grid_search.best_params(2:2:end)), ' / '), ...
    pm_gs.diagnostics.grid_search.best_score);
%%
%[text] ## 10b. Shrinkage sweeps (SVM C; LASSO-PCR components)
%[text] Each hyperparameter has an interior optimum. The score AT the chosen point is optimistically biased (you picked it on these folds) — nest the search for an honest estimate (see the .md).
Cgrid = logspace(-3, 2, 6);
accC  = nan(size(Cgrid));
for i = 1:numel(Cgrid)
    pmc = crossval(predictive_model('algorithm','svm','task','classification', ...
              'modeloptions', {'BoxConstraint', Cgrid(i)}), X, Y, 'groups', id, ...
              'cv', cv_splitter.stratified_group_kfold(5), 'scoring', 'balanced_accuracy');
    accC(i) = pmc.error_metrics.balanced_accuracy.value;
end

bmrk3 = load_image_set('bmrk3', 'noverbose');
Xr = double(bmrk3.dat');  Yr = bmrk3.Y;  idr = bmrk3.metadata_table.subject_id;
lnum = [5 10 20 40 80];                  % # components kept before the OLS refit
r2L  = nan(size(lnum));
for i = 1:numel(lnum)
    pml = crossval(predictive_model('algorithm','lassopcr','task','regression', ...
              'modeloptions', {'lasso_num', lnum(i)}), Xr, Yr, 'groups', idr, ...
              'cv', cv_splitter.group_kfold(5), 'scoring', 'r2');
    r2L(i) = pml.error_metrics.predicted_r2.value;
end

create_figure('shrinkage sweeps', 1, 2);
subplot(1,2,1); semilogx(Cgrid, accC, 'o-', 'LineWidth', 2);
xlabel('SVM C (BoxConstraint)'); ylabel('CV balanced accuracy'); title('DPSP Hot/Warm');
subplot(1,2,2); plot(lnum, r2L, 'o-', 'LineWidth', 2);
xlabel('LASSO-PCR components kept'); ylabel('predicted R^2'); title('bmrk3 pain');
snapnow
%%
%[text] ## 11. Univariate feature selection
pm_fs = predictive_model('algorithm','svm','task','classification', ...
                          'modeloptions', {'KernelFunction','linear'}, ...
                          'random_state', 2026);
pm_fs = select_features(pm_fs, X, Y, 'k', 1000, 'verbose', false);
fprintf('select_features: kept %d / %d features (top-1000 by univariate t-test)\n', ...
    pm_fs.diagnostics.feature_selection.n_selected, ...
    pm_fs.diagnostics.feature_selection.n_input);
%%
%[text] ## 12. Rejector vs Friend (same pipeline)
X  = double(rf_obj.dat');
Y  = rf_obj.Y;
id = grp2idx(rf_obj.metadata_table.subj_id);

pm_rf = predictive_model('algorithm','svm','task','classification', ...
                          'modeloptions', {'KernelFunction','linear'}, ...
                          'random_state', 2026);
pm_rf = crossval(pm_rf, X, Y, 'cv', cv_splitter.stratified_group_kfold(5), ...
                  'groups', id, 'scoring', 'balanced_accuracy');
fprintf('\nRejector-vs-Friend cv balanced accuracy: %.1f%%n', ...
    100 * pm_rf.error_metrics.balanced_accuracy.value);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
