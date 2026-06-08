%% Multivariate decoding — Part 3: the sklearn-style predictive_model API
%
% Part 3 of the multivariate-decoding series (Part 1 classification basics;
% Part 2 classification & regression; Part 4 cross-classification; Part 5
% algorithms & tuning). See the matching .md for the full narrative,
% nested-CV section, and code map.
%
% All of Part 1's workflow now in the new sklearn-style
% @predictive_model API: construct -> fit -> cross-validate ->
% bootstrap -> threshold -> visualise, plus calibration, permutation
% testing, hyperparameter search, and feature selection. Same DPSP
% Hot-vs-Warm dataset.
%
% Open this file in MATLAB and "Save As" .mlx to get a Live Script.

%% 1. Load DPSP
hw_obj = load_image_set('DPSP_hotwarm', 'noverbose');
rf_obj = load_image_set('DPSP_rejectorfriend', 'noverbose');

X  = double(hw_obj.dat');
Y  = hw_obj.Y;
id = grp2idx(hw_obj.metadata_table.subj_id);

%%


[cverr, stats, optout, pm_svm] = predict(hw_obj, 'algorithm_name', 'cv_svm', 'nfolds', 5, 'newapi');

confusion_matrix(pm_svm);

montage(pm_svm);
surface(pm_svm);

%%
[cverr, stats, optout, pm_pcr] = predict(hw_obj, 'algorithm_name', 'cv_pcr', 'nfolds', 5, 'newapi');

%% 2. Construct a predictive_model
pm = predictive_model( ...
    'algorithm',     'svm', ...
    'task',          'classification', ...
    'modeloptions',  {'KernelFunction', 'linear'}, ...
    'random_state',  2026);

%% 3. fit / predict / score
pm = fit(pm, X, Y, 'id', id);
[yhat, scores] = predict(pm, X);
in_sample_acc = score(pm, X, Y);
fprintf('In-sample balanced accuracy: %.3f\n', in_sample_acc);

%% 4. Cross-validation
pm = crossval(pm, X, Y, ...
    'cv',      cv_splitter.stratified_group_kfold(5), ...
    'groups',  id, ...
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

%% 4b. Repeated cross-validation (more stable estimate)
% A single k-fold split wobbles with the luck of the partition; repeated
% k-fold runs the whole CV n times and averages (cost scales linearly).
pm_rep = crossval(clone(pm), X, Y, ...
    'cv',      cv_splitter.repeated_kfold(5, 5), ...   % 5 folds x 5 repeats
    'groups',  id, ...
    'scoring', 'balanced_accuracy');
fprintf('Repeated-CV balanced accuracy: %.1f%% (sd across folds %.3f)\n', ...
    100 * pm_rep.error_metrics.balanced_accuracy.value, ...
    std(pm_rep.error_metrics.balanced_accuracy_per_fold.value));

%% 5. Visualise the raw (pre-bootstrap) weight map
% weight_map_object caches the @statistic_image on pm (so montage(pm) /
% surface(pm) work later with no source); the second output is the image.
[pm, si_raw] = weight_map_object(pm, hw_obj);
create_figure('SVM weights -- raw'); axis off
montage(pm);            % uses the cached map
snapnow
create_figure('SVM weights -- surface'); axis off
surface(pm);
snapnow

%% 6. Bootstrap inference
% Use a smaller nboot here for tutorial speed; bump to 1000+ for real use.
pm = bootstrap(pm, X, Y, 'nboot', 200, 'groups', id);

fprintf('After bootstrap (nboot=200):\n');
fprintf('  z range = [%g .. %g]\n', min(pm.weights.z), max(pm.weights.z));
fprintf('  empirical p range = [%g .. %g]\n', min(pm.weights.p), max(pm.weights.p));
fprintf('  FDR threshold = %.4g\n', pm.weights.fdr_thr);
fprintf('  FDR-significant voxels = %d / %d\n', sum(pm.weights.fdr_sig), numel(pm.weights.w));

%% 7. Visualise the FDR-thresholded weight map
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

%% 8. Permutation test (is the cv score significantly better than chance?)
% Small nperm for tutorial; bump to 1000+ for real use.
% DPSP is paired (Hot + Warm per subject) — 'auto' detects this and
% uses 'within_subjects' (the gold standard for paired designs).
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

%% 9. Calibrated probabilities
pm_cal = calibrate(pm, X, Y);
P = predict_proba(pm_cal, X);
fprintf('Calibrated probability: mean = %.3f, range [%.3f .. %.3f]\n', ...
    mean(P), min(P), max(P));

%% 10. Hyperparameter search (tiny grid for speed)
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

%% 10b. Shrinkage sweeps (SVM C; LASSO-PCR components)
% Each hyperparameter has an interior optimum. The score AT the chosen point
% is optimistically biased (you picked it on these folds) — nest the search
% for an honest estimate (see the .md).
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

%% 11. Univariate feature selection
pm_fs = predictive_model('algorithm','svm','task','classification', ...
                          'modeloptions', {'KernelFunction','linear'}, ...
                          'random_state', 2026);
pm_fs = select_features(pm_fs, X, Y, 'k', 1000, 'verbose', false);
fprintf('select_features: kept %d / %d features (top-1000 by univariate t-test)\n', ...
    pm_fs.diagnostics.feature_selection.n_selected, ...
    pm_fs.diagnostics.feature_selection.n_input);

%% 12. Rejector vs Friend (same pipeline)
X  = double(rf_obj.dat');
Y  = rf_obj.Y;
id = grp2idx(rf_obj.metadata_table.subj_id);

pm_rf = predictive_model('algorithm','svm','task','classification', ...
                          'modeloptions', {'KernelFunction','linear'}, ...
                          'random_state', 2026);
pm_rf = crossval(pm_rf, X, Y, 'cv', cv_splitter.stratified_group_kfold(5), ...
                  'groups', id, 'scoring', 'balanced_accuracy');
fprintf('\nRejector-vs-Friend cv balanced accuracy: %.1f%%\n', ...
    100 * pm_rf.error_metrics.balanced_accuracy.value);
