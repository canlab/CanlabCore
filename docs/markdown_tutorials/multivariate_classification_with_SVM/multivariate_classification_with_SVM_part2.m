%% Multivariate classification with SVM — Part 2
%
% All of Part 1's workflow now in the new sklearn-style
% @predictive_model API: construct -> fit -> cross-validate ->
% bootstrap -> threshold -> visualise, plus calibration, permutation
% testing, hyperparameter search, and feature selection. Same DPSP
% Hot-vs-Warm dataset.
%
% Open this file in MATLAB and "Save As" .mlx to get a Live Script.

%% 1. Load DPSP
hw_obj = load_image_set('DPSP_hotwarm');
rf_obj = load_image_set('DPSP_rejectorfriend');

X  = double(hw_obj.dat');
Y  = hw_obj.Y;
id = grp2idx(hw_obj.metadata_table.subj_id);

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

%% 5. Visualise the raw (pre-bootstrap) weight map
si_raw = weight_image(pm, hw_obj);
create_figure('SVM weights -- raw'); axis off
montage(si_raw);
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
si_thr = weight_image(pm, hw_obj, 'use', 'thresh_fdr');
create_figure('SVM weights -- FDR thresholded'); axis off
montage(si_thr);
snapnow

% Or via statistic_image's own threshold():
si = weight_image(pm, hw_obj);
si = threshold(si, .05, 'fdr');
create_figure('SVM weights -- statistic_image threshold(.05, fdr)'); axis off
montage(region(si));
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
