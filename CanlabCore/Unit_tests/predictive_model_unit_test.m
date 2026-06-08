function predictive_model_unit_test()
% predictive_model_unit_test  End-to-end test of the new sklearn-style
% @predictive_model API, using DPSP_hotwarm as the test bed.
%
% Small nboot / nperm values are used so the test runs in ~minutes,
% not hours. The point is to exercise every public method, not to
% produce final-quality inference.
%
% Run: cd to this directory, then `predictive_model_unit_test`.

    fprintf('=== predictive_model_unit_test ===\n');

    % --- Load DPSP Hot vs Warm ---
    hw_obj = load_image_set('DPSP_hotwarm', 'noverbose');
    X  = double(hw_obj.dat');
    Y  = hw_obj.Y;
    id = grp2idx(hw_obj.metadata_table.subj_id);
    fprintf('  Loaded DPSP_hotwarm: %d obs, %d features, %d subjects\n', ...
        numel(Y), size(X, 2), numel(unique(id)));

    % --- 1. Construction + hyperparameters ---
    pm = predictive_model('algorithm', 'svm', ...
                          'task',      'classification', ...
                          'random_state', 0);
    assert(strcmp(pm.algorithm, 'svm'),             'algorithm');
    assert(strcmp(pm.task,      'classification'),  'task');
    assert(~pm.is_fitted,                           'starts unfitted');
    assert(pm.is_classifier,                        'is_classifier');
    fprintf('  1. Construction + clone OK\n');

    % --- 2. fit (in-sample) ---
    pm = fit(pm, X, Y, 'id', id);
    assert(pm.is_fitted,                            'is_fitted after fit');
    assert(strcmp(pm.fit_type, 'insample'),         'fit_type=insample');
    assert(numel(pm.weights.w) <= size(X, 2),       'weights.w shape');
    fprintf('  2. fit OK (fit_type=%s, |w|=%d)\n', pm.fit_type, numel(pm.weights.w));

    % --- 3. predict ---
    [yhat, scores] = predict(pm, X);
    assert(numel(yhat) == numel(Y),                 'predict yhat length');
    assert(size(scores, 1) == numel(Y),             'predict scores length');
    fprintf('  3. predict OK (train acc=%.2f)\n', mean(yhat == Y));

    % --- 4. score ---
    pm.scorer = cv_scorer.balanced_accuracy();
    v = score(pm, X, Y);
    assert(v >= 0 && v <= 1,                        'score in [0,1]');
    fprintf('  4. score (balanced_accuracy) = %.3f\n', v);

    % --- 5. crossval with stratified_group_kfold ---
    pm = predictive_model('algorithm','svm','task','classification', ...
                          'modeloptions', {'KernelFunction','linear'}, ...
                          'random_state', 0);
    pm = crossval(pm, X, Y, ...
                  'cv',      cv_splitter.stratified_group_kfold(5), ...
                  'groups',  id, ...
                  'scoring', 'balanced_accuracy');
    assert(strcmp(pm.fit_type, 'crossval'),                              'fit_type=crossval');
    assert(pm.cv_partition.nfolds == 5,                                  '5 folds');
    assert(numel(pm.fitted_values.yfit) == numel(Y),                     'cv yfit length');
    assert(~isnan(pm.error_metrics.balanced_accuracy.value),             'metric populated');
    fprintf('  5. crossval (stratified_group_kfold(5)) bal_acc = %.3f\n', ...
        pm.error_metrics.balanced_accuracy.value);

    % --- 6. bootstrap (small nboot for speed) ---
    pm = bootstrap(pm, X, Y, 'nboot', 25, 'groups', id, 'verbose', false);
    assert(size(pm.weights.boot_w, 2) == 25,                             'nboot=25 samples');
    assert(numel(pm.weights.z)  == numel(pm.weights.w),                  'z length matches w');
    assert(numel(pm.weights.p)  == numel(pm.weights.w),                  'p length matches w');
    assert(all(pm.weights.p >= 0 & pm.weights.p <= 1),                   'p in [0,1]');
    assert(all(abs(pm.weights.z) < 1e6),                                 'z floored not 1e15');
    p_floor = 2 / (25 + 1);
    assert(min(pm.weights.p) >= p_floor - eps,                           'empirical p floor');
    fprintf('  6. bootstrap (nboot=25) OK\n');
    fprintf('       z range = [%g .. %g]; p floor = %.4f\n', ...
        min(pm.weights.z), max(pm.weights.z), min(pm.weights.p));

    % --- 7. weight_map_object: pre- and post-bootstrap thresholded, + caching ---
    [pm_w, si_full] = weight_map_object(pm, hw_obj);
    [~,    si_thr]  = weight_map_object(pm, hw_obj, 'use', 'thresh_fdr');
    assert(isa(si_full, 'statistic_image'),                              'si_full class');
    assert(size(si_full.dat, 1) == size(hw_obj.dat, 1),                  'si_full voxel count');
    assert(size(si_thr.dat,  1) == size(hw_obj.dat, 1),                  'si_thr voxel count');
    assert(~isempty(pm_w.weights.weight_obj),                            'weight_obj cached on pm');
    assert(isa(pm_w.weights.weight_obj, 'statistic_image'),             'cached weight_obj class');
    fprintf('  7. weight_map_object (full + thresh_fdr + cache) OK\n');
    fprintf('       full-w nonzero: %d, FDR-thresh nonzero: %d\n', ...
        sum(si_full.dat ~= 0), sum(si_thr.dat ~= 0));

    % --- 8. permutation_test (small nperm) ---
    pm_for_perm = predictive_model('algorithm','svm','task','classification', ...
                                   'modeloptions', {'KernelFunction','linear'}, ...
                                   'random_state', 0);
    pm_for_perm.cv = cv_splitter.stratified_group_kfold(3);
    pm_for_perm = permutation_test(pm_for_perm, X, Y, ...
                                   'nperm', 10, 'groups', id, ...
                                   'verbose', false);
    assert(numel(pm_for_perm.permutation_results.null_scores) == 10,     'nperm=10 null scores');
    assert(pm_for_perm.permutation_results.p_value > 0,                  'p_value > 0');
    assert(pm_for_perm.permutation_results.p_value <= 1,                 'p_value <= 1');
    % DPSP is a paired design (Hot+Warm per subject) — auto should pick within_subjects.
    assert(strcmp(pm_for_perm.permutation_results.permutation, 'within_subjects'), ...
        'auto should pick within_subjects on paired DPSP');
    fprintf('  8. permutation_test (nperm=10, %s) observed=%.3f, p=%.3f\n', ...
        pm_for_perm.permutation_results.permutation, ...
        pm_for_perm.permutation_results.observed, ...
        pm_for_perm.permutation_results.p_value);

    % --- 9. calibrate + predict_proba ---
    pm_cal = calibrate(pm, X, Y);
    assert(isfield(pm_cal.fitted_values, 'calibrator'),                  'calibrator present');
    P = predict_proba(pm_cal, X);
    assert(all(P >= 0 & P <= 1),                                         'predict_proba in [0,1]');
    fprintf('  9. calibrate + predict_proba: P mean=%.3f\n', mean(P));

    % --- 10. select_features ---
    pm_fs = predictive_model('algorithm','svm','task','classification');
    pm_fs = select_features(pm_fs, X, Y, 'k', 500, 'verbose', false);
    assert(pm_fs.diagnostics.feature_selection.n_selected == 500,        'k=500 selected');
    fprintf(' 10. select_features OK (kept 500)\n');

    % --- 11. visualisation methods (smoke test; no display assertions) ---
    ROC = rocplot(pm, 'noplot');
    assert(isfield(ROC, 'AUC') && ROC.AUC >= 0 && ROC.AUC <= 1,          'rocplot AUC in [0,1]');
    h = plot(pm); close(h);                    % classification dispatcher (violin + ROC)
    confusionchart(pm); close(gcf);
    fprintf(' 11. rocplot/plot/confusionchart OK (AUC=%.3f)\n', ROC.AUC);

    % --- 12. montage / surface delegates (smoke test) ---
    montage(pm, hw_obj); close all force;
    fprintf(' 12. montage(pm, source) OK\n');

    % --- 13. grid_search (tiny grid) ---
    pm_gs = predictive_model('algorithm','svm','task','classification', ...
                             'cv', cv_splitter.stratified_group_kfold(3));
    pm_gs = grid_search(pm_gs, X, Y, struct('BoxConstraint', [0.1 1]), ...
                        'groups', id, 'verbose', false);
    assert(isfield(pm_gs.diagnostics.grid_search, 'best_score'),         'grid_search best_score');
    fprintf(' 13. grid_search OK (best %s=%.3f)\n', ...
        pm_gs.scorer.name, pm_gs.diagnostics.grid_search.best_score);

    % --- 14. stability_selection (small nboot) ---
    pm_ss = stability_selection( ...
        predictive_model('algorithm','linear_svm','task','classification'), ...
        X, Y, 'nboot', 20, 'k', 1000, 'threshold', 0.6, 'groups', id, 'verbose', false);
    ss = pm_ss.diagnostics.stability_selection;
    assert(numel(ss.selection_freq) == size(X, 2),                       'selection_freq length');
    assert(all(ss.selection_freq >= 0 & ss.selection_freq <= 1),         'selection_freq in [0,1]');
    fprintf(' 14. stability_selection OK (%d stable voxels)\n', ss.n_stable);

    % --- 15. @pipeline: PCA -> svm, leakage-free crossval + weight back-projection ---
    est  = predictive_model('algorithm','svm','task','classification');
    pipe = pipeline({ {'pca','k',20} }, est);
    pipe = crossval(pipe, X, Y, 'groups', id, 'cv', cv_splitter.stratified_group_kfold(5));
    assert(strcmp(pipe.fit_type, 'crossval'),                            'pipeline crossval fit_type');
    assert(numel(pipe.weights.w) == size(X, 2),                          'pipeline back-projected weight length');
    si_pipe = weight_map_object(pipe, hw_obj);
    assert(isa(si_pipe, 'statistic_image'),                              'pipeline weight_map_object class');
    fprintf(' 15. pipeline (PCA->svm) OK (cv bal_acc=%.3f)\n', ...
        pipe.error_metrics.balanced_accuracy.value);

    % --- 16. fmri_data.predict 'newapi' routing ---
    [cverr, st, ~, pm_api] = predict(hw_obj, 'algorithm_name', 'cv_svm', ...
                                     'nfolds', 5, 'newapi', 'verbose', 0);
    assert(isa(pm_api, 'predictive_model'),                              'newapi returns predictive_model');
    assert(strcmp(pm_api.fit_type, 'crossval'),                          'newapi fit_type=crossval');
    assert(isfield(st, 'dist_from_hyperplane_xval'),                     'newapi stats has xval distances');
    fprintf(' 16. fmri_data.predict ''newapi'' OK (cverr=%.3f)\n', cverr);

    fprintf('\npredictive_model_unit_test: PASS\n');
end
