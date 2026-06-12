function xval_regression_multisubject_featureselect_unit_test()
% xval_regression_multisubject_featureselect_unit_test  DPSP smoke test.
%
% Verifies that xval_regression_multisubject_featureselect returns a
% populated @predictive_model object after the per-fold vox_weights
% padding fix. Predicts a synthetic continuous outcome from
% Hot - Warm DPSP contrast maps.

    fprintf('=== xval_regression_multisubject_featureselect_unit_test ===\n');

    canlabcore_dir = fileparts(fileparts(which('fmri_data')));
    sd = fullfile(canlabcore_dir, 'Sample_datasets', 'DPSP_pain_rejection_participant_maps');
    H = load(fullfile(sd, 'DPSP_single_subject_images_hot.mat'));
    W = load(fullfile(sd, 'DPSP_single_subject_images_warm.mat'));
    hvw = image_math(H.single_subject_images_hot, W.single_subject_images_warm, 'minus');

    rng(0);
    [p, n] = size(hvw.dat);
    keep_vox = randsample(p, min(2000, p));
    X = double(hvw.dat(keep_vox, :))';
    b_true = zeros(size(X, 2), 1);
    b_true(1:30) = randn(30, 1);
    Y = X * b_true + 0.5 * std(X * b_true) * randn(n, 1);

    pm = xval_regression_multisubject_featureselect('ols', {Y}, {X}, ...
        'holdout_method', 'balanced4', 'noverbose');

    assert(isa(pm, 'predictive_model'), 'class mismatch: %s', class(pm));
    assert(pm.is_fitted, 'is_fitted=false');
    fprintf('  class = %s, is_fitted = %d\n', class(pm), pm.is_fitted);

    % Canonical-path access (legacy flat aliases removed).
    % mean_vox_weights is only populated when pcsquash is enabled; this no-PCA test doesn't exercise that path.
    assert(~isempty(pm.fitted_values.subjfit),               'subjfit empty');
    assert(~isempty(pm.weights.w_perfold),                   'weights.w_perfold empty (was vox_weights)');
    assert(~isempty(pm.error_metrics.pred_err.value),        'pred_err empty');
    assert(~isempty(pm.error_metrics.pred_err_null.value),   'pred_err_null empty');
    assert(~isempty(pm.error_metrics.var_reduction.value),   'var_reduction empty');
    assert(~isempty(pm.Y),                                   'Y empty');
    fprintf('  Canonical-path access OK\n');

    % Verify weights.w_perfold matches the full feature space (the bug we fixed
    % was that this was the *selected*-features count instead).
    assert(size(pm.weights.w_perfold, 1) == size(X, 2), ...
        'weights.w_perfold should span full feature count (%d), got %d', ...
        size(X, 2), size(pm.weights.w_perfold, 1));
    fprintf('  w_perfold shape OK: %d x %d (features x folds)\n', ...
        size(pm.weights.w_perfold, 1), size(pm.weights.w_perfold, 2));

    r2 = pm.error_metrics.r_squared.value;
    assert(~isempty(r2) && r2(1) > 0, 'expected r_squared > 0; got %.3f', r2(1));
    fprintf('  r_squared = %.3f, pred_err = %.3f, pred_err_null = %.3f\n', ...
        r2(1), ...
        pm.error_metrics.pred_err.value(1), ...
        pm.error_metrics.pred_err_null.value(1));

    pm.validate_object('noverbose');
    fprintf('xval_regression_multisubject_featureselect_unit_test: PASS\n');
end
