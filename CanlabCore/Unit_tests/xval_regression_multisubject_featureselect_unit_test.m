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

    % Categorised <-> legacy alias agreement on key fields.
    % (mean_vox_weights is only populated when pcsquash is enabled;
    % this no-PCA test doesn't exercise that path.)
    assert(isequaln(pm.subjfit,        pm.fitted_values.subjfit),               'subjfit');
    % vox_weights renamed to weights.w_perfold in consolidated layout
    assert(isequaln(pm.vox_weights,    pm.weights.w_perfold),                   'vox_weights / w_perfold');
    % error_metrics entries are (value, descrip) tuples.
    assert(isequaln(pm.pred_err,       pm.error_metrics.pred_err.value),        'pred_err');
    assert(isequaln(pm.pred_err_null,  pm.error_metrics.pred_err_null.value),   'pred_err_null');
    assert(isequaln(pm.var_reduction,  pm.error_metrics.var_reduction.value),   'var_reduction');
    assert(~isempty(pm.Y_orig), 'Y_orig empty (should alias Y)');
    fprintf('  Categorised <-> legacy alias round-trip OK\n');

    % Verify vox_weights matches the full feature space (the bug we fixed
    % was that this was the *selected*-features count instead).
    assert(size(pm.vox_weights, 1) == size(X, 2), ...
        'vox_weights should span full feature count (%d), got %d', ...
        size(X, 2), size(pm.vox_weights, 1));
    fprintf('  vox_weights shape OK: %d x %d (features x folds)\n', ...
        size(pm.vox_weights, 1), size(pm.vox_weights, 2));

    assert(~isempty(pm.r_squared) && pm.r_squared(1) > 0, ...
        'expected r_squared > 0; got %.3f', pm.r_squared(1));
    fprintf('  r_squared = %.3f, pred_err = %.3f, pred_err_null = %.3f\n', ...
        pm.r_squared(1), pm.pred_err(1), pm.pred_err_null(1));

    pm.validate_object('noverbose');
    fprintf('xval_regression_multisubject_featureselect_unit_test: PASS\n');
end
