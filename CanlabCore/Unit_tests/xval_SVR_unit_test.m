function xval_SVR_unit_test()
% xval_SVR_unit_test  Smoke test of xval_SVR on the DPSP Hot-Warm contrast.
%
% Predicts a synthetic continuous outcome from Hot-Warm contrast maps
% (one per subject) via cross-validated linear SVR. Asserts that the
% wrapper returns a @predictive_model object with consistent
% categorised <-> legacy alias values.
%
% Run: cd to this directory, then `xval_SVR_unit_test`.

    fprintf('=== xval_SVR_unit_test ===\n');

    canlabcore_dir = fileparts(fileparts(which('fmri_data')));
    sample_dir = fullfile(canlabcore_dir, 'Sample_datasets', ...
                          'DPSP_pain_rejection_participant_maps');
    H = load(fullfile(sample_dir, 'DPSP_single_subject_images_hot.mat'));
    W = load(fullfile(sample_dir, 'DPSP_single_subject_images_warm.mat'));
    hot_vs_warm = image_math(H.single_subject_images_hot, ...
                             W.single_subject_images_warm, 'minus');

    rng(0);
    [p, n] = size(hot_vs_warm.dat);
    keep_vox = randsample(p, min(3000, p));
    X = double(hot_vs_warm.dat(keep_vox, :))';
    % Synthetic predictable continuous outcome.
    b_true = zeros(size(X, 2), 1);
    b_true(1:30) = randn(30, 1);
    Y = X * b_true + 0.5 * std(X * b_true) * randn(n, 1);
    id = (1:n)';

    pmodel_obj = xval_SVR(X, Y, id, 'nooptimize', 'norepeats', ...
                                    'nobootstrap', 'noverbose', 'noplot');

    assert(isa(pmodel_obj, 'predictive_model'), ...
        'Expected predictive_model, got %s', class(pmodel_obj));
    assert(pmodel_obj.is_fitted, 'is_fitted false');
    assert(pmodel_obj.is_regressor, 'Y is continuous, expected is_regressor true');
    fprintf('  class = %s, is_fitted = %d, is_regressor = %d\n', ...
        class(pmodel_obj), pmodel_obj.is_fitted, pmodel_obj.is_regressor);

    % Round-trip categorised vs. legacy aliases (regression-relevant fields only).
    assert(isequaln(pmodel_obj.Y,                    pmodel_obj.inputs.Y),                            'Y');
    assert(isequaln(pmodel_obj.yfit,                 pmodel_obj.fitted_values.yfit),                  'yfit');
    assert(isequaln(pmodel_obj.w,                    pmodel_obj.weights.w),                           'w');
    assert(isequaln(pmodel_obj.prediction_outcome_r, pmodel_obj.error_metrics.prediction_outcome_r),  'prediction_outcome_r');
    assert(isequaln(pmodel_obj.pred_outcome_r,       pmodel_obj.error_metrics.prediction_outcome_r),  'pred_outcome_r alias');
    assert(isequaln(pmodel_obj.d_singleinterval,     pmodel_obj.error_metrics.d_singleinterval),      'd_singleinterval');
    assert(isequaln(pmodel_obj.SVRModel,             pmodel_obj.ml_model),                            'SVRModel');
    fprintf('  Categorised <-> legacy alias round-trip OK\n');

    fprintf('  cv: r = %.3f, d = %.2f\n', pmodel_obj.prediction_outcome_r, pmodel_obj.d_singleinterval);

    pmodel_obj.validate_object('noverbose');
    fprintf('xval_SVR_unit_test: PASS\n');
end
