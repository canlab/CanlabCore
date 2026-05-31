function xval_regression_multisubject_unit_test()
% xval_regression_multisubject_unit_test  Smoke test using DPSP sample data.
%
% Verifies that xval_regression_multisubject returns a populated
% @predictive_model object whose categorised sub-structs and legacy
% flat aliases agree, on a real (small) brain dataset.
%
% Data: Sample_datasets/DPSP_pain_rejection_participant_maps
%       Single-subject Hot - Warm contrast maps, regressed against a
%       synthetic continuous outcome to exercise the regression code path.
%
% Run: cd to this directory, then `xval_regression_multisubject_unit_test`.

    fprintf('=== xval_regression_multisubject_unit_test ===\n');

    % --- Load DPSP Hot - Warm contrast across subjects ---
    canlabcore_dir = fileparts(fileparts(which('fmri_data')));
    sample_dir = fullfile(canlabcore_dir, 'Sample_datasets', ...
                          'DPSP_pain_rejection_participant_maps');
    H = load(fullfile(sample_dir, 'DPSP_single_subject_images_hot.mat'));
    W = load(fullfile(sample_dir, 'DPSP_single_subject_images_warm.mat'));
    hot_obj  = H.single_subject_images_hot;
    warm_obj = W.single_subject_images_warm;
    hot_vs_warm = image_math(hot_obj, warm_obj, 'minus');

    % Keep the test light: ~5000 random voxels, all subjects.
    rng(0);
    [p, n] = size(hot_vs_warm.dat);
    keep_vox = randsample(p, min(5000, p));
    X = double(hot_vs_warm.dat(keep_vox, :))';       % n x v

    % Synthetic continuous outcome predictable from X with sparse weights.
    b_true = zeros(size(X, 2), 1);
    b_true(1:50) = randn(50, 1);
    Y = X * b_true + 0.5 * std(X * b_true) * randn(n, 1);

    % --- Run xval_regression_multisubject with OLS + PCA on a single dataset ---
    pmodel_obj = xval_regression_multisubject('ols', {Y}, {X}, ...
        'pca', 'ndims', 10, 'holdout_method', 'balanced4', 'noverbose');

    % --- Class invariants ---
    assert(isa(pmodel_obj, 'predictive_model'), ...
        'Expected output of class predictive_model, got %s', class(pmodel_obj));
    assert(pmodel_obj.is_fitted, 'Returned object reports is_fitted = false');
    fprintf('  class = %s, is_fitted = %d\n', class(pmodel_obj), pmodel_obj.is_fitted);

    % --- Categorised <-> legacy alias agreement ---
    assert(isequal(pmodel_obj.subjfit,            pmodel_obj.fitted_values.subjfit),               'subjfit');
    assert(isequal(pmodel_obj.subjbetas,          pmodel_obj.weights.subjbetas),                   'subjbetas');
    assert(isequal(pmodel_obj.mean_vox_weights,   pmodel_obj.weights.mean_vox_weights),            'mean_vox_weights');
    % error_metrics entries are (value, descrip) tuples; legacy alias unwraps .value.
    assert(isequal(pmodel_obj.pred_err,           pmodel_obj.error_metrics.pred_err.value),        'pred_err');
    assert(isequal(pmodel_obj.pred_err_null,      pmodel_obj.error_metrics.pred_err_null.value),   'pred_err_null');
    assert(isequal(pmodel_obj.var_reduction,      pmodel_obj.error_metrics.var_reduction.value),   'var_reduction');
    assert(isequal(pmodel_obj.r_each_subject,     pmodel_obj.error_metrics.r_each_subject.value),  'r_each_subject');
    % After consolidation: Y/Y_orig are top-level; INPUTS is a Dependent
    % alias for the top-level inputParameters struct.
    assert(~isempty(pmodel_obj.Y_orig),          'Y_orig empty (should alias Y)');
    assert(~isempty(pmodel_obj.inputParameters), 'inputParameters empty');
    fprintf('  Categorised <-> legacy alias round-trip OK\n');

    % --- Shape / sanity ---
    assert(numel(pmodel_obj.subjfit) == 1,                       'subjfit cell length should equal num datasets');
    assert(numel(pmodel_obj.subjfit{1}) == n,                    'subjfit{1} length should equal n');
    assert(size(pmodel_obj.mean_vox_weights, 1) == size(X, 2),   'mean_vox_weights should have one row per voxel');
    assert(~isempty(pmodel_obj.r_squared) && pmodel_obj.r_squared(1) >  0, ...
        'OLS+PCA should explain non-trivial variance on this synthetic outcome');
    fprintf('  Shape/sanity OK: r_squared = %.3f, pred_err = %.3f, pred_err_null = %.3f\n', ...
        pmodel_obj.r_squared(1), pmodel_obj.pred_err(1), pmodel_obj.pred_err_null(1));

    % --- validate_object accepts the returned object ---
    pmodel_obj.validate_object('noverbose');
    fprintf('  validate_object OK\n');

    fprintf('xval_regression_multisubject_unit_test: PASS\n');
end
