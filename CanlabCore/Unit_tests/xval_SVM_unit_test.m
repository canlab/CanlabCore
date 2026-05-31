function xval_SVM_unit_test()
% xval_SVM_unit_test  Smoke test of xval_SVM on the DPSP Hot vs. Warm task.
%
% Reproduces the binary-classification setup from
% docs/markdown_tutorials/multivariate_classification_with_SVM/
%   multivariate_classification_with_SVM_part1.mlx
% and asserts that the wrapper returns a populated @predictive_model
% object whose categorised sub-structs agree with the legacy flat aliases.
%
% Data: Sample_datasets/DPSP_pain_rejection_participant_maps
%       Single-subject Hot and Warm condition maps.
%
% Run: cd to this directory, then `xval_SVM_unit_test`.

    fprintf('=== xval_SVM_unit_test ===\n');

    % --- Load DPSP Hot and Warm ---
    canlabcore_dir = fileparts(fileparts(which('fmri_data')));
    sample_dir = fullfile(canlabcore_dir, 'Sample_datasets', ...
                          'DPSP_pain_rejection_participant_maps');
    H = load(fullfile(sample_dir, 'DPSP_single_subject_images_hot.mat'));
    W = load(fullfile(sample_dir, 'DPSP_single_subject_images_warm.mat'));
    hot  = H.single_subject_images_hot;
    warm = W.single_subject_images_warm;

    n_hot  = size(hot.dat,  2);
    n_warm = size(warm.dat, 2);

    % Stack into one (subjects+subjects) x voxels matrix
    rng(0);
    p = size(hot.dat, 1);
    keep_vox = randsample(p, min(5000, p));
    X = double([hot.dat(keep_vox, :) warm.dat(keep_vox, :)])';
    Y = [ones(n_hot, 1); -ones(n_warm, 1)];
    id = [(1:n_hot)'; (1:n_warm)'];   % each subject contributes both conditions

    % --- Run xval_SVM (fast: no optimize / no repeats / no bootstrap) ---
    pmodel_obj = xval_SVM(X, Y, id, 'nooptimize', 'norepeats', 'nobootstrap', ...
                                    'noverbose', 'noplot');

    % --- Class invariants ---
    assert(isa(pmodel_obj, 'predictive_model'), ...
        'Expected output of class predictive_model, got %s', class(pmodel_obj));
    assert(pmodel_obj.is_fitted, 'Returned object reports is_fitted = false');
    assert(pmodel_obj.is_classifier, 'Y is binary, expected is_classifier true');
    fprintf('  class = %s, is_fitted = %d, is_classifier = %d\n', ...
        class(pmodel_obj), pmodel_obj.is_fitted, pmodel_obj.is_classifier);

    % --- Categorised <-> legacy alias agreement ---
    % Y, id are now top-level (was pmodel_obj.inputs.Y / .id) — sanity-check they're present.
    assert(~isempty(pmodel_obj.Y),  'Y empty');
    assert(~isempty(pmodel_obj.id), 'id empty');
    assert(isequal(pmodel_obj.yfit,                           pmodel_obj.fitted_values.yfit),                          'yfit');
    assert(isequal(pmodel_obj.dist_from_hyperplane_xval,      pmodel_obj.fitted_values.dist_from_hyperplane_xval),     'dist_from_hyperplane_xval');
    assert(isequal(pmodel_obj.class_probability_xval,         pmodel_obj.fitted_values.class_probability_xval),        'class_probability_xval');
    assert(isequal(pmodel_obj.scorediff,                      pmodel_obj.fitted_values.scorediff),                     'scorediff');
    assert(isequal(pmodel_obj.w,                              pmodel_obj.weights.w),                                   'w');
    % error_metrics entries are now (value, descrip) tuples; legacy alias unwraps .value.
    assert(isequal(pmodel_obj.crossval_accuracy,              pmodel_obj.error_metrics.crossval_accuracy.value),       'crossval_accuracy');
    assert(isequaln(pmodel_obj.classification_d_singleinterval, pmodel_obj.error_metrics.d_singleinterval.value),      'classification_d_singleinterval');
    assert(isequaln(pmodel_obj.crossval_accuracy_within,      pmodel_obj.error_metrics.crossval_accuracy_within.value),'crossval_accuracy_within');
    assert(isequaln(pmodel_obj.classification_d_within,       pmodel_obj.error_metrics.d_within.value),                'classification_d_within');
    assert(isequal(pmodel_obj.SVMModel,                       pmodel_obj.ml_model),                                    'SVMModel');
    assert(isequal(pmodel_obj.ClassificationModel,            pmodel_obj.ml_model),                                    'ClassificationModel');
    assert(isequal(pmodel_obj.trIdx,                          pmodel_obj.cv_partition.trIdx),                          'trIdx');
    assert(isequal(pmodel_obj.teIdx,                          pmodel_obj.cv_partition.teIdx),                          'teIdx');
    assert(isequal(pmodel_obj.nfolds,                         pmodel_obj.cv_partition.nfolds),                         'nfolds');
    fprintf('  Categorised <-> legacy alias round-trip OK\n');

    % --- Shape / sanity ---
    n = numel(Y);
    assert(numel(pmodel_obj.yfit) == n,                       'yfit length mismatch');
    assert(numel(pmodel_obj.dist_from_hyperplane_xval) == n,  'dist_from_hyperplane_xval length mismatch');
    assert(size(pmodel_obj.w, 1) == size(X, 2),               'w length should match number of features');
    assert(pmodel_obj.crossval_accuracy >= 0 && pmodel_obj.crossval_accuracy <= 100, ...
        'crossval_accuracy out of [0,100]');
    assert(pmodel_obj.diagnostics.mult_obs_within_person == true, ...
        'Should detect multiple obs per id (each subj has Hot+Warm)');
    fprintf('  Shape/sanity OK: cv_acc = %.1f%%, d_single = %.2f, cv_acc_within = %.1f%%, d_within = %.2f\n', ...
        pmodel_obj.crossval_accuracy, pmodel_obj.classification_d_singleinterval, ...
        pmodel_obj.crossval_accuracy_within, pmodel_obj.classification_d_within);

    % --- fit_type + omitted markers (Phase B) ---
    assert(strcmp(pmodel_obj.fit_type, 'crossval'), ...
        'fit_type should be ''crossval'', got ''%s''', pmodel_obj.fit_type);
    assert(islogical(pmodel_obj.omitted_cases),    'omitted_cases must be logical');
    assert(islogical(pmodel_obj.omitted_features), 'omitted_features must be logical');
    assert(numel(pmodel_obj.omitted_cases)    == numel(Y),         'omitted_cases length should match original Y');
    assert(numel(pmodel_obj.omitted_features) == size(X, 2),       'omitted_features length should match original feature count');
    fprintf('  fit_type=%s, omitted_cases=%d, omitted_features=%d\n', ...
        pmodel_obj.fit_type, sum(pmodel_obj.omitted_cases), sum(pmodel_obj.omitted_features));

    % --- validate_object accepts the returned object ---
    pmodel_obj.validate_object('noverbose');
    fprintf('  validate_object OK\n');

    % --- Clone preserves hyperparameters, clears fitted state ---
    pmodel_obj2 = clone(pmodel_obj);
    assert(~pmodel_obj2.is_fitted,                                       'clone did not clear fitted state');
    assert(isequal(pmodel_obj2.modeloptions, pmodel_obj.modeloptions),   'clone lost modeloptions');
    fprintf('  clone() OK\n');

    fprintf('xval_SVM_unit_test: PASS\n');
end
