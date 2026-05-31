function xval_discriminant_classifier_unit_test()
% xval_discriminant_classifier_unit_test  Smoke test on DPSP Hot vs Warm.
%
% Drives cross-validated LDA via fitcdiscr on the DPSP single-subject
% Hot vs Warm maps. Asserts @predictive_model return and
% categorised <-> legacy alias agreement.

    fprintf('=== xval_discriminant_classifier_unit_test ===\n');

    canlabcore_dir = fileparts(fileparts(which('fmri_data')));
    sample_dir = fullfile(canlabcore_dir, 'Sample_datasets', ...
                          'DPSP_pain_rejection_participant_maps');
    H = load(fullfile(sample_dir, 'DPSP_single_subject_images_hot.mat'));
    W = load(fullfile(sample_dir, 'DPSP_single_subject_images_warm.mat'));
    hot  = H.single_subject_images_hot;
    warm = W.single_subject_images_warm;

    % Stratify: 100 random voxels (LDA needs n > p)
    rng(0);
    p = size(hot.dat, 1);
    keep_vox = randsample(p, 100);
    X = double([hot.dat(keep_vox, :) warm.dat(keep_vox, :)])';
    labels = int32([ones(size(hot.dat,  2), 1); ...
                    2*ones(size(warm.dat, 2), 1)]);

    pmodel_obj = xval_discriminant_classifier(X, labels, 'nFolds', 5, ...
                                              'verbose', false, 'doplot', false);

    assert(isa(pmodel_obj, 'predictive_model'), ...
        'Expected predictive_model, got %s', class(pmodel_obj));
    assert(pmodel_obj.is_fitted, 'is_fitted false');
    fprintf('  class = %s, is_fitted = %d\n', class(pmodel_obj), pmodel_obj.is_fitted);

    % Categorised <-> legacy alias agreement
    assert(~isempty(pmodel_obj.Y), 'Y empty');
    assert(isequal(pmodel_obj.yfit,              pmodel_obj.fitted_values.yfit),                'yfit');
    assert(isequal(pmodel_obj.predictions,       pmodel_obj.fitted_values.predictions),         'predictions');
    % trueLabels was per-fold cell, now routed to fitted_values.Y_per_fold
    assert(isequal(pmodel_obj.trueLabels,        pmodel_obj.fitted_values.Y_per_fold),          'trueLabels');
    % error_metrics entries are (value, descrip) tuples.
    assert(isequal(pmodel_obj.accuracy,          pmodel_obj.error_metrics.accuracy.value),        'accuracy');
    assert(isequal(pmodel_obj.overallAccuracy,   pmodel_obj.error_metrics.overallAccuracy.value), 'overallAccuracy');
    assert(isequal(pmodel_obj.trIdx,             pmodel_obj.cv_partition.trIdx),                'trIdx');
    assert(isequal(pmodel_obj.teIdx,             pmodel_obj.cv_partition.teIdx),                'teIdx');
    assert(isequal(pmodel_obj.nfolds,            pmodel_obj.cv_partition.nfolds),               'nfolds');
    % `models` is routed to fold_models (cell of per-fold trained classifiers)
    assert(numel(pmodel_obj.fold_models) == pmodel_obj.nfolds, ...
        'fold_models length should equal nfolds');
    fprintf('  Categorised <-> legacy alias round-trip OK\n');

    fprintf('  overall_acc = %.1f%%, mean_fold_acc = %.1f%%\n', ...
        pmodel_obj.overallAccuracy, mean(pmodel_obj.accuracy));

    % --- fit_type + omitted markers (Phase B) ---
    assert(strcmp(pmodel_obj.fit_type, 'crossval'), 'fit_type should be crossval');
    assert(islogical(pmodel_obj.omitted_cases),    'omitted_cases must be logical');
    assert(islogical(pmodel_obj.omitted_features), 'omitted_features must be logical');
    fprintf('  fit_type=%s, omitted_cases=%d, omitted_features=%d\n', ...
        pmodel_obj.fit_type, sum(pmodel_obj.omitted_cases), sum(pmodel_obj.omitted_features));

    pmodel_obj.validate_object('noverbose');
    fprintf('xval_discriminant_classifier_unit_test: PASS\n');
end
