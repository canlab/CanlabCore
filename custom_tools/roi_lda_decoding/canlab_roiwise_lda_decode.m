function results = canlab_roiwise_lda_decode(C_obj, Group_labels, varargin)
% ROI-wise LDA decoding across subject-level contrast maps.
%
% results = canlab_roiwise_lda_decode(C_obj, Group_labels)
%
% This function:
% 1. Loads or accepts an atlas object
% 2. Aligns the atlas to the contrast-object voxel space
% 3. Builds one subject-by-voxel matrix per atlas ROI
% 4. Runs cross-validated LDA
% 5. Optionally estimates empirical p-values via study-level permutations
% 6. Returns ROI-wise summaries and parcel-expanded maps

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'C_obj', @(x) iscell(x) && ~isempty(x));
addRequired(p, 'Group_labels', @(x) numel(x) == numel(C_obj));
addParameter(p, 'atlas', 'canlab2024_coarse_fmriprep20_2mm', @(x) ischar(x) || isstring(x) || isa(x, 'atlas'));
addParameter(p, 'nFolds', 5, @(x) isnumeric(x) && isscalar(x) && x >= 2);
addParameter(p, 'leaveStudyOut', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'minVoxels', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'nPermutations', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'randomSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'verbose', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'doplot', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'saveModels', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'useParallel', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'parallelMode', 'permutations', @(x) ischar(x) || isstring(x));
parse(p, C_obj, Group_labels, varargin{:});

atlas_input = p.Results.atlas;
nFoldsRequested = double(p.Results.nFolds);
leaveStudyOutRequested = p.Results.leaveStudyOut;
minVoxels = double(p.Results.minVoxels);
nPermutations = double(p.Results.nPermutations);
randomSeed = p.Results.randomSeed;
verbose = p.Results.verbose;
doplot = p.Results.doplot;
saveModels = p.Results.saveModels;
useParallelRequested = p.Results.useParallel;
parallelMode = lower(char(p.Results.parallelMode));

if ~isempty(randomSeed)
    rng(randomSeed);
end

prep = canlab_custom_prepare_roi_lda_data(C_obj, Group_labels, atlas_input, 'verbose', verbose);

classes = unique(prep.subject_labels);
nClasses = numel(classes);
if nClasses < 2
    error('%s: at least two classes are required for decoding.', mfilename);
end

leaveStudyOut = local_resolve_cv_mode(prep, leaveStudyOutRequested, nClasses, verbose);

if leaveStudyOut
    counts_for_folds = accumarray(prep.study_labels_numeric, 1);
else
    counts_for_folds = accumarray(prep.subject_labels, 1);
end

maxAllowedFolds = min(counts_for_folds);
if maxAllowedFolds < 2
    error('%s: each class must contribute at least two units to cross-validation.', mfilename);
end

nFolds = min(nFoldsRequested, maxAllowedFolds);
if verbose && nFolds < nFoldsRequested
    fprintf('%s: using %d folds (requested %d).\n', mfilename, nFolds, nFoldsRequested);
end

observed_cv = local_make_cv_partition(prep.subject_labels, nFolds, leaveStudyOut, prep.study_id_subject);

[useParallel, pool_info] = local_prepare_parallel(useParallelRequested, verbose);
if useParallel && saveModels
    useParallel = false;
    if verbose
        warning('%s: saveModels=true forces serial execution for the observed-data pass.', mfilename);
    end
end

nRegions = numel(prep.region_ids);
mean_fold_accuracy = NaN(nRegions, 1);
overall_accuracy = NaN(nRegions, 1);
p_value = NaN(nRegions, 1);
q_value = NaN(nRegions, 1);
fdr_p_value = NaN(nRegions, 1);
nvox_total = zeros(nRegions, 1);
nvox_used = zeros(nRegions, 1);
status = repmat({'not_run'}, nRegions, 1);
notes = repmat({''}, nRegions, 1);
roi_col_index = cell(nRegions, 1);

if saveModels
    roi_stats = cell(nRegions, 1);
else
    roi_stats = [];
end

for r = 1:nRegions
    roi_mask = prep.atlas_index == prep.region_ids(r);
    Xroi = prep.X(:, roi_mask);
    nvox_total(r) = size(Xroi, 2);

    valid_vox = all(isfinite(Xroi), 1) & (std(Xroi, 0, 1) > 0);
    roi_col_index{r} = find(roi_mask);
    roi_col_index{r} = roi_col_index{r}(valid_vox);
    nvox_used(r) = numel(roi_col_index{r});

    if nvox_used(r) < minVoxels
        status{r} = 'skipped';
        notes{r} = sprintf('Fewer than %d usable voxels after filtering.', minVoxels);
    end
end

ok_candidates = find(strcmp(status, 'not_run'));
ok_candidate_cols = roi_col_index(ok_candidates);

if useParallel && strcmp(parallelMode, 'regions') && ~saveModels
    mean_local = NaN(numel(ok_candidates), 1);
    overall_local = NaN(numel(ok_candidates), 1);
    status_local = repmat({''}, numel(ok_candidates), 1);
    notes_local = repmat({''}, numel(ok_candidates), 1);
    Xfull = prep.X;
    y = prep.subject_labels;
    cv_partition = observed_cv;

    parfor idx = 1:numel(ok_candidates)
        Xroi = Xfull(:, ok_candidate_cols{idx});
        try
            [mean_local(idx), overall_local(idx)] = local_lightweight_lda_accuracy(Xroi, y, cv_partition);
            status_local{idx} = 'ok';
            notes_local{idx} = '';
        catch ME
            status_local{idx} = 'failed';
            notes_local{idx} = ME.message;
        end
    end

    mean_fold_accuracy(ok_candidates) = mean_local;
    overall_accuracy(ok_candidates) = overall_local;
    status(ok_candidates) = status_local;
    notes(ok_candidates) = notes_local;
else
    for idx = 1:numel(ok_candidates)
        r = ok_candidates(idx);
        Xroi = prep.X(:, roi_col_index{r});

        try
            if saveModels
                classify_args = {'nFolds', nFolds, 'verbose', false, 'doplot', doplot};
                if leaveStudyOut
                    classify_args = [classify_args {'id', prep.study_id_subject}];
                end
                S = xval_classify(Xroi, prep.subject_labels, classify_args{:});
                mean_fold_accuracy(r) = mean(S.accuracy);
                overall_accuracy(r) = S.overallAccuracy;
                roi_stats{r} = S;
            else
                [mean_fold_accuracy(r), overall_accuracy(r)] = local_lightweight_lda_accuracy(Xroi, prep.subject_labels, observed_cv);
            end
            status{r} = 'ok';
        catch ME
            status{r} = 'failed';
            notes{r} = ME.message;
        end
    end
end

ok_mask = strcmp(status, 'ok');
ok_regions = find(ok_mask);
ok_region_cols = roi_col_index(ok_regions);

if nPermutations > 0 && any(ok_mask)
    if verbose
        fprintf('%s: running %d study-level permutations for empirical p-values.\n', mfilename, nPermutations);
    end

    null_ge_count = zeros(nRegions, 1);
    study_labels_numeric = prep.study_labels_numeric;
    study_nsubjects = prep.study_nsubjects;
    observed_acc = overall_accuracy(ok_regions);

    if useParallel && strcmp(parallelMode, 'permutations')
        perm_counts = zeros(numel(ok_regions), nPermutations);
        Xfull = prep.X;
        leaveStudyOut_local = leaveStudyOut;
        study_id_local = prep.study_id_subject;

        parfor permIdx = 1:nPermutations
            permuted_study_labels = study_labels_numeric(randperm(numel(study_labels_numeric)));
            permuted_subject_labels = repelem(permuted_study_labels(:), study_nsubjects(:));
            perm_cv = local_make_cv_partition(permuted_subject_labels, nFolds, leaveStudyOut_local, study_id_local);
            perm_hit = zeros(numel(ok_regions), 1);

            for j = 1:numel(ok_regions)
                Xroi = Xfull(:, ok_region_cols{j});
                try
                    [~, perm_acc] = local_lightweight_lda_accuracy(Xroi, permuted_subject_labels, perm_cv);
                    perm_hit(j) = double(perm_acc >= observed_acc(j));
                catch
                end
            end

            perm_counts(:, permIdx) = perm_hit;
        end

        null_ge_count(ok_regions) = sum(perm_counts, 2);
    else
        for permIdx = 1:nPermutations
            permuted_study_labels = study_labels_numeric(randperm(numel(study_labels_numeric)));
            permuted_subject_labels = repelem(permuted_study_labels(:), study_nsubjects(:));
            perm_cv = local_make_cv_partition(permuted_subject_labels, nFolds, leaveStudyOut, prep.study_id_subject);

            for j = 1:numel(ok_regions)
                r = ok_regions(j);
                Xroi = prep.X(:, roi_col_index{r});
                try
                    [~, perm_acc] = local_lightweight_lda_accuracy(Xroi, permuted_subject_labels, perm_cv);
                    null_ge_count(r) = null_ge_count(r) + double(perm_acc >= overall_accuracy(r));
                catch
                end
            end

            if verbose && (permIdx == 1 || mod(permIdx, max(1, floor(nPermutations / 10))) == 0 || permIdx == nPermutations)
                fprintf('  permutation %d/%d complete\n', permIdx, nPermutations);
            end
        end
    end

    p_value(ok_mask) = (1 + null_ge_count(ok_mask)) ./ (1 + nPermutations);
    q_value(ok_mask) = local_bh_fdr(p_value(ok_mask));
    fdr_p_value = q_value;
end

roi_table = table( ...
    prep.region_ids, ...
    prep.region_labels(:), ...
    nvox_total, ...
    nvox_used, ...
    mean_fold_accuracy, ...
    overall_accuracy, ...
    p_value, ...
    q_value, ...
    fdr_p_value, ...
    status, ...
    notes, ...
    'VariableNames', {'region_id', 'region_label', 'nvox_total', 'nvox_used', ...
    'mean_fold_accuracy', 'overall_accuracy', 'p_value', 'q_value', 'fdr_p_value', 'status', 'notes'});

results = struct();
results.atlas_obj = prep.atlas_obj;
results.region_ids = prep.region_ids;
results.region_labels = prep.region_labels;
results.roi_table = roi_table;
results.mean_fold_accuracy = mean_fold_accuracy;
results.overall_accuracy = overall_accuracy;
results.p_value = p_value;
results.q_value = q_value;
results.fdr_p_value = fdr_p_value;
results.mean_fold_accuracy_map = parcel_data2fmri_data(prep.atlas_obj, mean_fold_accuracy);
results.overall_accuracy_map = parcel_data2fmri_data(prep.atlas_obj, overall_accuracy);
results.p_value_map = parcel_data2fmri_data(prep.atlas_obj, p_value);
results.q_value_map = parcel_data2fmri_data(prep.atlas_obj, q_value);
results.fdr_p_value_map = parcel_data2fmri_data(prep.atlas_obj, fdr_p_value);
results.nvox_total = nvox_total;
results.nvox_used = nvox_used;
results.status = status;
results.notes = notes;
results.subject_labels = prep.subject_labels;
results.study_labels_numeric = prep.study_labels_numeric;
results.study_id_subject = prep.study_id_subject;
results.class_names = prep.class_names;
results.leaveStudyOutRequested = leaveStudyOutRequested;
results.leaveStudyOutUsed = leaveStudyOut;
results.nFolds = nFolds;
results.nPermutations = nPermutations;
results.nsubjects_total = prep.nsubjects_total;
results.nvoxels_shared = prep.nvoxels_shared;
results.useParallelRequested = useParallelRequested;
results.useParallelUsed = useParallel;
results.parallelMode = parallelMode;
results.parallelPool = pool_info;
results.observed_cv = observed_cv;

if saveModels
    results.roi_stats = roi_stats;
end

end

function cv_partition = local_make_cv_partition(labels, nFolds, leaveStudyOut, study_id)
if leaveStudyOut
    [trIdxCell, teIdxCell] = xval_stratified_holdout_leave_whole_subject_out(labels, study_id, 'nfolds', nFolds, 'doverbose', false, 'doplot', false);
else
    folds = stratified_holdout_set(labels, nFolds);
    trIdxCell = folds.trIdx;
    teIdxCell = folds.teIdx;
end

cv_partition = struct('trIdx', {trIdxCell}, 'teIdx', {teIdxCell}, 'nFolds', nFolds);
end

function [mean_fold_acc, overall_acc] = local_lightweight_lda_accuracy(X, labels, cv_partition)
labels = labels(:);
yfit = NaN(size(labels));
fold_acc = NaN(cv_partition.nFolds, 1);

for f = 1:cv_partition.nFolds
    trainMask = logical(cv_partition.trIdx{f});
    testMask = logical(cv_partition.teIdx{f});
    model = fitcdiscr(X(trainMask, :), labels(trainMask));
    ypred = predict(model, X(testMask, :));
    yfit(testMask) = ypred;
    fold_acc(f) = 100 * mean(ypred == labels(testMask));
end

mean_fold_acc = mean(fold_acc);
overall_acc = 100 * mean(yfit == labels);
end

function leaveStudyOut = local_resolve_cv_mode(prep, leaveStudyOutRequested, nClasses, verbose)
leaveStudyOut = leaveStudyOutRequested;
if ~leaveStudyOutRequested
    return
end

study_counts = accumarray(prep.study_labels_numeric, 1);

if nClasses > 2
    leaveStudyOut = false;
    if verbose
        warning('%s: multiclass problems use subject-level CV in this function.', mfilename);
    end
    return
end

if any(study_counts < 2)
    leaveStudyOut = false;
    if verbose
        warning('%s: at least one class has fewer than two studies, so using subject-level CV instead.', mfilename);
    end
end
end

function [useParallel, pool_info] = local_prepare_parallel(useParallelRequested, verbose)
useParallel = false;
pool_info = '';

if ~useParallelRequested
    return
end

hasPct = license('test', 'Distrib_Computing_Toolbox');
if ~hasPct
    if verbose
        warning('%s: Parallel Computing Toolbox not available. Using serial execution.', mfilename);
    end
    return
end

pool = gcp('nocreate');
if isempty(pool)
    try
        pool = parpool('threads');
    catch ME
        if verbose
            warning('%s: could not start a thread-based pool (%s). Using serial execution.', mfilename, ME.message);
        end
        return
    end
end

useParallel = true;
pool_info = sprintf('pool with %d workers', pool.NumWorkers);
if verbose
    fprintf('%s: using %s.\n', mfilename, pool_info);
end
end

function q = local_bh_fdr(p)
p = p(:);
q = NaN(size(p));
valid = isfinite(p);
if ~any(valid)
    return
end

pv = p(valid);
[p_sorted, sort_idx] = sort(pv, 'ascend');
m = numel(p_sorted);
rank = (1:m)';
q_sorted = p_sorted .* m ./ rank;
q_sorted = min(1, q_sorted);
q_sorted = flipud(cummin(flipud(q_sorted)));
q_valid = NaN(size(pv));
q_valid(sort_idx) = q_sorted;
q(valid) = q_valid;
end
