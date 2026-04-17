# ROI-wise LDA decoding

Custom helper functions for parcel-wise decoding from Canlab contrast objects, added without modifying the core toolbox.

## Main entry point

`canlab_roiwise_lda_decode.m`

## Typical use

```matlab
addpath(genpath('C:\Work\Toolboxes\CanlabCore\custom_tools\roi_lda_decoding'));

results = canlab_roiwise_lda_decode(C_obj, Group_labels, ...
    'atlas', 'canlab2024_coarse_fmriprep20_2mm', ...
    'nFolds', 2, ...
    'leaveStudyOut', false, ...
    'nPermutations', 100, ...
    'useParallel', true, ...
    'parallelMode', 'permutations', ...
    'doplot', false);
```

## Outputs

- `results.roi_table`: one row per ROI with voxel counts, accuracy, `p_value`, and `q_value`
- `results.mean_fold_accuracy_map`: parcel-expanded `fmri_data` object
- `results.overall_accuracy_map`: parcel-expanded `fmri_data` object
- `results.p_value_map`: parcel-expanded `fmri_data` object
- `results.q_value_map`: parcel-expanded `fmri_data` object

## Performance notes

- The function now stores only ROI column indices, not a full copy of each ROI matrix.
- Permutations use a lightweight internal LDA path that returns only accuracies and does not keep model objects in memory.
- `useParallel=true` tries to start a thread-based pool, which is much more memory-friendly than process workers.
- `parallelMode='permutations'` is usually the best starting point.
- `saveModels=true` disables the parallel observed-data pass because model objects are relatively heavy.

## Statistical notes

- `Group_labels` are study-level labels and are expanded to the subject level internally.
- Voxels with any NaNs across subjects, or zero variance within an ROI, are excluded before LDA.
- The default is subject-level stratified CV because some classes may include only one study.
- If `leaveStudyOut=true`, the function will use study-wise grouping only when each class has at least two studies and the problem is binary; otherwise it falls back to subject-level CV.
- Permutation testing always shuffles labels at the study level, then expands them to subjects.
- `q_value` uses Benjamini-Hochberg FDR across ROIs.
