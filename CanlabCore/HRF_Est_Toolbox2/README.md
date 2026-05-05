# HRF Estimation Toolbox 2: Simple End-to-End Pipeline

This folder now includes a user-facing pipeline that goes from a **4D fMRI NIfTI** and **BIDS events.tsv** directly to model-based HRF estimates.

## Quick start

```matlab
addpath(genpath('/path/to/CanlabCore/CanlabCore/HRF_Est_Toolbox2'));

results = run_hrf_pipeline( ...
    '/path/to/sub-01_task-pain_bold.nii.gz', ...
    '/path/to/sub-01_task-pain_events.tsv', ...
    'TR', 0.8, ...
    'MaskNii', '/path/to/brain_mask.nii.gz', ...
    'Conditions', {'pain', 'neutral'}, ...
    'WindowSeconds', 30, ...
    'Models', {'logit', 'sfir', 'canonical'}, ...
    'OutputMat', '/path/to/sub-01_hrf_results.mat');
```

## What this pipeline does

1. Loads a 4D fMRI volume and extracts a mean timeseries (whole brain or mask).
2. Loads and validates events from BIDS `events.tsv`.
3. Builds one stick function per condition (`trial_type`).
4. Runs selected HRF models using existing toolbox fitters.
5. Returns a single `results` struct (and can save `.mat`).

## New API

- `run_hrf_pipeline.m` — main entry point.
- `hrf_extract_timeseries_from_nii.m` — 4D NIfTI to z-scored timeseries.
- `hrf_load_events_tsv.m` — BIDS events loader/validator.
- `hrf_build_stick_functions.m` — converts events to condition-wise sticks.
- `hrf_fit_all_models.m` — unified interface to model fitting.

## Notes

- `events.tsv` must include: `onset`, `duration`, `trial_type`.
- If TR is not passed, the pipeline attempts to use NIfTI header `PixelDimensions(4)`.
- Model wrappers call legacy toolbox methods (`Fit_Logit2`, `Fit_sFIR`, `Fit_Canonical_HRF`, `Fit_Spline`, `Fit_NLgamma`).

## Trial averaging + condition comparison

```matlab
% Run full pipeline first
results = run_hrf_pipeline(fmri_nii, events_tsv, 'TR', 0.8);

% Average one condition across trials (20 s window)
pain = hrf_average_condition_trials(results.timeseries, results.events, 'pain', 0.8, 20, 'BaselineSeconds', 2);
neutral = hrf_average_condition_trials(results.timeseries, results.events, 'neutral', 0.8, 20, 'BaselineSeconds', 2);

% Compare conditions
cmp = hrf_compare_conditions(pain, neutral, 'Alpha', 0.05);

% Plot means with SEM
figure;
subplot(1,2,1); hold on;
plot(pain.time, pain.mean, 'r', 'LineWidth', 1.5);
plot(neutral.time, neutral.mean, 'b', 'LineWidth', 1.5);
legend({'pain','neutral'}); xlabel('Seconds'); ylabel('Signal');
title('Condition-averaged responses');

subplot(1,2,2);
plot(cmp.time, cmp.mean_diff, 'k', 'LineWidth', 1.5); hold on;
yline(0,'--');
xlabel('Seconds'); ylabel('Pain - Neutral');
title(sprintf('Mean difference (AUC diff = %.3f)', cmp.auc.diff));
```


## Signature-based (interpretable) time-series

You can now estimate HRFs from a signature-expression time series (via `apply_all_signatures`) instead of a mean BOLD time series:

```matlab
results_sig = run_hrf_pipeline( ...
    fmri_nii, events_tsv, ...
    'TR', 0.8, ...
    'SignalSource', 'signature', ...
    'SimilarityMetric', 'dotproduct', ...
    'ImageSet', 'all', ...
    'SignatureName', 'NPS');

disp(results_sig.signature_meta)
```

If `SignatureName` is omitted, the first signature from `apply_all_signatures` is used.
Note: if you pass `SignatureName` while `SignalSource` is left at default (`'mean'`), the pipeline now auto-switches to signature mode and warns you.


## Quick plotting helper for new results structure

```matlab
% Basic: one model, selected conditions
plot_hrf_results(results, 'Model', 'sfir', 'Conditions', [4 9 10 11]);
plot_hrf_results(results, 'Model', 'sfir', 'Conditions', {'pain','neutral'});

% Signature-specific (when SignalSource='signature')
plot_hrf_results(results_sig, 'Model', 'sfir', 'Signature', 'NPS', 'Conditions', [1 2 3]);
```

`results_sig.signature_meta` now includes selected and available signature names, and
`results_sig.fits_by_signature` stores fitted models for each signature.


`plot_hrf_results` now supports condition names (cellstr) and draws standard-error style shading when multiple signature fits are available.

## Speeding up all-signature fitting

If fitting all signatures is slow, you can:

```matlab
results_sig = run_hrf_pipeline( ...
    fmri_nii, events_tsv, ...
    'TR', 0.8, ...
    'SignalSource', 'signature', ...
    'ImageSet', 'all', ...
    'UseParallel', true, ...           % parfor across signatures
    'MaxSignatures', 10, ...           % only fit first 10 signatures
    'SignatureNames', {'NPS','SIIPS'}  % or fit an explicit subset
    );
```

For fastest runs, pass a single `SignatureName` (e.g., `'NPS'`) so the pipeline skips full all-signature loading/fitting.

## Apply canonical brain map sets (Buckner, Margulies, Hansen)

You can use `SignalSource='imageset'` to apply maps from `load_image_set` directly:

```matlab
% Buckner 7-network maps
res_buckner = run_hrf_pipeline(fmri_nii, events_tsv, ...
    'SignalSource', 'imageset', 'ImageSet', 'bucknerlab_wholebrain', ...
    'TR', 0.8, 'UseParallel', true);

% Margulies principal gradient
res_marg = run_hrf_pipeline(fmri_nii, events_tsv, ...
    'SignalSource', 'imageset', 'ImageSet', 'marg', 'TR', 0.8);

% Hansen receptor maps (subset)
res_hansen = run_hrf_pipeline(fmri_nii, events_tsv, ...
    'SignalSource', 'imageset', 'ImageSet', 'hansen22', ...
    'MapNames', {'D1', '5HT1a'}, 'TR', 0.8);
```

## Study-level pipeline (multiple images/subjects)

```matlab
study = run_hrf_study_pipeline(fmri_files, events_files, subject_ids, ...
    'TR', 0.8, 'SignalSource', 'imageset', 'ImageSet', 'bucknerlab_wholebrain');

% Long-format subject summary table
study.summary

% Visualize separated by subject
plot_hrf_study_by_subject(study, 'Model', 'sfir', 'Condition', 'pain');
```

## Multilevel time-unfolding significance testing

```matlab
study = run_hrf_study_pipeline(fmri_files, events_files, subject_ids, ...
    'TR', 0.8, 'SignalSource', 'mean');

% Subject-level mean extracted timeseries (subjects x time)
study.mean_timeseries

% Within-subject condition contrast, then group-level test across time bins
stats = hrf_time_unfolding_stats(study, ...
    'Model', 'sfir', ...
    'ConditionA', 'pain', ...
    'ConditionB', 'neutral', ...
    'Alpha', 0.05);

plot_hrf_time_unfolding_stats(stats);
```

Optional between-group testing is available by supplying `Group` labels (two-group t-test per time bin).

## Efficient ROI- and pattern-based HRF fitting

Yes—this is often more efficient than voxelwise fitting. You can now do:

```matlab
% 1) Atlas ROI means (SignalSource='atlas')
at = load_atlas('canlab2018');
res_roi = run_hrf_pipeline(fmri_nii, events_tsv, ...
    'TR', 0.8, 'SignalSource', 'atlas', ...
    'AtlasObj', at, 'Regions', {'ACC','vmPFC'});

% 2) Pattern/map expression (SignalSource='imageset')
res_maps = run_hrf_pipeline(fmri_nii, events_tsv, ...
    'TR', 0.8, 'SignalSource', 'imageset', ...
    'ImageSet', 'bucknerlab_wholebrain', 'UseParallel', true);
```

## Fast montage animation over time (betas or thresholded t maps)

```matlab
% Animate a 4D image quickly by skipping frames
hrf_make_montage_animation(fmri_nii, 'timeseries.mp4', ...
    'FrameStep', 3, 'FPS', 10, 'TitlePrefix', 'BOLD');

% Thresholded t-map animation (if your img is already t-values)
hrf_make_montage_animation(tmap_4d, 'tvals_thresh.mp4', ...
    'Threshold', 2.0, 'FrameStep', 2, 'FPS', 12, 'TitlePrefix', 't-value');
```
