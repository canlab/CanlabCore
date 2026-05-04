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

## Quick plotting helper for new results structure

```matlab
% Basic: one model, selected conditions
plot_hrf_results(results, 'Model', 'sfir', 'Conditions', [4 9 10 11]);

% Signature-specific (when SignalSource='signature')
plot_hrf_results(results_sig, 'Model', 'sfir', 'Signature', 'NPS', 'Conditions', [1 2 3]);
```

`results_sig.signature_meta` now includes selected and available signature names, and
`results_sig.fits_by_signature` stores fitted models for each signature.
