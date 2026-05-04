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
