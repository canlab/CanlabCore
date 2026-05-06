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

- `run_hrf_pipeline.m` - main entry point.
- `hrf_extract_timeseries_from_nii.m` - 4D NIfTI to z-scored timeseries.
- `hrf_load_events_tsv.m` - BIDS events loader/validator.
- `hrf_build_stick_functions.m` - converts events to condition-wise sticks.
- `hrf_fit_all_models.m` - unified interface to model fitting.
- `hrf_write_slurm_study_script.m` - writes a SLURM array script, manifest, and MATLAB worker for study-wide whole-brain HRF fitting.

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

Use `'UseParallelSubjects', true` to run subjects in a local MATLAB parallel pool. For
clusters, the same per-file call is suitable for a SLURM array job.
Use `'WholeBrainOutputDir', '/path/to/hrf_outputs'` with `'WriteWholeBrain', true`
to write subject-specific beta/T outputs into one second-level input directory.

## Study-wide whole-brain maps plus all signatures/imagesets

For a full study, write the 4D beta/T maps once per subject and then apply all
requested signatures and image sets to those 4D maps. This is usually faster
than refitting or reloading signatures separately for each summary.

```matlab
output_dir = '/path/to/hrf_outputs';

image_sets = { ...
    'bucknerlab_wholebrain', ...
    'marg', ...
    'hansen22'};

study = run_hrf_study_pipeline(fmri_files, events_files, subject_ids, ...
    'TR', 0.8, ...
    'Conditions', {'pain', 'neutral'}, ...
    'WindowSeconds', 30, ...
    'WriteWholeBrain', true, ...
    'WholeBrainOutputDir', output_dir, ...
    'WholeBrainPThresh', 0.005, ...
    'WholeBrainThreshType', 'unc', ...
    'UseParallelSubjects', true);

for i = 1:numel(study.results)
    if ~study.success(i), continue; end

    prefix = fullfile(output_dir, [regexprep(subject_ids{i}, '[^\w.-]', '_') '_hrf']);

    hrf_apply_maps_to_wholebrain(study.results{i}.wholebrain, ...
        'Object', 'beta', ...
        'SignatureSets', {'all'}, ...
        'ImageSets', image_sets, ...
        'SimilarityMetric', 'dotproduct', ...
        'OutputCsv', [prefix '_beta_map_scores.csv']);

    hrf_apply_maps_to_wholebrain(study.results{i}.wholebrain, ...
        'Object', 't', ...
        'SignatureSets', {'all'}, ...
        'ImageSets', image_sets, ...
        'SimilarityMetric', 'dotproduct', ...
        'OutputCsv', [prefix '_t_map_scores.csv']);
end

second_level_inputs = hrf_collect_wholebrain_outputs(output_dir, ...
    'OutputCsv', fullfile(output_dir, 'second_level_inputs.csv'));
```

For a cluster, generate a SLURM array job that runs the same per-subject work.
Each array task writes beta/T maps, metadata, and beta/T map-score CSVs:

```matlab
slurm_paths = hrf_write_slurm_study_script(fmri_files, events_files, subject_ids, ...
    'OutputDir', '/path/to/hrf_outputs', ...
    'CanlabRoot', '/path/to/CanlabCore/CanlabCore', ...
    'PipelineArgs', { ...
        'TR', 0.8, ...
        'Conditions', {'pain', 'neutral'}, ...
        'WindowSeconds', 30, ...
        'WholeBrainPThresh', 0.005, ...
        'WholeBrainThreshType', 'unc'}, ...
    'SignatureSets', {'all'}, ...
    'ImageSets', image_sets, ...
    'ScoreObjects', {'beta', 't'}, ...
    'JobName', 'hrf_wholebrain', ...
    'Time', '24:00:00', ...
    'Mem', '32G', ...
    'CpusPerTask', 1, ...
    'ModuleLoad', 'matlab');
```

Submit the generated script from the cluster shell:

```bash
sbatch /path/to/hrf_outputs/run_hrf_study.sbatch
```

If you update the subject list or generator options, rerun
`hrf_write_slurm_study_script` before submitting so the manifest, worker,
config `.mat`, and `.sbatch` file stay in sync. The generated `.sbatch`
also defaults `SLURM_ARRAY_TASK_ID` to `1`, so you can run it with
`bash run_hrf_study.sbatch` inside an interactive allocation for a quick
task-1 smoke test.

After the array completes:

```matlab
second_level_inputs = hrf_collect_wholebrain_outputs('/path/to/hrf_outputs', ...
    'OutputCsv', '/path/to/hrf_outputs/second_level_inputs.csv');
```

## Whole-brain 4D HRF beta and T maps

```matlab
results = run_hrf_pipeline(fmri_nii, events_tsv, ...
    'TR', 0.8, ...
    'WindowSeconds', 30, ...
    'WriteWholeBrain', true, ...
    'WholeBrainOutputPrefix', '/path/to/sub-01_hrf', ...
    'WholeBrainPThresh', 0.005, ...
    'WholeBrainThreshType', 'unc');

% Written files:
%   /path/to/sub-01_hrf_beta.nii
%   /path/to/sub-01_hrf_t.nii
%   /path/to/sub-01_hrf_t_thresh.nii
%   /path/to/sub-01_hrf_metadata.csv
```

The whole-brain output is also available in memory:

```matlab
beta_obj = results.wholebrain.b;  % statistic_image, .dat = beta maps
t_obj = results.wholebrain.t;     % statistic_image, .dat = T maps, .p/.ste/.sig set
```

Apply signatures or image sets after writing/fitting the 4D maps:

```matlab
scores = hrf_apply_maps_to_wholebrain(results.wholebrain, ...
    'Object', 'beta', ...
    'SignatureSets', {'all'}, ...
    'ImageSets', {'bucknerlab_wholebrain', 'hansen22'}, ...
    'SimilarityMetric', 'dotproduct', ...
    'OutputCsv', '/path/to/sub-01_hrf_map_scores.csv');
```

Make a quick condition-lag animation with CANlab montage:

```matlab
hrf_animate_wholebrain_stats(results.wholebrain, ...
    'Object', 't', ...
    'Condition', 'pain', ...
    'OutputFile', '/path/to/sub-01_pain_tmaps.mp4');
```

Collect written files for second-level scripts:

```matlab
input_table = hrf_collect_wholebrain_outputs('/path/to/hrf_outputs', ...
    'OutputCsv', '/path/to/hrf_outputs/second_level_inputs.csv');
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
