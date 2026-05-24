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
    'Models', {'fir', 'sfir', 'canonical'}, ...
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
- `hrf_resolve_condition_patterns.m` - shared exact/wildcard/regex condition matching.
- `hrf_fit_all_models.m` - unified interface to model fitting.
- `plot_hrf_by_condition.m` - single subject/run plotting by condition with explicit source/model labels and fit SE when available.
- `hrf_write_slurm_study_script.m` - writes a SLURM array script, manifest, and MATLAB worker for study-wide whole-brain HRF fitting.
- `hrf_load_wholebrain_stats.m` - rebuilds beta/T `statistic_image` objects from written NIfTI + metadata sidecars.
- `hrf_analyze_second_level_inputs.m` - analyzes signature/imageset score CSVs across subjects from `second_level_inputs.csv`.
- `hrf_second_level_inputs_to_study.m` - converts collected beta/T map-score CSVs into a study-like structure for `hrf_time_unfolding_stats`.
- `hrf_input_table_to_study.m` - rebuilds a study/results structure from `second_level_inputs.csv`, including whole-brain objects and optional map-score curves.

## Notes

- `events.tsv` must include: `onset`, `duration`, `trial_type`.
- If TR is not passed, the pipeline attempts to use NIfTI header `PixelDimensions(4)`.
- Model wrappers call legacy toolbox methods (`Fit_Logit2`, `Fit_sFIR`, `Fit_Canonical_HRF`, `Fit_Spline`, `Fit_NLgamma`). `fir` is the unsmoothed `Fit_sFIR(..., mode=0)` fit; `sfir` is the smoothed `Fit_sFIR(..., mode=1)` fit.
- `canonical` and `nlgamma` require SPM on the MATLAB path; `spline` requires the FDA package (`create_bspline_basis`, `eval_basis`). By default, unavailable optional models are skipped with a warning. Use `'ModelDependencyPolicy', 'error'` to fail fast instead.
- Subject/run-level fit uncertainty is stored in `.se`, `.t`, `.p`, `.dfe`,
  and `.uncertainty_source` for linear fits (`fir`, `sfir`, `canonical`, `spline`).
  Nonlinear fits (`logit`, `nlgamma`) currently store the curve but not a
  model-based SE/p curve.

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

## Condition Wildcards And Regex

Any HRF toolbox function with `Condition`, `Conditions`, `ConditionA`, or
`ConditionB` accepts exact names, wildcards, or regex patterns. Exact names are
preferred when they exist. Wildcards use `*` and `?`; regex can be written as
`'regex:<pattern>'`, `'/<pattern>/'`, or a pattern containing common regex
operators such as `.*`, `|`, `[]`, `()`, `^`, `$`, or `+`.

When one condition pattern matches multiple condition names, those conditions
are averaged and the plotted label records what was averaged:

```matlab
% Average pain_hot and pain_warm into one fitted curve.
plot_hrf_by_condition(results, ...
    'Model', 'sfir', ...
    'Conditions', {'pain_*', 'neutral'});

% Same idea with regex.
stats = hrf_time_unfolding_stats(study, ...
    'Model', 'sfir', ...
    'ConditionA', 'regex:^pain_', ...
    'ConditionB', 'neutral');

% Build a model with grouped event sticks from the beginning.
results = run_hrf_pipeline(fmri_nii, events_tsv, ...
    'Conditions', {'pain_*', 'neutral'}, ...
    'Models', {'sfir'});
```

`results.condition_groups` stores the original condition names that contributed
to each grouped condition.


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
% Single subject/run fitted curves by condition.
% For fir/sfir/canonical/spline fits, ribbons are model-based SE from the run.
plot_hrf_by_condition(results, 'Model', 'sfir', 'Conditions', {'pain','neutral'});

% Compare fitted model families on the same axes.
plot_hrf_by_condition(results, ...
    'Model', {'fir','sfir','canonical'}, ...
    'Condition', 'pain');

% Raw event-locked trial means from the selected 1D time series.
% Ribbons are SEM across repeated events/trials within the run.
plot_hrf_by_condition(results, ...
    'PlotType', 'trialmean', ...
    'Conditions', {'pain','neutral'}, ...
    'BaselineSeconds', 2);

% Signature-specific fitted curves (when SignalSource='signature')
plot_hrf_by_condition(results_sig, ...
    'Model', 'sfir', ...
    'Signature', 'NPS', ...
    'Conditions', {'pain','neutral'});

% Across a study, one curve per subject for a selected condition.
% Ribbons are within-run/within-subject SE when available; the black
% group curve shows mean +/- SEM across plotted subjects or runs.
plot_hrf_study_by_subject(study, ...
    'Model', 'sfir', ...
    'Condition', 'pain', ...
    'Unit', 'subject');

% Compare study-level group means across fitted models.
plot_hrf_study_by_subject(study, ...
    'Model', {'fir','sfir','canonical'}, ...
    'Condition', 'pain', ...
    'Unit', 'subject');

% Across a study, event-locked trial means use repeated event instances
% within each run for SEM at each TR.
plot_hrf_study_by_subject(study, ...
    'PlotType', 'trialmean', ...
    'Condition', 'pain', ...
    'Unit', 'subject');
```

`results_sig.signature_meta` now includes selected and available signature names, and
`results_sig.fits_by_signature` stores fitted models for each signature.

`plot_hrf_results` is kept as a backward-compatible wrapper around
`plot_hrf_by_condition`. It no longer uses across-signature variability as a
standard-error ribbon because that is not subject-level fit uncertainty. Plot
titles now explicitly label the curve source, selected model
(`fir`, `sfir`, `canonical`, `nlgamma`, etc.), and whether the ribbon is within-run
model SE, trial SEM, or unavailable.

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
Use `'ReuseWholeBrainOutputs', true` to reuse existing
`*_beta.nii`/`*_t.nii` files in that directory instead of refitting the
whole-brain 4D maps.

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
    'Models', {'sfir','canonical','spline'}, ...
    'WriteWholeBrain', true, ...
    'WholeBrainOutputDir', output_dir, ...
    'ReuseWholeBrainOutputs', true, ...
    'WholeBrainPThresh', 0.005, ...
    'WholeBrainThreshType', 'unc', ...
    'UseParallelSubjects', true);

for i = 1:numel(study.results)
    if ~study.success(i), continue; end

    prefix = fullfile(output_dir, [regexprep(subject_ids{i}, '[^\w.-]', '_') '_hrf']);
    model_names = fieldnames(study.results{i}.wholebrain_by_model);

    for m = 1:numel(model_names)
        model_name = model_names{m};
        wholebrain = study.results{i}.wholebrain_by_model.(model_name);

        hrf_apply_maps_to_wholebrain(wholebrain, ...
            'Object', 'beta', ...
            'SignatureSets', {'all'}, ...
            'ImageSets', image_sets, ...
            'SimilarityMetric', 'dotproduct', ...
            'OutputCsv', sprintf('%s_%s_beta_map_scores.csv', prefix, model_name));

        hrf_apply_maps_to_wholebrain(wholebrain, ...
            'Object', 't', ...
            'SignatureSets', {'all'}, ...
            'ImageSets', image_sets, ...
            'SimilarityMetric', 'dotproduct', ...
            'OutputCsv', sprintf('%s_%s_t_map_scores.csv', prefix, model_name));
    end
end

second_level_inputs = hrf_collect_wholebrain_outputs(output_dir, ...
    'OutputCsv', fullfile(output_dir, 'second_level_inputs.csv'));
```

For a cluster, generate a SLURM array job that runs the same per-subject work.
Each array task writes beta/T/SE/P maps, metadata, and beta/T map-score CSVs:

```matlab
slurm_paths = hrf_write_slurm_study_script(fmri_files, events_files, subject_ids, ...
    'OutputDir', '/path/to/hrf_outputs', ...
    'CanlabRoot', '/path/to/CanlabCore/CanlabCore', ...
    'PipelineArgs', { ...
        'TR', 0.8, ...
        'Conditions', {'pain', 'neutral'}, ...
        'WindowSeconds', 30, ...
        'Models', {'sfir','canonical','spline'}, ...
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

If `subject_ids` contains repeated subjects, the SLURM writer derives a run
label from each BIDS fMRI filename (`ses-*`, `task-*`, `run-*`, etc.) and
writes run-unique prefixes such as
`SID000002_ses-20_task-distractmap_run-01_hrf`. This prevents multiple runs
from overwriting the same `SID000002_hrf_*` files. To supply your own labels,
pass `'RunLabels', run_labels`.

To run `canonical` or `nlgamma`, add SPM to the worker path; to run `spline`,
also add FDA:

```matlab
slurm_paths = hrf_write_slurm_study_script(fmri_files, events_files, subject_ids, ...
    'OutputDir', '/path/to/hrf_outputs', ...
    'CanlabRoot', '/path/to/CanlabCore/CanlabCore', ...
    'ExtraMatlabPaths', {'/path/to/spm12', '/path/to/fda'}, ...
    'PipelineArgs', {'TR', 0.8, 'Models', {'canonical', 'spline', 'nlgamma'}, ...
        'ModelDependencyPolicy', 'error'}, ...
    'ModuleLoad', 'matlab');
```

`Models` controls the 1D fitters saved in each result MAT file and, when
whole-brain writing is enabled with the default `'WholeBrainMode','auto'`,
also controls which linear 4D whole-brain map families are written. Supported
whole-brain models are `fir`, `sfir`, `canonical`, and `spline`; nonlinear
`logit` and `nlgamma` remain 1D-only. A request such as
`'Models', {'sfir','canonical','spline'}` writes separate files such as
`sub-01_hrf_sfir_beta.nii`, `sub-01_hrf_canonical_beta.nii`, and
`sub-01_hrf_spline_beta.nii`, with matching `_t`, `_se`, `_p`, metadata, and
map-score CSV files. To force only one map family, pass
`'WholeBrainMode','FIR'`, `'sFIR'`, `'canonical'`, or `'spline'` explicitly.

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

`run_hrf_pipeline` saves result MAT files with `-v7.3` by default. When
whole-brain maps are written, the MAT file stores paths and metadata instead
of embedding the large 4D `statistic_image` objects unless
`'SaveWholeBrainInMat', true` is supplied. This keeps cluster result MAT files
loadable while the NIfTI outputs remain the source of truth for 4D maps.

After the array completes:

```matlab
second_level_inputs = hrf_collect_wholebrain_outputs('/path/to/hrf_outputs', ...
    'OutputCsv', '/path/to/hrf_outputs/second_level_inputs.csv');

analysis = hrf_analyze_second_level_inputs(second_level_inputs, ...
    'Object', 'beta', ...
    'ConditionA', 'pain', ...
    'ConditionB', 'neutral', ...
    'OutputSummaryCsv', '/path/to/hrf_outputs/beta_score_group_summary.csv');

% Strongest significant signature/map effects across subjects
analysis.interpretation
```

## Whole-brain 4D HRF beta and T maps

```matlab
results = run_hrf_pipeline(fmri_nii, events_tsv, ...
    'TR', 0.8, ...
    'WindowSeconds', 30, ...
    'Models', {'sfir','canonical','spline'}, ...
    'WriteWholeBrain', true, ...
    'WholeBrainOutputPrefix', '/path/to/sub-01_hrf', ...
    'WholeBrainPThresh', 0.005, ...
    'WholeBrainThreshType', 'unc');

% Written files:
%   /path/to/sub-01_hrf_sfir_beta.nii
%   /path/to/sub-01_hrf_sfir_t.nii
%   /path/to/sub-01_hrf_sfir_se.nii
%   /path/to/sub-01_hrf_sfir_p.nii
%   /path/to/sub-01_hrf_sfir_t_thresh.nii
%   /path/to/sub-01_hrf_sfir_metadata.csv
%   ...and matching canonical/spline files
```

The whole-brain outputs are also available in memory:

```matlab
beta_obj = results.wholebrain_by_model.sfir.b;  % statistic_image beta maps
t_obj = results.wholebrain_by_model.sfir.t;     % statistic_image T maps, .p/.ste/.sig set
```

`Models` and `WholeBrainMode` are related but not identical. `Models` controls
the 1D curve fitters (`logit`, `fir`, `sfir`, `canonical`, `spline`,
`nlgamma`) used for mean, signature, image-set, or atlas time series.
`WholeBrainMode='auto'` writes one 4D condition-lag map set for each requested
supported linear model (`fir`, `sfir`, `canonical`, `spline`). The written
metadata table records the actual model/mode for each volume, and
`hrf_collect_wholebrain_outputs` includes a `model` column.

Existing whole-brain files are overwritten by default
(`'WholeBrainOverwrite', true`). Set `'ReuseWholeBrainOutput', true` to load
existing `_beta.nii`/`_t.nii` outputs instead of recomputing, or
`'WholeBrainOverwrite', false` to stop if outputs already exist.

You can also reconstruct the same object structure later without the result
MAT file:

```matlab
wholebrain = hrf_load_wholebrain_stats('/path/to/sub-01_hrf_sfir');

beta_obj = wholebrain.b;  % .ste and .p restored from _se.nii/_p.nii
t_obj = wholebrain.t;
```

Apply signatures or image sets after writing/fitting the 4D maps:

```matlab
scores = hrf_apply_maps_to_wholebrain(results.wholebrain_by_model.sfir, ...
    'Object', 'beta', ...
    'SignatureSets', {'all'}, ...
    'ImageSets', {'bucknerlab_wholebrain', 'hansen22'}, ...
    'SimilarityMetric', 'dotproduct', ...
    'OutputCsv', '/path/to/sub-01_hrf_sfir_beta_map_scores.csv');
```

Make a quick condition-lag animation with CANlab montage:

```matlab
hrf_animate_wholebrain_stats(results.wholebrain_by_model.sfir, ...
    'Object', 't', ...
    'Condition', 'pain', ...
    'OutputFile', '/path/to/sub-01_pain_tmaps.mp4');
```

Collect written files for second-level scripts:

```matlab
input_table = hrf_collect_wholebrain_outputs('/path/to/hrf_outputs', ...
    'OutputCsv', '/path/to/hrf_outputs/second_level_inputs.csv');
```

Rebuild a `study.results` structure from the collected file index:

```matlab
study = hrf_input_table_to_study(input_table);

% One subject/run/model in the familiar run_hrf_pipeline-style shape
results = study.results{1};
wholebrain = results.wholebrain;
```

If the SLURM/result MAT files are available and you want event-level trial
SEM plots, load those MAT files into the study too:

```matlab
study_runs = hrf_input_table_to_study(input_table, ...
    'LoadWholeBrain', false, ...
    'IncludeMapScores', false, ...
    'LoadResultMat', true);

plot_hrf_study_by_subject(study_runs, ...
    'PlotType', 'trialmean', ...
    'Condition', 'nback-stimblock', ...
    'Unit', 'subject');
```

Use the reconstructed whole-brain object for animation:

```matlab
hrf_animate_wholebrain_stats(study.results{1}.wholebrain, ...
    'Object', 't', ...
    'Condition', 'pain', ...
    'OutputFile', '/path/to/sub-01_pain_tmaps.mp4');
```

Use the reconstructed whole-brain object to apply signatures or image sets.
This writes score CSVs that can later be plotted or analyzed across subjects:

```matlab
for i = 1:numel(study.results)
    if ~study.wholebrain_success(i), continue; end

    prefix = input_table.prefix{i};
    hrf_apply_maps_to_wholebrain(study.results{i}.wholebrain, ...
        'Object', 'beta', ...
        'SignatureSets', {'all'}, ...
        'ImageSets', {'bucknerlab_wholebrain'}, ...
        'SimilarityMetric', 'dotproduct', ...
        'OutputCsv', [prefix '_beta_map_scores.csv']);
end
```

Analyze signature/imageset map scores across subjects:

```matlab
analysis = hrf_analyze_second_level_inputs(input_table, ...
    'Object', 'beta', ...
    'ConditionA', 'pain', ...
    'ConditionB', 'neutral', ...
    'LagSeconds', 6, ...
    'Alpha', 0.05);

analysis.subject_table   % one row per subject x lag x signature/map
analysis.summary         % group mean/SEM/t/p per lag x signature/map
analysis.interpretation  % strongest significant effects
```

Use map-score CSVs with `hrf_time_unfolding_stats`:

```matlab
input_table = hrf_collect_wholebrain_outputs('/path/to/hrf_outputs');

% Backfill missing canonical/sfir/spline map-score CSVs after an older run.
input_table = hrf_score_wholebrain_input_table(input_table, ...
    'SourceModel', 'sfir', ...
    'ScoreObjects', {'beta','t'}, ...
    'SignatureSets', {'all'}, ...
    'ImageSets', image_sets, ...
    'OutputCsv', '/path/to/hrf_outputs/second_level_inputs.csv');

study_scores = hrf_input_table_to_study(input_table, ...
    'LoadWholeBrain', false, ...
    'Object', 'beta', ...
    'SourceModel', 'sfir');

% Default Signature is mean_mapscore, the average over score columns.
% Model='sfir' selects rows whose source whole-brain model was sfir.
plot_hrf_study_by_subject(study_scores, ...
    'Model', 'sfir', ...
    'Condition', 'pain', ...
    'Unit', 'subject');

% Plot one specific score column instead of the average.
plot_hrf_study_by_subject(study_scores, ...
    'Model', 'sfir', ...
    'Signature', 'sig_all_NPS', ...
    'Condition', 'pain', ...
    'Unit', 'subject');

stats_nps = hrf_time_unfolding_stats(study_scores, ...
    'Model', 'sfir', ...
    'Signature', 'sig_all_NPS', ...
    'ConditionA', 'pain', ...
    'ConditionB', 'neutral', ...
    'Unit', 'subject', ...
    'Alpha', 0.05);

plot_hrf_time_unfolding_stats(stats_nps);
```

Map-score-only studies rebuilt from CSVs do not contain the original event
table or time series, so `PlotType='trialmean'` is unavailable. To inspect
run-level event-locked trial means, rebuild from the saved result MAT files:

```matlab
study_runs = hrf_input_table_to_study(input_table, ...
    'LoadResultMat', true, ...
    'LoadWholeBrain', false, ...
    'IncludeMapScores', false);

plot_hrf_study_by_subject(study_runs, ...
    'PlotType', 'trialmean', ...
    'Condition', 'nback-stimblock', ...
    'Unit', 'run', ...
    'TrialOutlierPolicy', 'huber');
```

`hrf_score_wholebrain_input_table` validates existing score CSVs against the
matching metadata sidecar. If row counts, condition labels, or lag indices do
not match, stale score CSVs are regenerated by default (`'OverwriteStale',
true`). If the metadata row count does not match the number of volumes in the
4D beta/T NIfTI, regenerate the whole-brain maps and metadata together before
scoring; the score table cannot safely infer which volumes belong to omitted
conditions such as `*_ttl_*`.

For map-score studies, `Model='sfir'`, `Model='canonical'`, or
`Model='spline'` selects the source whole-brain model recorded by
`hrf_collect_wholebrain_outputs`. Prefer filtering during study construction
with `'SourceModel','sfir'`/`'canonical'`/`'spline'`. For compatibility,
`'ModelName','canonical'` is also interpreted as source-model selection when
the input table has a whole-brain `model` column; it does not relabel FIR rows
as canonical. You can still use `Model='mapscore'` with
`'SourceModel','sfir'` if you prefer to be explicit. Column names such as
`sig_all_NPS` mean "signature set `all`, signature `NPS`." Rows in the score
CSV identify the HRF condition and lag. These values are therefore NPS
expression of HRF beta/T maps, not an NPS time series computed directly from
the original 4D BOLD volumes.

When `Signature` is omitted, the selected curve is `mean_mapscore`, an average
over the loaded score columns. Pass a specific signature or map column, such
as `'Signature','sig_all_NPS'`, when you want that one score instead of the
average.

Map-score curves can look more jagged than `run_hrf_pipeline` 1D HRF fits
because each point is a spatial pattern score from a separate condition-lag
beta/T image. When beta scores are written with a matching `_se.nii` image,
`hrf_apply_maps_to_wholebrain` adds propagated score-SE columns such as
`sig_all_NPS_se`, computed as `sqrt(sum((signature_weight .* beta_se).^2))`
under a diagonal voxel-covariance approximation. `hrf_input_table_to_study(...,
'Object','beta')` uses those propagated SE columns for `mapscore.se`. For
older score CSVs without propagated SE columns, matching `t_scores_file`
values are still used as a fallback approximation, `abs(beta_score / t_score)`;
set `'ApproxSEFromT', false` to disable that fallback. Use
`hrf_time_unfolding_stats` for across-subject mean/SEM and p-values, or refit
the source 1D mean/signature/ROI time series with `run_hrf_pipeline` when you
need full model covariance.

Inspect noisy runs/subjects before interpreting large ribbons or strange
averages:

```matlab
qc = hrf_qc_study_curves(study_scores, ...
    'Model', 'sfir', ...
    'Signature', 'sig_all_NPS', ...
    'Condition', 'nback-stimblock', ...
    'Unit', 'run', ...
    'Weighting', 'huber', ...
    'Plot', true);

qc.table(qc.table.is_outlier, :)

plot_hrf_study_by_subject(study_scores, ...
    'Model', 'sfir', ...
    'Signature', 'sig_all_NPS', ...
    'Condition', 'nback-stimblock', ...
    'Unit', 'run', ...
    'CurveWeights', qc.table.weight);
```

For true trial-level inspection/down-weighting, use the original run-level
`results` structs rather than map-score-only studies:

```matlab
plot_hrf_by_condition(results, ...
    'PlotType', 'trialmean', ...
    'Condition', 'nback-stimblock', ...
    'TrialOutlierPolicy', 'huber');
```

For event-related designs with multiple event instances in a run, use
`PlotType='trialmean'` to plot raw event-locked means with SEM across repeated
events at each TR. That SEM is a different quantity from map-score or FIR beta
SE: it summarizes trial-to-trial variability in the extracted 1D signal, not
uncertainty in the condition-lag beta map.

## Multilevel time-unfolding significance testing

```matlab
study = run_hrf_study_pipeline(fmri_files, events_files, subject_ids, ...
    'TR', 0.8, 'SignalSource', 'mean');

% Subject-level mean extracted timeseries (subjects x time)
study.mean_timeseries

% Within-subject condition contrast, then group-level test across time bins
stats = hrf_time_unfolding_stats(study, ...
    'Model', 'sfir', ...
    'Unit', 'subject', ...
    'ConditionA', 'pain', ...
    'ConditionB', 'neutral', ...
    'Alpha', 0.05);

plot_hrf_time_unfolding_stats(stats);
```

Duplicate subject IDs are averaged before testing by default, so runs are
summarized at the subject level. Use `'Unit', 'run'` to test run-level rows
instead. Empty or failed results are skipped with warnings by default; use
`'MissingPolicy', 'error'` to stop on the first missing model/condition/result.

Signature-specific study-level testing and plotting use the same signature
names stored in `results.fits_by_signature`:

```matlab
stats_nps = hrf_time_unfolding_stats(study, ...
    'Model', 'sfir', ...
    'Signature', 'NPS', ...
    'ConditionA', 'pain', ...
    'ConditionB', 'neutral', ...
    'Unit', 'subject');

plot_hrf_time_unfolding_stats(stats_nps);

plot_hrf_study_by_subject(study, ...
    'Model', 'sfir', ...
    'Signature', 'NPS', ...
    'Condition', 'pain', ...
    'Unit', 'subject');
```

Optional between-group testing is available by supplying `Group` labels
(two-group t-test per time bin).

## Efficient ROI- and pattern-based HRF fitting

Yes - this is often more efficient than voxelwise fitting. You can now do:

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
