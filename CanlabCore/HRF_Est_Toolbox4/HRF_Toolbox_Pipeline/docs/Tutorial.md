# HRF Toolbox Pipeline — Tutorial

A walkthrough for going from raw fMRI data to thresholded group HRF maps using the `HRF_Toolbox_Pipeline`. Aimed at users who have done some neuroimaging in MATLAB but are new to this toolbox.

For the architectural overview, see [Architecture.md](Architecture.md). For the canonical Lindquist-lab fitters this pipeline wraps, see the parent folder `HRF_Est_Toolbox4/`.

---

## Contents

1. [Setup](#1-setup)
2. [Quick start: one subject, one run](#2-quick-start-one-subject-one-run)
3. [Inspecting the results struct](#3-inspecting-the-results-struct)
4. [Working with the new OOP classes](#4-working-with-the-new-oop-classes)
5. [Study scale: writing a SLURM array job](#5-study-scale-writing-a-slurm-array-job)
6. [Auditing and repairing a finished SLURM run](#6-auditing-and-repairing-a-finished-slurm-run)
7. [Post-SLURM: collecting outputs and building a study struct](#7-post-slurm-collecting-outputs-and-building-a-study-struct)
8. [Group-level analysis](#8-group-level-analysis)
9. [Visualization: plots, montages, animations](#9-visualization-plots-montages-animations)
10. [Common patterns and FAQ](#10-common-patterns-and-faq)

---

## 1. Setup

You need three folders on the MATLAB path: CanlabCore (the parent class library), this pipeline subfolder, and the legacy Lindquist fitters at `HRF_Est_Toolbox4/`. Adding CanlabCore recursively with `genpath` covers all three:

```matlab
addpath(genpath('/path/to/CanlabCore/CanlabCore'));
```

That's it. The `genpath` walks into `HRF_Est_Toolbox4/` (legacy fitters), `HRF_Est_Toolbox4/HRF_Toolbox_Pipeline/` (this pipeline), and all the class folders (`@fmri_data`, `@statistic_image`, `@fmri_hrf`, `@statistic_hrf`, etc.).

Sanity check:

```matlab
which run_hrf_pipeline                 % should resolve into HRF_Toolbox_Pipeline/
which fmri_data                        % should resolve into @fmri_data/
which fmri_hrf                         % should resolve into HRF_Toolbox_Pipeline/@fmri_hrf/
which Fit_sFIR                         % should resolve into HRF_Est_Toolbox4/ (sibling)
```

If any of these come up empty, your path is missing something.

---

## 2. Quick start: one subject, one run

The shortest path from raw data to fitted HRFs. Single subject, single run, three models:

```matlab
results = run_hrf_pipeline( ...
    '/path/to/sub-01_task-pain_bold.nii.gz', ...     % 4D fMRI
    '/path/to/sub-01_task-pain_events.tsv', ...      % BIDS events
    'TR', 0.8, ...
    'MaskNii', '/path/to/brain_mask.nii.gz', ...
    'Conditions', {'pain', 'neutral'}, ...
    'WindowSeconds', 30, ...
    'Models', {'fir', 'sfir', 'canonical'}, ...
    'OutputMat', '/path/to/sub-01_hrf_results.mat');
```

What this does, step by step:

1. Loads the 4D NIfTI and extracts a mean-within-mask timeseries.
2. Reads the BIDS `events.tsv`; validates required columns (`onset`, `duration`, `trial_type`).
3. Builds one stick function per condition.
4. Fits each requested model (here: FIR, sFIR, SPM canonical) using the Lindquist-lab fitters.
5. Returns a single `results` struct and (optionally) saves it as `.mat`.

A few practical defaults:

- If you don't pass `TR`, the pipeline reads it from the NIfTI header (`PixelDimensions(4)`). Pass it explicitly when you can — header values are sometimes wrong.
- `WindowSeconds` is the HRF estimation window (the time after stimulus onset over which the HRF is estimated). 30s covers a canonical BOLD response with margin. For event-related rapid designs you can shrink to 24s.
- `Conditions` filters which `trial_type` values from the events file you actually fit. Wildcards work: `{'pain_*'}` matches `pain_high`, `pain_low`, etc.

---

## 3. Inspecting the results struct

The returned `results` struct has the same shape regardless of how many models you fit. Top-level fields:

```matlab
results.config              % the parameters the pipeline saw
results.timeseries          % whole-brain (or in-mask) mean timeseries it actually fit
results.stick_functions     % one column per condition
results.fits_by_model       % per-model fit struct (the heart of the output)
results.wholebrain_by_model % only if you ran whole-brain fitting (see Section 5)
```

The per-model fit struct (`results.fits_by_model.sfir`, say) carries:

```matlab
results.fits_by_model.sfir.hrf            % time x conditions HRF estimates
results.fits_by_model.sfir.hrf_se         % time x conditions standard errors
results.fits_by_model.sfir.fit            % fitted timeseries (overall)
results.fits_by_model.sfir.residuals      % residuals
results.fits_by_model.sfir.amplitude      % per-condition amplitude estimate
results.fits_by_model.sfir.peak_lag       % per-condition peak lag
results.fits_by_model.sfir.metadata_table % one row per (condition, lag) pair
```

For a quick visual sanity check:

```matlab
plot_hrf_by_condition(results, 'Models', {'sfir', 'canonical'});
```

That gives you per-condition curves with error bands, one panel per condition, one curve per model.

---

## 4. Working with the new OOP classes

The `@fmri_hrf` and `@statistic_hrf` classes (added in Phase 3 of the toolbox refactor) wrap whole-brain HRF maps in objects that inherit from CanlabCore's `fmri_data` and `statistic_image`. This means all the canlabCore methods you already use — `threshold`, `multi_threshold`, `orthviews`, `apply_mask`, `extract_roi_averages`, `convert2mask`, `conjunction` — work on these objects for free.

### Building a paired (beta, t) pair from one fit

The canonical entry point is `make_fmri_stat_hrf`. It takes one of three source shapes and returns the matched pair:

```matlab
% From a results struct (most common — right after run_hrf_pipeline with whole-brain on):
[Hb, Ht] = make_fmri_stat_hrf(results.wholebrain_by_model.sfir, ...
    'Subject', 'sub-01', ...
    'RunLabel', 'task-pain_run-01', ...
    'ModelName', 'sfir', ...
    'TR', 0.8);

% From a saved result.mat path:
[Hb, Ht] = make_fmri_stat_hrf('/path/to/sub-01_hrf_results.mat', ...
    'ModelName', 'sfir', ...
    'Subject', 'sub-01', ...
    'TR', 0.8);

% From a NIfTI prefix (loads beta + t + se + metadata from disk):
[Hb, Ht] = make_fmri_stat_hrf('/path/to/sub-01_hrf_sfir', ...
    'Subject', 'sub-01', ...
    'ModelName', 'sfir', ...
    'TR', 0.8);
```

`Hb` is the beta side (`fmri_hrf < fmri_data`), `Ht` is the inferential side (`statistic_hrf < statistic_image`). They share metadata (subject, run_label, model_name, conditions, lags, TR) by construction.

### Inspecting an object

The classes have HRF-aware `disp` methods:

```matlab
>> Hb
  fmri_hrf  subject=sub-01  run=task-pain_run-01  model=sfir
            99837 voxels x 60 volumes  (2 conditions x 30 lags)  TR=0.8s
            conditions: pain, neutral

>> Ht
  statistic_hrf  subject=sub-01  run=task-pain_run-01  model=sfir  type=T
                 99837 voxels x 60 volumes  (2 conditions x 30 lags)  TR=0.8s
                 conditions: pain, neutral
```

### Thresholding (inherited from statistic_image)

No new code needed — `threshold`, `multi_threshold`, `orthviews` all work:

```matlab
Ht_thr = threshold(Ht, 0.001, 'unc');
orthviews(Ht_thr);
```

### Concatenating across subjects/runs

`cat` is overridden to align HRF metadata across (subject, run) and stack:

```matlab
% After running the pipeline for every subject in your study:
Hb_study = cat(1, Hb_sub01, Hb_sub02, Hb_sub03, ..., Hb_subN);
Ht_study = cat(1, Ht_sub01, Ht_sub02, Ht_sub03, ..., Ht_subN);
```

`cat` requires all inputs to share `model_name`, `TR`, and the condition set. Subject and run_label can (and should) differ — the resulting `metadata_table` gets `subject` and `run_label` columns added so you can recover the source axis with table indexing later.

### Pairing an existing fmri_hrf with an SE image

If you already have an `fmri_hrf` with betas and want to derive its statistic_hrf sibling on the fly:

```matlab
Ht = Hb.to_statistic_hrf('SE', se_obj);
% or, with an already-computed t image:
Ht = Hb.to_statistic_hrf('TStat', t_obj);
```

### What's not yet on these objects

v0 ships the constructor, `disp`, `cat`, and `to_statistic_hrf`. Methods deferred to Phase 4 (the misspecification pipeline):

- `plot(Hb)` — condition × lag curves
- `montage(Ht)` — animated over lag
- `residuals(Hb)`
- `misspec_metrics(Hb, 'Reference', ...)`
- `cv_amplitude(Hb)`
- `ttest(Hb)` — returns a `statistic_hrf`
- `group_hrf(Hb)` — HMHRFest-style data-driven group HRF

For now, the older procedural functions (`plot_hrf_by_condition`, `plot_hrf_study_by_subject`, `hrf_make_average_montage_animations`) cover the same ground.

---

## 5. Study scale: writing a SLURM array job

For a real study (dozens to hundreds of subjects), you don't run `run_hrf_pipeline` in a loop locally — you generate a SLURM array script that fans out one task per (subject, run, model) onto a cluster.

```matlab
% Build the input lists. One entry per (subject, run).
fmri_files = { ...
    '/dartfs/.../sub-01_task-pain_bold.nii.gz', ...
    '/dartfs/.../sub-02_task-pain_bold.nii.gz', ...
    ... };
events_files = { ...
    '/dartfs/.../sub-01_task-pain_events.tsv', ...
    '/dartfs/.../sub-02_task-pain_events.tsv', ...
    ... };
subject_ids = {'sub-01','sub-02', ...};

hrf_write_slurm_study_script( ...
    fmri_files, events_files, subject_ids, ...
    'JobName', 'hrf_pain_study', ...
    'OutputDir', '/dartfs/.../hrf_outputs', ...
    'TR', 0.8, ...
    'Conditions', {'pain', 'neutral'}, ...
    'WindowSeconds', 30, ...
    'Models', {'sfir', 'canonical'}, ...
    'WriteWholeBrain', true, ...
    'SignatureSets', {'all'}, ...
    'PartitionTime', '04:00:00');
```

The generator writes three things into `OutputDir`:

- `hrf_pain_study.sh` — the SLURM array script you `sbatch`.
- `hrf_pain_study_worker.m` — the per-task MATLAB worker. Each array task runs this with its task index.
- `hrf_pain_study_manifest.csv` — one row per array task. The audit step uses this.

Each task does the work of `run_hrf_pipeline` plus, if `'SignatureSets'` is non-empty, scoring the resulting whole-brain maps against the signature set(s) and writing `<prefix>_<object>_map_scores.csv`. The scoring step is wrapped in `try/catch` so a single failing signature doesn't kill the whole task — the audit step (next section) catches anything genuinely missing.

To submit:

```bash
sbatch /dartfs/.../hrf_outputs/hrf_pain_study.sh
```

---

## 6. Auditing and repairing a finished SLURM run

After the array job finishes, the audit step reconciles the manifest against on-disk outputs and reports any gaps. Two modes:

```matlab
% Report only — fast, read-only.
[audit, summary] = hrf_audit_slurm_outputs('/dartfs/.../hrf_outputs');
summary
audit(~audit.complete, :)

% Report AND repair missing score CSVs in place.
[audit, summary] = hrf_audit_slurm_outputs('/dartfs/.../hrf_outputs', ...
    'RepairMissing', true, ...
    'Verbose', true);
summary.repair
```

`audit` is a table with one row per (task, model). Key columns: `core_complete` (all NIfTIs + metadata present?), `score_complete` (all signature score CSVs present?), `complete` (both). `RepairMissing=true` finds every row where `core_complete && ~score_complete` and re-runs scoring for that prefix — it calls `hrf_score_one_prefix` (the same helper the SLURM worker uses), so the post-repair files are byte-identical to what the worker would have written.

If `score_complete` is still false after repair, look at `summary.repair.errors` — the helper logged why (typically: missing signature in `load_image_set`, voxel-space mismatch, metadata file unreadable).

---

## 7. Post-SLURM: collecting outputs and building a study struct

Once the array job is clean, collect the outputs into a single table for downstream work:

```matlab
input_table = hrf_collect_wholebrain_outputs('/dartfs/.../hrf_outputs');
```

`input_table` has one row per (subject, run, model) with paths to beta/t/se/p NIfTIs, metadata CSV, score CSVs, and the run's `result.mat`. This is the "second-level input table" — most downstream analytics consume it.

To rebuild an in-memory study struct (the same shape `run_hrf_study_pipeline` would produce, but assembled after the fact from disk):

```matlab
study = hrf_input_table_to_study(input_table, 'LoadWholeBrain', true);
```

`study.results{i}` now has `.wholebrain` (`fmri_data` / `statistic_image` objects loaded from the NIfTIs) and `.fits_by_signature` (signature scores from the CSVs). Memory cost is non-trivial — `LoadWholeBrain=false` skips the NIfTI loads if you only need the score CSVs.

For the extended version with per-condition contrasts derived from metadata, use:

```matlab
study = hrf_second_level_inputs_to_study(input_table, ...
    'Conditions', {'pain','neutral'});
```

---

## 8. Group-level analysis

Three entry points, depending on what you're testing:

### Time-point-wise t-tests across subjects

```matlab
T = hrf_time_unfolding_stats(study, ...
    'Conditions', {'pain', 'neutral'}, ...
    'Models', {'sfir'}, ...
    'Signatures', {'NPS', 'SIIPS'});

T.timepoint_tstats        % per-timepoint t-tests per (condition, signature)
T.fdr_corrected           % FDR-thresholded version
```

### 2×2 factorial (two between-subject groups × two conditions)

```matlab
F = hrf_2x2_study_score_stats(study, group_labels, ...
    'ConditionPair', {'pain', 'neutral'}, ...
    'Signature', 'NPS');

F.main_effect_group
F.main_effect_condition
F.interaction
```

### Flexible long-form dataframe for custom modeling

```matlab
df = hrf_analyze_second_level_inputs(study);
% df is a wide table: rows = (subject, condition, lag, signature); columns include
%   score, model, run, group, etc. Hand off to fitlme / fitglme / R.
```

---

## 9. Visualization: plots, montages, animations

### Per-subject HRF curves

```matlab
plot_hrf_study_by_subject(study, ...
    'Conditions', {'pain'}, ...
    'Models', {'sfir'}, ...
    'Signatures', {'NPS'});
```

One panel per subject; one curve per condition × model × signature. The biggest plotter in the toolbox; expects a study struct from Section 7.

### Group-level animated montages (the big one)

The standout deliverable for whole-brain visualization. Averages whole-brain maps across runs within subjects, then across subjects, per condition, and writes an MP4 that animates the group t-map across HRF lags:

```matlab
avg = hrf_make_average_montage_animations(input_table, ...
    'Model', 'sfir', ...
    'Object', 'beta', ...
    'Conditions', {'pain', 'neutral'}, ...
    'OutputDir', '/dartfs/.../animations', ...
    'GroupStatistic', 't', ...
    'GroupCorrection', 'fdr', ...
    'GroupAlpha', 0.05, ...
    'FrameStep', 1, ...
    'FPS', 8);
```

This is I/O heavy — it reads the per-subject whole-brain NIfTIs from disk. The pipeline caches per-prefix loads so each NIfTI is read at most once across all conditions; without that cache, a 4-condition run reads each subject's maps 4 times. The `'Verbose'` flag (on by default) prints per-subject timing and a `[cache hit]` indicator from condition 2 onward — useful sanity check.

If the run is unexpectedly slow, the cache tag tells you whether the bottleneck is the network (`[load N/N]` per subject in conditions 2+) or the rendering itself (`[cache hit]` per subject but condition wall-clock is huge).

### Other visualizations

- `plot_hrf_by_condition` — quick per-condition curves from a single `results` struct.
- `plot_hrf_2x2_study_score_stats` — visualizes a 2×2 factorial result.
- `plot_hrf_time_unfolding_stats` — group-level timepoint-wise t-tests with FDR overlay.
- `hrf_animate_wholebrain_stats` — animate a single subject's whole-brain HRF maps over lag.

---

## 10. Common patterns and FAQ

### "I have a result.mat from an old run. How do I get to the new objects?"

```matlab
[Hb, Ht] = make_fmri_stat_hrf('/path/to/old_results.mat', ...
    'ModelName', 'sfir', ...
    'Subject', 'sub-01', ...
    'TR', 0.8);
```

### "How do I find an event I lost?"

The audit step is the canonical "what's there, what's missing" tool. Read-only mode is fast:

```matlab
[audit, summary] = hrf_audit_slurm_outputs(output_dir);
summary
audit.failed_reason(~audit.complete)
```

### "My signature scoring is slow even with the cache."

The cache eliminates redundant *loads*. If the signatures themselves are expensive (image-similarity computation, large signature sets), per-signature time still adds up. Run on a smaller signature subset to confirm:

```matlab
'SignatureSets', {'NPS', 'SIIPS'}     % instead of 'all'
```

### "I need to add a new HRF model."

The legacy fitters live at `HRF_Est_Toolbox4/` (parent folder, not this pipeline). Add your `Fit_<name>.m` there. Then add a dispatch case in `hrf_fit_all_models` so `'Models', {'<name>'}` routes to it. The pipeline entry points pick it up automatically.

### "What's an `output_prefix` vs a `model_prefix`?"

For a single-model fit, they're the same. For a multi-model fit (you passed `'Models', {'sfir', 'canonical'}` to the pipeline), each model writes to its own prefix:

- `output_prefix` = task-level prefix, e.g. `/path/sub-01_hrf`
- `model_prefix` = `<output_prefix>_<model>`, e.g. `/path/sub-01_hrf_sfir`

The audit table carries both columns so you can index by either.

### "Where do the threshold methods live for `statistic_hrf`?"

Inherited from `@statistic_image/`. `threshold`, `multi_threshold`, `orthviews`, `convert2mask`, `conjunction` all work unchanged. The Phase 3 OOP refactor was specifically built so you don't have to reimplement them.

### "How do I cite this?"

The legacy fitters (`Fit_sFIR`, `Fit_Logit2`, `HMHRFest`, `ResidScan`, `PowerLoss`) are Lindquist-lab and should cite their original papers. The pipeline glue and OOP refactor are part of Michael Sun's HRF Toolbox work (Sun et al., in preparation).

---

For the architectural breakdown of how all these pieces fit together, see [Architecture.md](Architecture.md). For the four-phase refactor plan, see [HRF_Pipeline_Phased_Plan.md](HRF_Pipeline_Phased_Plan.md).
