# HRF Toolbox: Phased Refactor & Extension Plan

**Author:** Michael Sun (drafted with Claude, May 2026)
**Audience:** Tor (review), Michael (execution)
**Scope:** Four interlocking workstreams across `HRF_Est_Toolbox2`, `HRF_Est_Toolbox4`, and a new `HRF_Toolbox_Pipeline` subfolder.

---

## At a glance

Four workstreams, sequenced so each phase unblocks the next:

| Phase | Workstream | Deliverable | Effort | Blocks |
|---|---|---|---|---|
| 1 | Audit/score gap fix | Patched `hrf_audit_slurm_outputs` + shared scoring helper | Small (1–2 days) | nothing — ship now |
| 2 | Repo housekeeping + Toolbox2/4 merge | Cleaned `MichaelSun` branch; new `HRF_Toolbox_Pipeline/` in Toolbox4 | Medium (3–5 days) | requires parent‑directory access |
| 3 | `fmri_hrf` + `statistic_hrf` OOP refactor | New `@fmri_hrf` (< fmri_data) and `@statistic_hrf` (< statistic_image) classes with `plot`, `montage`, `cat`, residual/misspec methods | Medium‑large (1–2 weeks for v0) | Phase 2 |
| 4 | Misspecification pipeline | Diagnostics + CV + group t‑tests, packaged as `fmri_hrf` methods | Medium (1–2 weeks) | Phase 3 |

The phases are sequenced so that the OOP refactor (Phase 3) becomes the substrate for the misspecification pipeline (Phase 4) — Tor's "maximally reuse canlabCore" goal is much easier to honour when the new diagnostics are written as `fmri_hrf` methods than as a parallel set of free functions.

---

## Phase 1 — Audit/Score Gap Fix (immediate)

### Diagnosis

The pipeline already attempts per‑task scoring. Looking at the generated SLURM worker in `hrf_write_slurm_study_script.m` (lines 282–301), each SLURM array task calls `hrf_apply_maps_to_wholebrain` after writing the beta/T/SE/P NIfTIs, conditional on `config.signature_sets` or `config.image_sets` being non‑empty. So in principle the score CSVs should already exist by the time `hrf_audit_slurm_outputs` runs.

Three concrete reasons they're often missing:

1. **Scoring is not wrapped in `try/catch` in the SLURM worker.** Per worker lines 282–301, the inline scoring loop has no error handling. If `load_image_set` fails for a single signature, if a signature's voxel space cannot be resampled, or if metadata/volume counts disagree, the entire SLURM task fails *after* `beta/t/se/p/metadata` are written. The audit then sees those NIfTIs but no score CSVs.
2. **`hrf_audit_slurm_outputs` only reports, it never repairs.** `local_score_status` (audit.m:380–396) checks file existence and surfaces the gap in `failed_reason`, but nothing in the audit takes action.
3. **The post‑hoc backfill (`hrf_score_wholebrain_input_table`) is slow because it re‑reads whole‑brain NIfTIs from disk per row.** Lines 467 and 505 of that file each call `statistic_image(fmri_data(image_file, 'noverbose'))` separately for beta and SE — three full whole‑brain reads per row (beta, t, se), iterated serially with no parallelization.

The duplication is also a maintenance hazard: there are now **three** places that orchestrate scoring (the SLURM worker, the audit's existence check, the post‑hoc backfill), and only the SLURM worker actually calls `hrf_apply_maps_to_wholebrain`.

### Proposed fix (smallest viable change)

Three coordinated edits:

**1.1 — Factor a single per‑prefix scoring helper.**
Create `hrf_score_one_prefix.m` that, given an `output_prefix`, `model_name`, `n_models`, `score_objects`, `signature_sets`, `image_sets`, and metadata‑resolution policy, writes the missing `*_<object>_map_scores.csv` for one (task, model). This is the only function that knows how to: (a) resolve metadata via the current "csv → result.mat → sibling csv" fallback chain in `hrf_score_wholebrain_input_table`; (b) load the beta/t/se NIfTIs; (c) call `hrf_apply_maps_to_wholebrain`. It returns a status struct.

**1.2 — Replace inline scoring in the SLURM worker with a call to the helper, wrapped in `try/catch`.**
The worker emits a per‑task scoring loop today. Replace it with `hrf_score_one_prefix(output_prefix, model_name, …)` wrapped in `try`/`catch`, log the error, and continue. A signature failure should no longer kill a whole SLURM task; the audit will catch any genuinely missing CSVs.

**1.3 — Add `RepairMissing` mode to `hrf_audit_slurm_outputs`.**
New name‑value parameter `'RepairMissing', false` (default `false`, preserving current behaviour). When `true`, for every row where `core_complete && ~score_complete`, call `hrf_score_one_prefix` to backfill the missing CSVs in place. Optionally accept `'UseParallel', true` to run the repair loop under `parfor`. After repair, re‑evaluate `score_complete` and write the audit table with the repaired status.

`hrf_score_wholebrain_input_table` becomes a thin wrapper over `hrf_score_one_prefix` — keep it for back‑compat but stop carrying the heavy logic.

### Why this design rather than "just fold scoring into the audit"

Folding scoring straight into the audit (the literal reading of "the audit should produce the score files") would solve the symptom but lock the duplication in. The factored helper costs the same effort and lets the **fastest** path (per‑task scoring inside SLURM) and the **slow** path (post‑hoc repair) share the same code, so any future optimization — e.g., memory‑mapped NIfTI reads, GPU‑accelerated dot products, batched signature application — benefits all three call sites at once.

### Verification (Phase 1 acceptance)

Smoke test on one finished study directory:

```matlab
% Before fix
[A1, S1] = hrf_audit_slurm_outputs(output_dir);
sum(~A1.score_complete & A1.core_complete)   % expected: some > 0

% After fix
[A2, S2] = hrf_audit_slurm_outputs(output_dir, 'RepairMissing', true);
sum(~A2.score_complete & A2.core_complete)   % expected: 0
```

Then diff one repaired `*_map_scores.csv` against an `hrf_score_wholebrain_input_table` re‑run on the same row — they should be identical to floating‑point tolerance.

---

## Phase 2 — Housekeeping & Toolbox2 ↔ Toolbox4 Merge

### Cruft inventory (from `HRF_Est_Toolbox2/`)

A first pass over the working folder surfaces these candidates for cleanup. I have not yet diffed against `HRF_Est_Toolbox4` (no access to it from my current mount), so the merge decisions in the right column are provisional.

| Item | Likely status | Provisional disposition |
|---|---|---|
| `hrf_write_slurm_study_script.asv` | MATLAB autosave file | Delete |
| `Old_stuff/` (15 files, 140 KB) | Pre‑refactor logit/sFIR variants | Archive in git history, then delete from working tree |
| `Old_stuff/More_recent_old_stuff/` | Deeper legacy, including `.mat` fixtures | Archive in git history, then delete from working tree |
| `New_fminsearchbnd/` | Single‑purpose third‑party utility | Keep but move to a shared deps folder; check whether canlabCore already vendors it |
| `Example.m`, `Example2.m` | Legacy demos that reference old fitters | Either modernize as `examples/` or remove if superseded by README's quick‑start |
| `HRF Analysis Code Map.pptx`, `README.pptx` | Documentation artifacts | Move to `docs/` |
| Top‑level non‑`hrf_` `.m` files (Fit_Canonical_HRF, Fit_Logit2, Fit_sFIR, Fit_NLgamma, Fit_Spline, etc.) | The legacy fitters that the new `run_hrf_pipeline` wraps; load‑bearing | **Keep** — these are the canonical Lindquist‑lab estimators and the new pipeline depends on them |

I'd like to walk through this list with you before deleting anything — particularly the `Old_stuff/More_recent_old_stuff/*.mat` files, which look like reference fixtures rather than code.

### Toolbox2 ↔ Toolbox4 merge — what I need from you

To do the merge credibly I need access to **the parent `CanlabCore` directory** (or at least `HRF_Est_Toolbox4`). Your selected folder right now is just `HRF_Est_Toolbox2`, so I cannot:

- diff Toolbox2 vs Toolbox4 method‑by‑method,
- check which functions are newer / diverged / unique to each side,
- inventory the canlabCore `@fmri_data` / `@statistic_image` class folders that Phase 3 must inherit from.

Two options, either works:

- **(a)** Re‑select the workspace folder one level up at `MichaelSun/CanlabCore/CanlabCore/` so I can see Toolbox2, Toolbox4, and the canlabCore class folders side by side.
- **(b)** Keep me here and run a `git log --stat HRF_Est_Toolbox2 HRF_Est_Toolbox4` for me, plus paste me the file listing of `HRF_Est_Toolbox4/`. I can then propose the merge map.

(a) is much faster.

### Merge approach (once I can see both)

1. **Build a function‑by‑function comparison table.** For every `*.m` in either toolbox, classify as `only_in_2`, `only_in_4`, `diverged`, `identical`. For diverged files, prefer the newer mtime *but* manually inspect anything where Toolbox4 supersedes a Lindquist‑lab core fitter — those are load‑bearing and changes should be reviewed.
2. **Scaffold `HRF_Est_Toolbox4/HRF_Toolbox_Pipeline/`** as the new home for the user‑facing pipeline (the `run_hrf_pipeline`, `hrf_write_slurm_study_script`, audit, collect, score, plot, time‑unfolding family). The legacy fitters (`Fit_*`, `Anneal_*`, `ResidScan`, `PowerLoss`, etc.) stay at the Toolbox4 root next to the canonical Lindquist code.
3. **Move the new pipeline files** from Toolbox2 to `HRF_Toolbox_Pipeline/`, preserving git history with `git mv`.
4. **Delete Toolbox2** (or leave an empty `HRF_Est_Toolbox2/README.md` redirecting to Toolbox4) after one full pass of regression tests on the pipeline.

### Verification (Phase 2 acceptance)

- `run_hrf_pipeline` runs end‑to‑end on the smoke‑test subject from the new location.
- `hrf_write_slurm_study_script` generates a script and worker that point at the new paths.
- A `git log --follow` on any moved file shows the pre‑merge history.

---

## Phase 3 — The `fmri_hrf` and `statistic_hrf` Objects (OOP refactor)

### The two-class refinement (Michael's design)

Rather than one polymorphic `hrf` object, build **two sibling classes** that mirror the canlabCore split between data containers and inferential containers:

- **`fmri_hrf < fmri_data`** — holds beta / amplitude maps. These are continuous, unthresholded numbers; the natural operations are amplitude/AUC/peak-lag extraction, ROI averaging, signature application, time-series reconstruction.
- **`statistic_hrf < statistic_image`** — holds t / p / se maps. Inherits `threshold`, `multi_threshold`, `orthviews`, `convert2mask`, `conjunction`, etc., so thresholding "on the fly" comes free.

Both classes share the same HRF-specific metadata (subject, run, condition, lag, TR, design info) and the same `cat` semantics across runs/subjects. They're paired (one `fmri_hrf` and one `statistic_hrf` per fit) and convertible into each other (e.g. `statistic_hrf(fmri_hrf_obj, 'SE', se_obj)` to assemble inferential maps from a beta object + uncertainty).

This is significantly cleaner than a single class with optional inferential fields: it puts `threshold()` in exactly the place the rest of CanlabCore puts it, and avoids the awkward situation of an `hrf.threshold()` method that does nothing for raw amplitudes.

### Motivating example

```matlab
% Build from one run
[Hb, Ht] = make_fmri_stat_hrf(results.wholebrain_by_model.sfir, ...
    'Subject', 'sub-01', ...
    'Run', 'task-pain_run-01', ...
    'Conditions', {'pain','neutral'});
% Hb is fmri_hrf, Ht is statistic_hrf

% Concatenate across runs / subjects (both classes work identically)
Hb_sub  = cat(1, Hb_run01, Hb_run02, Hb_run03);
Hb_study = cat(1, Hb_sub01, Hb_sub02, ..., Hb_subN);
Ht_study = cat(1, Ht_sub01, Ht_sub02, ..., Ht_subN);

% Threshold the t-side on the fly (inherited from statistic_image)
Ht_thr = threshold(Ht_study, 0.001, 'unc');

% Combined plot: brain montage from one side + time-series line plot from the other
plot(Hb_study, 'Condition', 'pain');                                  % amplitude curves
plot(Ht_study, 'Condition', 'pain', 'Threshold', [0.05 'unc']);       % thresholded t montage
montage(Ht_study, 'Condition', 'pain');                               % animated over lag

% Misspecification & CV operate on the amplitude side
residuals(Hb)                                                          % single-run residual diagnostics
misspec_metrics(Hb, 'Reference', 'canonical')
cv_amplitude(Hb_study, 'KFold', 5);
ttest(Hb_study, 'ConditionA','pain','ConditionB','neutral')           % returns a statistic_hrf
```

### Class skeletons (proposed)

```
@fmri_hrf/
  fmri_hrf.m         % constructor; accepts wholebrain.b, an fmri_data + metadata,
                     %               a result.mat path, or a NIfTI prefix.
  cat.m              % vertical concatenation across (subject, run) axes
  plot.m             % time-series-style condition x lag plot, optional montage
  residuals.m        % requires design_matrix; returns residual time-series
  misspec_metrics.m  % residual-task-corr, AUC ratio, peak-lag bias, power loss
  cv_amplitude.m     % k-fold CV of amplitude / AUC / sqrt-variance
  ttest.m            % returns a statistic_hrf
  group_hrf.m        % HMHRFest-style data-driven group HRF
  display.m, isempty.m, size.m, end.m
  to_statistic_hrf.m % helper: build a paired statistic_hrf from a se obj
  private/
    init_from_wholebrain.m
    init_from_resultmat.m
    align_for_concat.m
    hrf_metadata_validate.m  % shared with @statistic_hrf
```

```
@statistic_hrf/
  statistic_hrf.m    % constructor; accepts wholebrain.t, statistic_image + metadata,
                     %               or a NIfTI prefix
  cat.m              % parallel to @fmri_hrf/cat
  plot.m             % thresholded montage + condition x lag t-curves
  montage.m          % animated over lag (calls @fmridisplay)
  threshold.m        % override-or-delegate to @statistic_image/threshold
                     % (likely just inherit and let the parent handle it)
  display.m, isempty.m, size.m, end.m
```

### Shared HRF metadata (lives on both classes)

```matlab
% Same fields on fmri_hrf and statistic_hrf
properties
    metadata_table   % volume_index, condition, lag_index, lag_seconds,
                     %   image_label, N, dfe, TR, mode (one row per 4D volume)
    subject          % string
    run_label        % string
    model_name       % 'fir' | 'sfir' | 'canonical' | 'spline'
    conditions       % cellstr
    TR               % seconds
    design_matrix    % NaN if not loaded; required by residuals()
    design_info      % output_lift, fir_columns, intercept/nuisance columns
    source_paths     % beta/t/se/p/metadata/result_mat files
    hrf_meta_version % bump on schema change
end
```

A small `hrf_metadata_validate.m` lives under both class folders' `private/` (or in a shared helper folder) and is the single source of truth for "are these fields well-formed and compatible for concatenation."

### canlabCore reuse — what each class inherits vs adds

| Behaviour | Lives in |
|---|---|
| Voxel storage, masking, write/read, removed_voxels semantics | `fmri_data` (parent of `fmri_hrf`) and `statistic_image` (parent of `statistic_hrf`). |
| `threshold`, `multi_threshold`, `orthviews`, `convert2mask`, `conjunction` | `statistic_image` — free for `statistic_hrf`. |
| `apply_mask`, `extract_roi_averages`, `pattern_expression`, `image_similarity_plot` | `fmri_data` / `image_vector` — free for `fmri_hrf`. |
| HRF-side: condition-lag indexing, residual diagnostics, misspec metrics, group HRF, k-fold CV | New methods on `fmri_hrf` (data-side analytics live there). |
| HRF-side: thresholded montage over lag, condition x lag t-curves | New methods on `statistic_hrf`. |
| `cat` across (subject, run) rows | New on both classes — different from `statistic_image/cat` which concatenates *images*. |
| Pairing | A small helper `[Hb, Ht] = make_fmri_stat_hrf(...)` builds the matched pair from one fit. `Hb.to_statistic_hrf()` and `statistic_hrf(Hb, 'SE', se_obj)` round-trip. |

### Confirmed from the CanlabCore mount

With the second mount in place I can see the parent classes' actual surface area: `@fmri_data/` provides `cat.m`, `hrf_fit.m`, `extract_roi_averages.m`, `denoise_timeseries_pipeline.m`, `fitlme_voxelwise.m`, and the data-side methods we want for `fmri_hrf`. `@statistic_image/` provides `cat.m`, `threshold.m`, `multi_threshold.m`, `orthviews.m`, `convert2mask.m`, `conjunction.m`, `estimateBayesFactor.m` — every method `statistic_hrf` needs for thresholded inference is already there.

The one wrinkle: `@statistic_image/cat.m` already exists, and our `@statistic_hrf/cat.m` must call it (for the underlying voxel data) and then align the HRF metadata. Same on the fmri_data side.

### Open questions for Tor

1. **One pair per (subject, run, model), or one pair per (subject, run) with multiple models?** Multi-model pairing makes misspec comparison (Phase 4) cleaner — `Hb_sfir` vs `Hb_canonical` for the same subject are siblings rather than separately-concatenated objects. My instinct: single-model classes plus an `hrf_modelset` container that holds `{model_name → {fmri_hrf, statistic_hrf}}`. Misspec methods take an `hrf_modelset` and return cross-model comparisons.
2. **Does the 1D signature/ROI time series live on `fmri_hrf` or on a parallel `hrf_timeseries`?** They share metadata (subject, run, conditions, TR) but the storage is completely different (`fmri_hrf` is voxels × condition·lag, the time series is time × signatures).
3. **Naming** — `fmri_hrf` is clean and fits the `fmri_data` / `fmri_mask_image` / `fmri_timeseries` family already in CanlabCore. `statistic_hrf` parallels `statistic_image`. Both are unambiguous in CanlabCore but `hrf` alone would shadow many SPM/legacy variables.

### Verification (Phase 3 acceptance)

- Round‑trip: `H = hrf(prefix); write(H, new_prefix); H2 = hrf(new_prefix); assert(isequal(H, H2))`.
- `plot(H)` renders correctly for a single subject, single run, all 4 supported models.
- `cat(1, H1, H2, ..., HN)` for N=10 subjects renders a study‑level montage in under a minute on a typical workstation.
- All canlabCore inheritance paths exercised by at least one test.

---

## Phase 4 — Misspecification Pipeline

### Goal

Move from "we fit several HRF models and look at them" to "we quantify, for each model, how *wrong* it is, where the wrongness is, and whether a better model is reachable from the data". The deliverable is a set of `fmri_hrf` methods (operating on amplitude maps, with diagnostic outputs that can be lifted into `statistic_hrf` for thresholded visualization) that work at run / subject / study level and emit interpretable metrics.

### Drawing from Lindquist's tools

`HRF_Est_Toolbox2` already vendors several Lindquist‑lab functions that are directly relevant. Once I can see Toolbox4 I'll list its additions, but from Toolbox2 alone the relevant pieces are:

- **`ResidScan.m`** — scans residuals for task‑locked structure. This is the foundation: if residuals at the canonical HRF still show a task‑locked component, the canonical is misspecified for this voxel/region. The misspec pipeline should run this per voxel (or per ROI) and aggregate.
- **`PowerLoss.m`** — quantifies power loss from HRF misspecification. The literature this implements (Loh, Lindquist, Wager 2008; Lindquist et al. 2009) gives the canonical metric: efficiency loss as a function of how far the assumed HRF is from the true one. We can use this both as a metric and as a CV objective.
- **`PowerSim.m`** — Monte‑Carlo helper used by `PowerLoss`. Useful for null distributions of misspec metrics.
- **`Anneal_Logit.m`, `Fit_Logit2.m`, `Det_Logit.m`** — the inverse‑logit HRF family. Provides a flexible HRF basis that's nearly canonical when well‑fit; useful as a "best achievable" reference against which to measure other models' misspec.
- **`HMHRFest.m`** — hierarchical mixed‑effects HRF estimation. The natural backbone for group‑level "optimized HRF agnostic of hypothesis" — this gives you a data‑driven group HRF that you can then test against the canonical.

### Metrics to expose (proposed)

For a fitted `fmri_hrf` object `Hb` with model `m`, compared against reference model `r`:

| Metric | Formula sketch | Interpretation |
|---|---|---|
| `residual_taskcorr(H)` | corr(residual, stick_function) per voxel/ROI | Non‑zero = task‑locked variance unexplained → misspec |
| `residual_acf(H, lag)` | residual autocorrelation, normalized | Distinguishes physiological noise from misspec |
| `misspec_R2(H, ref)` | 1 − SS_res(m)/SS_res(r) | "R² of model improvement over reference" |
| `power_loss(H, ref)` | from `PowerLoss.m` | Lindquist's canonical metric — degrees of effective sample size lost |
| `peak_lag_bias(H, ref)` | argmax_lag H − argmax_lag ref | Where in time is the misspec concentrated? |
| `auc_ratio(H, ref)` | ∫H / ∫ref | Amplitude/AUC bias |
| `curve_sd(H_runs)` | sqrt(var across runs) per lag | Stability of the estimated HRF |

These should all be expressible as voxel‑wise maps **and** as ROI summaries — the misspec map is in itself a useful product (it tells you *where* the model is wrong).

### Cross‑validation for amplitude/AUC

K‑fold CV across runs (or trials within run, for event‑related designs):

```matlab
cv = cv_amplitude(Hstudy, ...
    'KFold', 5, ...
    'Metric', {'amplitude','auc','sqrt_var'}, ...
    'Stratify', 'subject');

cv.amplitude.train_mean   % per fold
cv.amplitude.test_mean
cv.amplitude.generalization_gap
```

Use the test‑fold metric to compare models without circular fits.

### Convergence to an "optimal HRF" (data‑driven, hypothesis‑agnostic)

The pipeline Tor describes is essentially: take all subjects/runs, fit a flexible HRF (sFIR or inverse‑logit), pool across the study with `HMHRFest`‑style mixed effects, get a group HRF, then re‑evaluate misspec metrics against *that* HRF instead of (or in addition to) the canonical. Workflow:

1. `Hstudy_sfir = cat(1, H_run_sfir...)`
2. `H_opt = group_hrf(Hstudy_sfir, 'Method','HMHRFest')` — returns a single group‑level HRF per condition (or per ROI × condition)
3. `metrics = misspec_metrics(Hstudy_canonical, 'Reference', H_opt)` — now you can see how much worse canonical is than the *empirically optimal* HRF, not just how it differs from itself.

### Group t‑tests

```matlab
T = ttest(Hstudy, ...
    'ConditionA','pain', 'ConditionB','neutral', ...
    'Across','subject', ...
    'Lag', 6, ...        % or {} for all lags
    'Metric','amplitude'); % or 'auc', 'peak_value', 'peak_lag'

T.tmap         % statistic_image
T.pmap         % statistic_image
T.summary      % table per ROI
```

This is what the current `hrf_time_unfolding_stats` does, but with a unified object‑oriented entry point.

### Verification (Phase 4 acceptance)

- On synthetic data where the true HRF is known (use `Fit_Logit2` to generate ground truth, then add noise), `misspec_metrics(canonical_fit, 'Reference', true_hrf)` returns power loss within 10% of the analytic value.
- `cv_amplitude` on a real dataset gives `generalization_gap < 0.2 * train_mean` for sFIR fits on n ≥ 20 subjects.
- `group_hrf` converges within 50 iterations on a synthetic dataset and recovers the true HRF to within 5% RMS.

---

## Risks & open issues

- **Folder access for Phase 2.** Until I can see `HRF_Est_Toolbox4` and the canlabCore class folders, every estimate for Phases 2–4 is provisional. Easy to fix — re‑mount one folder up.
- **Naming collisions resolved.** The two-class split (`fmri_hrf` and `statistic_hrf`) avoids the shadow problem that a bare `hrf` class would have caused with SPM/legacy variables. Names fit the existing `fmri_*` / `statistic_*` family in CanlabCore. Confirm with Tor at Phase-3 kickoff.
- **`Old_stuff/*.mat` files.** Some look like reference fixtures used by tutorial/example scripts. Will check uses before deleting.
- **Backward compatibility of `hrf_score_wholebrain_input_table`.** The README documents this entry point. The Phase 1 refactor must keep its signature working, even if the heavy lifting moves to `hrf_score_one_prefix`.
- **SLURM `try/catch` change.** Wrapping per‑task scoring in `try/catch` will mean a task can succeed at the NIfTI level but quietly miss scores. The audit (now with `RepairMissing`) is the safety net for this — both pieces need to ship together.
- **`HMHRFest` performance.** Mixed‑effects HRF estimation at whole‑brain scale is slow. Phase 4 should pilot it on ROIs first and only attempt voxelwise once we know it's tractable.

---

## Sequencing recommendation

1. **This week:** implement Phase 1 against `MichaelSun/HRF_Est_Toolbox2`. Smoke‑test against one finished study. Ship.
2. **Next week:** re‑mount one folder up so I can see Toolbox4 and canlabCore classes; execute Phase 2 (housekeeping + merge into `HRF_Toolbox_Pipeline/`).
3. **Following 2–3 weeks:** Phase 3 OOP refactor, in `HRF_Toolbox_Pipeline/@hrf/`. Iterate naming and field layout with Tor.
4. **Final 2 weeks:** Phase 4 misspec pipeline as `fmri_hrf` methods (with thresholded visualizations on the paired `statistic_hrf`). Validate against synthetic ground truth before applying to real data.

Each phase ships independently — Phase 1 doesn't have to wait for the merge, and Phase 4 doesn't have to wait for everything to be perfect about the object model.
