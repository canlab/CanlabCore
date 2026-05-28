# HRF Toolbox Pipeline — Architecture & Code Map

How the pipeline is laid out, what calls what, and where to add new functionality. For usage examples, see [Tutorial.md](Tutorial.md). For the four-phase refactor history, see [HRF_Pipeline_Phased_Plan.md](HRF_Pipeline_Phased_Plan.md).

---

## 1. Two-zone layout

After the Phase 2 merge, all HRF code lives under `HRF_Est_Toolbox4/`:

```
HRF_Est_Toolbox4/                            (~28 files at root)
│
├── Fit_Canonical_HRF.m                      Lindquist-lab canonical HRF fitters
├── Fit_sFIR.m, Fit_sFIR_epochmodulation.m   ───────────────────────────────────
├── Fit_Spline.m, Fit_Logit2.m, ...          legacy single-voxel HRF estimators
├── Anneal_Logit.m, Det_Logit.m              (called by the pipeline orchestrator)
├── ResidScan.m, PowerLoss.m, PowerSim.m
├── HMHRFest.m                               hierarchical mixed-effects HRF est.
├── hrf_fit_one_voxel.m                      dispatcher: routes (FIR/IL/CHRF) -> fitters
├── New_fminsearchbnd/                       3rd-party bounded fminsearch
├── timecourse.mat                           legacy demo fixture
└── HRF_Toolbox_Pipeline/                    user-facing pipeline + new OOP
    │
    ├── run_hrf_pipeline.m                   single-subject, single-run entry
    ├── run_hrf_study_pipeline.m             multi-subject loop
    ├── hrf_*.m                              30 pipeline functions
    ├── plot_hrf_*.m                         5 plotter functions
    ├── make_fmri_stat_hrf.m                 pairing helper for the new classes
    ├── sFIR_multisubject_example_script.m   tutorial-by-script
    ├── hrf_smoke_test_phase1.m              Phase 1 audit/score smoke test
    ├── hrf_smoke_test_phase3.m              Phase 3 OOP smoke test
    │
    ├── @fmri_hrf/                           NEW (Phase 3 v0)
    │   ├── fmri_hrf.m                       classdef + constructor + disp
    │   ├── cat.m                            HRF-aware (subject, run) concatenation
    │   └── to_statistic_hrf.m               build paired statistic_hrf from beta + SE
    │
    ├── @statistic_hrf/                      NEW (Phase 3 v0)
    │   ├── statistic_hrf.m                  classdef + constructor + disp
    │   └── cat.m
    │
    ├── EstHRF_inAtlas/                      atlas/ROI HRF analysis helpers
    │
    └── docs/                                this folder
        ├── HRF_Pipeline_Phased_Plan.md
        ├── Tutorial.md
        ├── Architecture.md                  (you are here)
        └── HRF Analysis Code Map.pptx       (slide-deck companion to this file)
```

**Zone separation matters:** the Lindquist fitters at the root are stable, mostly-untouched canonical implementations. The pipeline subfolder is where Michael Sun's work lives, and where new methods get added.

---

## 2. Layer architecture

```
┌──────────────────────────────────────────────────────────────────────┐
│ L0  ENTRY POINTS                                                     │
│     run_hrf_pipeline                  single subject/run             │
│     run_hrf_study_pipeline            multi-subject loop             │
│     hrf_write_slurm_study_script      cluster array job              │
└─────────────────────────┬────────────────────────────────────────────┘
                          │
┌─────────────────────────▼────────────────────────────────────────────┐
│ L1  FITTING ORCHESTRATION                                            │
│     hrf_fit_all_models                dispatches per-voxel calls     │
│     hrf_fit_wholebrain_stats          whole-brain across voxels      │
│     hrf_load_events_tsv               BIDS events -> struct          │
│     hrf_build_stick_functions         events -> stick matrix         │
│     hrf_extract_timeseries_from_nii   NIfTI -> z-scored timeseries   │
└─────────────────────────┬────────────────────────────────────────────┘
                          │
┌─────────────────────────▼────────────────────────────────────────────┐
│ L2  LINDQUIST-LAB FITTERS  (parent folder HRF_Est_Toolbox4/)         │
│     Fit_sFIR, Fit_Canonical_HRF, Fit_Spline, Fit_Logit2,             │
│     Fit_NLgamma, Anneal_Logit, Det_Logit, ResidScan, PowerLoss       │
│     hrf_fit_one_voxel  (dispatcher: 'IL' | 'FIR' | 'CHRF')           │
└──────────────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────┐
│ L3  WHOLE-BRAIN I/O                                                  │
│     hrf_load_wholebrain_stats         load beta/t/se/p NIfTIs        │
│     hrf_collect_wholebrain_outputs    index SLURM outputs -> table   │
└─────────────────────────┬────────────────────────────────────────────┘
                          │
┌─────────────────────────▼────────────────────────────────────────────┐
│ L4  SCORING & SIGNATURE APPLICATION                                  │
│     hrf_score_one_prefix              ★ single source of truth       │
│     hrf_apply_maps_to_wholebrain      generic signature scorer       │
│     hrf_score_wholebrain_input_table  post-hoc backfill wrapper      │
│     hrf_extract_signature_timeseries  multi-signature extraction     │
│     hrf_extract_imageset_timeseries   custom image-set application   │
│     hrf_extract_roi_timeseries        ROI averaging                  │
│     hrf_extract_all_signature_timeseries                             │
└─────────────────────────┬────────────────────────────────────────────┘
                          │
┌─────────────────────────▼────────────────────────────────────────────┐
│ L5  AUDIT / REPAIR                                                   │
│     hrf_audit_slurm_outputs           ★ also drives RepairMissing    │
└─────────────────────────┬────────────────────────────────────────────┘
                          │
┌─────────────────────────▼────────────────────────────────────────────┐
│ L6  STUDY ASSEMBLY                                                   │
│     hrf_input_table_to_study          input table -> study struct    │
│     hrf_second_level_inputs_to_study  ... with per-condition contr.  │
└─────────────────────────┬────────────────────────────────────────────┘
                          │
              ┌───────────┴───────────────┐
              ▼                           ▼
┌──────────────────────────┐   ┌──────────────────────────────────────┐
│ L7  ANALYTICS            │   │ L8  VISUALIZATION                    │
│ hrf_time_unfolding_stats │   │ plot_hrf_by_condition                │
│ hrf_analyze_second_level │   │ plot_hrf_study_by_subject  (largest) │
│ hrf_2x2_study_score_stats│   │ plot_hrf_2x2_study_score_stats       │
│ hrf_average_condition_tr.│   │ plot_hrf_time_unfolding_stats        │
│ hrf_compare_conditions   │   │ plot_hrf_results                     │
│ hrf_qc_study_curves      │   │ hrf_animate_wholebrain_stats         │
│ hrf_summarize_study      │   │ hrf_make_montage_animation           │
└──────────────────────────┘   │ hrf_make_average_montage_animations  │
                               └──────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────┐
│ L9  OOP WRAPPERS  (Phase 3 v0)                                       │
│     @fmri_hrf < fmri_data                                            │
│     @statistic_hrf < statistic_image                                 │
│     make_fmri_stat_hrf      canonical paired constructor             │
│                                                                      │
│     Consumes from L3 (loaded objects) and L6 (study struct).         │
│     Built so Phase 4 misspec methods slot in as @fmri_hrf methods.   │
└──────────────────────────────────────────────────────────────────────┘
```

The `★` markers point at modules that consolidate dispersed logic — `hrf_score_one_prefix` is called from the SLURM worker, the audit-repair path, *and* the post-hoc backfill, so all three paths produce byte-identical score CSVs. `hrf_audit_slurm_outputs` is the only place that reasons about completeness.

---

## 3. Call graph (cross-layer arrows only)

```
L0 ─┬─→ L1 ─→ L2                       fit a subject/run
    ├─→ L4 (optional, after fit)       score wholebrain maps
    ├─→ L1 + writes a worker for L2    SLURM array job generation
    └─→ L4 (worker scoring)            try/catch wrapped

After the SLURM job lands its NIfTIs:
L3 ───→ tables / fmri_data objects
  └─→ L5 (audit)
       └─→ L4 (repair via hrf_score_one_prefix)
  └─→ L6 (build study struct)
       ├─→ L7 (analytics)
       ├─→ L8 (visualization)
       └─→ L9 (OOP wrap via make_fmri_stat_hrf)
```

The pipeline is largely DAG-shaped, not cyclic. The one back-edge: L5 (audit) repairs by calling L4 (score one prefix). Everything else flows downstream.

---

## 4. Data flow

```
INPUT
    fmri.nii.gz        ──► fmri_data()              [CanlabCore parent]
    events.tsv         ──► hrf_load_events_tsv()    ──► stick functions
                                                          (volumes × conditions)

FIT (per voxel)
    timecourse + sticks ──► hrf_fit_all_models()
                              ├─► Fit_sFIR(...)            ┐
                              ├─► Fit_Canonical_HRF(...)   │  Lindquist fitters
                              ├─► Fit_Spline(...)          │  (parent folder)
                              ├─► Fit_Logit2(...)          │
                              └─► Fit_NLgamma(...)         ┘
                          ──► .hrf (time × cond)
                              .hrf_se, .fit, .residuals
                              .amplitude, .peak_lag
                              .metadata_table (one row per cond × lag)

STACK across voxels ──► statistic_image objects:
                            .b   (beta, voxels × volumes)
                            .t   (t-stat)
                            .ste (standard error)
                            .p   (p-value)

WRITE TO DISK (whole-brain)
    <prefix>_beta.nii     <prefix>_t.nii
    <prefix>_se.nii       <prefix>_p.nii
    <prefix>_metadata.csv          (volume → cond × lag mapping)

SCORE (optional, immediately or via audit-repair)
    + signature sets / image sets
        ──► hrf_score_one_prefix()
              ──► <prefix>_<object>_map_scores.csv

SECOND-LEVEL ASSEMBLY
    hrf_collect_wholebrain_outputs() ──► input_table  (one row per (subj, run, model))
    hrf_input_table_to_study()      ──► study struct  (in-memory; LoadWholeBrain optional)

GROUP-LEVEL OUTPUTS
    hrf_time_unfolding_stats(study)  ──► per-timepoint t-tests
    hrf_make_average_montage_animations(input_table)
        ──► subject avg → group avg ──► thresholded t-map per condition
        ──► MP4 animation over HRF lags

OOP WRAPPING (for any of: result struct, result.mat path, NIfTI prefix)
    make_fmri_stat_hrf(source, ...) ──► [Hb (fmri_hrf), Ht (statistic_hrf)]
```

---

## 5. Dependencies

### canlabCore parent classes

| Parent | What pipeline uses it for |
|---|---|
| `@fmri_data` | All loaded fMRI/whole-brain data. Subclassed by `@fmri_hrf`. Methods used: constructor, `cat`, `apply_mask`, `extract_roi_averages`, `mean`, `write`, `resample_space`. |
| `@statistic_image` | All whole-brain inferential maps. Subclassed by `@statistic_hrf`. Methods used: constructor, `threshold`, `multi_threshold`, `orthviews`, `convert2mask`, `conjunction`. |
| `@image_vector` | Indirect (base of both above). Methods used: `apply_mask`, `apply_atlas`, `apply_parcellation`, `remove_empty`, `replace_empty`. |
| Free functions | `load_image_set`, `apply_all_signatures`, `apply_mask` (function form), `get_wh_image`, `tor_make_deconv_mtx3`. |

The Phase 3 OOP refactor is engineered so all of these inherit unchanged — `@fmri_hrf` and `@statistic_hrf` don't reimplement any of it. The win: thresholding, masking, ROI extraction, and `orthviews` work on the new classes the day they ship.

### Lindquist-lab fitters

All Fitters are called from **exactly one** dispatcher (`hrf_fit_all_models` for the pipeline, `hrf_fit_one_voxel` for direct single-voxel calls):

| Fitter | Dispatcher entry-point |
|---|---|
| `Fit_sFIR` (and `Fit_sFIR_epochmodulation`) | `method='FIR'` or `'sFIR'` in `hrf_fit_one_voxel`; routed via `hrf_fit_all_models` for the pipeline |
| `Fit_Canonical_HRF` | `method='CHRF'` |
| `Fit_Logit2` (with `Anneal_Logit` / `Det_Logit` / `Get_Logit`) | `method='IL'` |
| `Fit_Spline` | direct in `hrf_fit_all_models` |
| `Fit_NLgamma` | direct in `hrf_fit_all_models` |
| `ResidScan` / `PowerLoss` / `PowerSim` | model-quality utilities (post-fit) |
| `HMHRFest` | hierarchical mixed-effects HRF — used by Phase 4 misspec pipeline (not by Phase 1–3) |

The batch variants at root (`Fit_sFIR_all`, `Fit_Canonical_HRF_all`, `*_all_AR.m`) are **not called by the pipeline** — they're legacy multi-subject wrappers preserved for backward compatibility with old scripts. The pipeline fans out via SLURM array jobs instead.

### Third-party

- SPM (`spm_get_bf`, `spm_hrf`) — required by `Fit_Canonical_HRF` and `Fit_NLgamma`.
- FDA package (`create_bspline_basis`, `eval_basis`) — required by `Fit_Spline`. Error message points at the GitHub repo if missing.
- `New_fminsearchbnd/` — vendored bounded `fminsearch`; used by `Fit_NLgamma`.

---

## 6. Top 10 load-bearing files (by line count)

| File | Lines | Role |
|---|---|---|
| `plot_hrf_study_by_subject.m` | 906 | Multi-panel per-subject HRF plots; biggest single visualization function. |
| `hrf_make_average_montage_animations.m` | 841 | Group-level animated montages; per-prefix load cache (5bd381ae). |
| `hrf_audit_slurm_outputs.m` | 810 | Audit + RepairMissing; orchestrates score backfill. |
| `hrf_score_one_prefix.m` | 695 | **Single source of truth** for per-prefix scoring. 3 callers. |
| `hrf_second_level_inputs_to_study.m` | 690 | Study assembly with per-condition contrasts. |
| `hrf_fit_wholebrain_stats.m` | 540 | Whole-brain voxelwise fitting. |
| `plot_hrf_by_condition.m` | 482 | Per-condition curves with SE bands. |
| `hrf_write_slurm_study_script.m` | 480 | SLURM array job + worker generation. |
| `run_hrf_pipeline.m` | 426 | Single-subject entry point. |
| `hrf_time_unfolding_stats.m` | 424 | Group-level timepoint-wise t-tests. |

The top 4 are the modules to study first if you're trying to understand or extend the pipeline.

---

## 7. Where to add things

| If you want to add… | Put it here |
|---|---|
| A new HRF fitter | `HRF_Est_Toolbox4/Fit_<name>.m` (parent folder, alongside the Lindquist fitters). Then add a case in `hrf_fit_all_models`. |
| A new pipeline orchestration step | `HRF_Toolbox_Pipeline/hrf_<verb>_<noun>.m`. Follow the naming convention. |
| A new visualization | `HRF_Toolbox_Pipeline/plot_hrf_<thing>.m` (procedural) OR a method on `@fmri_hrf` / `@statistic_hrf` (OOP, preferred for new work). |
| A new signature / image set | This goes in CanlabCore proper (`load_image_set`'s registry), not here. |
| A new analysis (group-level test, etc.) | If it takes a study struct: `HRF_Toolbox_Pipeline/hrf_<analysis>.m`. If it takes an `fmri_hrf`: method under `@fmri_hrf/`. The Phase 4 misspec pipeline goes the second route. |
| A new OOP method | Inside the class folder: `@fmri_hrf/<method>.m` or `@statistic_hrf/<method>.m`. Use parent's method via direct field access — don't try to downcast through the parent constructor (it expects a filename as first arg). |

---

## 8. Conventions

- **Naming.** All pipeline functions start with `hrf_` (snake_case verb_noun). All plotters start with `plot_hrf_`. Class folders use `@<classname>/`. Helper-style classes are lowercase (`fmri_hrf`, `statistic_hrf`).
- **Prefix vs model_prefix.** For single-model fits, they're equal. For multi-model fits, `model_prefix = <output_prefix>_<model>` and that's where NIfTIs and metadata are written. The score CSV filename follows the same rule. The helper `local_image_prefix(prefix, model_name, n_models)` in `hrf_score_one_prefix.m` is the canonical resolver.
- **Metadata table.** Every whole-brain output (in-memory or on-disk) is paired with a metadata table — one row per 4D volume — with columns `volume_index`, `condition`, `lag_index`, `lag_seconds`, `image_label`, optionally `N`, `dfe`, `TR`, `mode`. The OOP classes carry this table as a property.
- **Three-way scoring consolidation.** `hrf_score_one_prefix.m` is invoked from (1) the SLURM worker, (2) the audit's RepairMissing path, and (3) the post-hoc `hrf_score_wholebrain_input_table` wrapper. All three produce byte-identical outputs. Any future scoring optimization (caching, batching, GPU) goes in the helper and benefits all three paths.

---

For walkthrough-style examples of using these layers, see [Tutorial.md](Tutorial.md).
