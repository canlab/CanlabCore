# CanlabCore RSA / RSM Tools

First-class Representational Similarity Analysis for CanlabCore. Build
representational similarity/dissimilarity matrices from `fmri_data` objects,
run the full inference battery (reliability, within/between contrasts,
drift, multi-level LME, formal RDM comparison), and project results onto the
brain as `statistic_image` maps.

Generalizes ~12 bespoke RSA workflows (Sun et al. 2026 pain-imagination
study and the WASABI / Novus / DistractMap / acceptance studies) into a
single coherent, reusable extension.

---

## Install / setup

The tools live inside CanlabCore. Just add CanlabCore to your path:

```matlab
addpath(genpath('/path/to/CanlabCore'));
```

Check optional dependencies (everything works without them via built-in
fallbacks):

```matlab
rsa_toolbox_status
```

| Dependency | Used for | Fallback if absent |
|---|---|---|
| Kriegeskorte **rsatoolbox** | rank-transformed plots, MDS, `compareRefRDM2candRDMs` | built-in plotting + `rsm.compare` engine |
| **ICC.m** (File Exchange) | reliability | built-in ICC(2,k)/ICC(3,k) |
| **Statistics Toolbox** `fitlme` | multi-level LME | **required** for LME methods |
| **Statistics Toolbox** `rangesearch` | fast searchlight | O(n²) loop |

> **Reloading after edits:** if you change any `@rsm`, `@fmri_data`, or
> `@image_vector` method, run `clear classes; rehash path` and re-cast any
> existing data objects: `dat = fmri_data_st(struct(dat));`

---

## The `@rsm` object

The container for similarity (RSM) or dissimilarity (RDM) matrices:

```
dat              k×k or k×k×N matrix (N = subjects/sessions/runs)
is_dissimilarity true => RDM, false => RSM
metric           'correlation' | 'spearman' | 'crossnobis' | ...
labels           {k×1} condition names
metadata_table   k-row per-condition attributes
groupings        struct: name -> indices (auto-built by compute_rsm)
replicate_table  N-row description of the 3rd dim
level            'subject' | 'session' | 'run' | 'collapsed'
additional_info  .qc diagnostics (missing-condition report)
```

Value semantics: every method returns a new `rsm`; nothing mutates in place.

---

## Quick start

```matlab
% 1. Build a per-subject RSM (24 conditions = 3 conditions x 8 bodysites)
R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, ...
                     'subject_var', 'sub', 'metric', 'spearman');

% 2. Visualize the group-mean geometry
figure; plot(mean(R), 'block_borders_by', 'condition');

% 3. Within vs between condition contrasts (FDR-corrected)
spec = { 'within_hot', 'hot', [];  'hot_vs_warm', 'hot', 'warm' };
T = R.ttest_contrasts(spec, 'tail', 'right', 'correction', 'fdr');

% 4. Reliability across sessions (per-subject ICC, then aggregate)
R_sess = compute_rsm(dat, 'group_by', {'condition','bodysite'}, ...
                          'level', 'session', 'subject_var', 'sub', ...
                          'session_var', 'sesno', 'metric', 'spearman');
T_rel = R_sess.reliability_by_grouping;

% 5. Multi-level LME
mdl = dat.rsa_lme('predictors', {'condition','bodysite','sesno'}, ...
                  'subject_var', 'sub');

% 6. Formal RDM comparison to model RDMs
result = R.compare({'condition','bodysite'}, 'correlation_type','kendall_taua');

% 7. Parcelwise brain maps
atlas = load_atlas('canlab2024');
res = dat.rsa_parcelwise('atlas', atlas, 'group_by', {'condition','bodysite'}, ...
                         'subject_var', 'sub', 'contrasts', spec);
res.maps.hot_vs_warm.threshold(0.05, 'unc').montage;
```

---

## Method reference

### Construction (Phase 1)

| Call | Purpose |
|---|---|
| `compute_rsm(dat, ...)` | Omnibus RSM builder. Metrics: `correlation, spearman, cosine, euclidean, mahalanobis, crossnobis`. Levels: `subject, session, run, collapsed`. Options: `whiten, diagonal_correction, parcellation, mask, nan_policy`. |
| `rsm.from_categorical(meta, cols)` | Same-vs-different model RDM(s) from metadata columns. |
| `rsm.from_metadata_distance(meta, col)` | Continuous-distance model RDM (e.g. session distance). |
| `rsm.from_design(X, ...)` | Model RDM(s) from a design matrix. |
| `R.plot(...)` | Heatmap (raw/rank), MDS, dendrogram, grid; block borders, matched-pair overlays. |
| `R.mean / fisher_z / to_rdm / to_rsm / subset / reorder` | Transforms. |
| `R.get_by_label(name)` | Look up a parcel in an rsm array by name. |

### Cell extraction & inference (Phase 2)

| Call | Purpose |
|---|---|
| `R.cells(A, B)` | Per-replicate scalar from within/between cells (Fisher-z mean). |
| `R.cells_table(spec)` | Multiple groupings -> per-replicate table. |
| `R.contrast(A, [B])` | Contrast scalar per replicate. |
| `R.ttest_contrasts(spec, ...)` | Battery of paired/one-sample t-tests + FDR. |
| `R.reliability(...)` | ICC across replicates (per-subject by default). |
| `R.reliability_by_grouping / reliability_per_condition` | Per-group / per-row ICC. |
| `R.drift(...)` | Within-condition stability over the replicate axis + linear fit. |
| `rsa_group_inference(matrix, ...)` | Standalone group stats for raw per-subject scalars. |
| `plot_rsm_contrast_bars(T_or_R, ...)` | Bar plot with within-subject lines + sig markers. |

### Multi-level modeling (Phase 3)

| Call | Purpose |
|---|---|
| `dat.rsa_lme(formula_or_nv, ...)` | Random-effects LME on within-subject RDM cell pairs. |
| `dat.rsa_lm(...)` | Fixed-effects model (all pairs incl. cross-subject). |
| `dat.rsa_lm_by_subject(...)` | Per-subject fits + partial R², for 2nd-level inference. |
| `dat.rsa_compare_models(formulas, ...)` | Nested LRT + AIC/BIC ladder (auto-ML). |
| `rsa_model_sequence(base)` | Build the nested-formula list imperatively. |
| `rsa_lme_icc(mdl) / rsa_lme_blups(mdl)` | Variance components / per-subject BLUPs. |
| `rsa_partial_r2(mdl, tbl)` | Per-predictor partial R² + Cohen's d. |
| `assemble_lme_table(dat, ...)` | The long-format pair table (engine; also exported). |

Interaction columns are named by joining the main-effect columns with `x`,
e.g. `SameConditionxSameBodysite`. Either naming form is accepted in
formulas; a clear error lists available columns on a typo.

### Brain maps (Phase 4)

| Call | Purpose |
|---|---|
| `dat.rsa_parcelwise('atlas', a, 'contrasts', spec, ...)` | Per-parcel contrasts -> `statistic_image` per contrast, FDR across parcels. |
| `dat.rsa_parcelwise('atlas', a, 'lme', formula, ...)` | Per-parcel LME -> `statistic_image` per term. |
| `dat.searchlight_rsa(model, ...)` | Searchlight RSM-vs-model correlation maps; group/subject level, permutation inference. |
| `assign_vals_to_atlas(atlas, names, vals, ...)` | Project per-parcel values onto a brain map. |

### Formal RDM comparison (Phase 5)

| Call | Purpose |
|---|---|
| `R.compare(models, ...)` | Nili et al. (2014) framework: per-subject correlation (Kendall τ-a / Spearman / Pearson), relatedness test (subject RFX signed-rank), pairwise differences test, noise ceiling, FDR. |
| `rsa_toolbox_status` | Report optional-dependency availability. |

---

## Workflow → method map

| Original workflow | Now |
|---|---|
| `generate_RSA*.m` (4 variants) | one `compute_rsm` call |
| Sun et al. Figs 7A–C | `compute_rsm` + `plot` |
| Sun et al. Figs 7D–F | `ttest_contrasts` + `plot_rsm_contrast_bars` |
| Sun et al. Figs 7G / S8 / S9 | `rsa_parcelwise` + `.montage` |
| `01282025 RSM Reliability` | `reliability_by_grouping` |
| `08192024 Representational Drift` | `R.drift` |
| `10252024 Whitening Walkthrough` | `compute_rsm(..., 'whiten','within_subject')` |
| `generate_RSA_accept_crossnobis` | `compute_rsm(..., 'metric','crossnobis')` |
| `08072024 Run-Level RDM Analysis` (LME) | `rsa_lme` + `rsa_compare_models` |
| RSAToolbox recipe (RDM comparison) | `R.compare` |

---

## Notes on statistics

- **Reliability ICC** uses pairs of conditions as items, replicates as
  raters. Per-subject ICC (across that subject's sessions) is the default;
  the pooled-across-all-replicates form is `'pool','replicate'`. Groupings
  with <5 cells (e.g. single-bodysite) produce unreliable ICCs — flagged
  with a warning.
- **LME rows** are within-subject condition-pair similarities; `Subject` is
  the random-effects grouping. Cross-subject pairs only enter `rsa_lm`.
- **compare** works in dissimilarity space; positive correlation = the
  brain's geometry matches the model. Kendall τ-a is the default per Nili
  et al. (2014) for categorical models with tied ranks.
- **searchlight** converts each sphere's similarity RSM to dissimilarity
  before comparing to model RDMs.

---

## See also

- `RSA_Pipeline_Phased_Plan.md` — the design document.
- `RSA_Phase3_LME_Design.md` — LME table contract and API spec.
- `examples/` — runnable example scripts.
- `tests/` — unit tests.
