# RSA/RSM Toolbox: Phased Extension Plan

**Author:** Michael Sun (drafted with Claude, May 2026)
**Audience:** Tor (review), Michael (execution)
**Scope:** New `RSA_tools/` subfolder + `@rsm` class + new methods on `@fmri_data`, `@image_vector`, `@statistic_image`. Generalizes ~12 bespoke RSA workflows (`2026_Sun_et_al_Pain_Imagination/.../Representational Similarity Analysis.mlx`, four `generate_RSA*.m` paradigm wrappers, plus the WASABI / Novus / DistractMap / acceptance studies represented by the dozen `*RSA*.mlx` notebooks in `C:\Temp`) into a single coherent CanlabCore extension.

---

## At a glance

Five phases, each shippable independently:

| Phase | Workstream | Deliverable | Effort | Blocks |
|---|---|---|---|---|
| 1 | `@rsm` class + RSM/RDM construction | Class skeleton; `compute_rsm` with metrics `{correlation, spearman, cosine, euclidean, mahalanobis, crossnobis}` + within-subject whitening; `rsm.from_design`, `rsm.from_categorical`, `rsm.from_table`; `rsm.plot` | Medium (5–7 days) | nothing — ship after sign-off |
| 2 | Cell extraction + group inference + reliability | `rsm.cells`, `rsm.contrast`, `rsm.ttest_contrasts`, `rsm.reliability` (ICC), `rsm.drift`, `rsa_group_inference` | Medium (5–7 days) | Phase 1 |
| 3 | Multi-level / LME modeling | `rsm.rsa_lm` (fixed-fx `fitlm`), `rsm.rsa_lme` (random-fx `fitlme`), `rsm.build_design_rdms`, model comparison + nested LRT helpers | Medium-large (1–1.5 weeks) | Phase 2 |
| 4 | Parcelwise + searchlight + brain maps | `fmri_data.rsa_parcelwise`, `image_vector.searchlight_rsa`, `assign_vals_to_atlas`, `region2statistic_image` | Medium-large (1–1.5 weeks) | Phase 3 |
| 5 | Kriegeskorte interop polish + docs + tests | Soft-dep probe centralization, README, examples (8+ reproducing the dozen workflows), unit tests | Medium (1 week) | Phase 4 |

The phases are sequenced as concentric rings of capability:
- **Phase 1** lets users *build and visualize* RSMs with any metric (including crossnobis).
- **Phase 2** adds *cell-level inference* (within/between contrasts, reliability, drift).
- **Phase 3** adds *multi-level modeling* — the LME machinery that handles "subjects ↑ sessions ↑ runs ↑ conditions ↑ bodysites" structure with random effects.
- **Phase 4** lifts everything to *atlas/searchlight space* and emits CanlabCore-native `statistic_image` maps.
- **Phase 5** polishes the optional Kriegeskorte rsatoolbox bridge and ships docs/tests/examples.

The existing [`@fmri_data/rsa_regression.m`](../@fmri_data/rsa_regression.m) (Kragel generalization-index RSA with within-study bootstrap) stays untouched. Reachable as `rsm.compare(model_rdm, 'method','rsa_regression')` so its bootstrap path stays alive.

---

## Phase 1 — `@rsm` class + RSM/RDM construction

### Diagnosis

Across the workflow corpus there are at least four bespoke RSM constructors (`generate_RSA.m`, `generate_RSA2.m`, `generate_RSA3.m`, `generate_RSA_accept_crossnobis.m`) plus the inline construction in `Representational Similarity Analysis.mlx`, `01282025 RSM Reliability`, `12111024 WASABI RSA Session Effects`, etc. They overlap heavily but each adds one twist:

| Workflow | What it adds beyond `corr(dat.dat)` |
|---|---|
| `generate_RSA.m` (distraction) | Metadata-aware condition aggregation (`mean(get_wh_image(...))`), trial-collapsing (`_trial` suffix), level=`collapsed`/`session`, optional ROI/GM masking + smoothing, `info` struct, matched-pair highlighting |
| `generate_RSA3.m` (Sun et al.) | Per-session RSMs, then average across sessions; **diagonal correction** by replacing same-bodysite diagonal entries with mean of across-session same-bodysite cells |
| `generate_RSA_accept_crossnobis.m` | **Crossnobis distance** with two-fold session split, mean-centering, diagonal Ledoit-Wolf whitening from session-difference noise |
| `10252024 RSA Whitening Walkthrough.mlx` | **Within-subject whitening** of vectorized RDMs across replicates via `rsa.stat.covdiag` + inverse-sqrt SVD |
| `01282025 RSM Reliability.mlx` | Stacks per-replicate RSMs into 3D matrix; computes `ICC(2,k)` reliability |

All five variants assume a 3D RSM-per-subject or array-of-RSMs-per-subject-per-replicate. **An `@rsm` class that standardizes the data + metadata + replicate axis** lets each twist become a single option, not a forked function.

### Proposed design

#### 1.1 — The `@rsm` class

```matlab
classdef rsm
properties
    dat              % k x k  or  k x k x N  (similarity OR dissimilarity)
                     %   N = subjects, replicates, or both (see 'level' below)
    is_dissimilarity % logical: true => RDM, false => RSM
    metric           % 'correlation' | 'spearman' | 'cosine' | 'euclidean' |
                     %   'mahalanobis' | 'crossnobis' | ...
    labels           % {k x 1} cellstr — condition names; row/col order of dat
    metadata_table   % table, k rows: per-condition attributes
                     %   (bodysite, side, condition, ...)
    groupings        % struct: name -> indices into 1:k, e.g. groupings.hot=1:8
    level            % 'subject' | 'session' | 'run' | 'collapsed' | 'group'
                     %   describes what each slice of dim 3 represents
    replicate_table  % N x p table describing dim 3:
                     %   (subject_id, session_number, run_number, fold, ...)
                     %   one row per slice along dim 3
    whitened         % struct describing whitening applied:
                     %   .level    'none'|'within_subject'|'across_subject'|'session_difference'
                     %   .shrinkage  scalar
                     %   .method   'covdiag'|'diag'|'none'
    source           % 'fmri_data' | 'parcel:<name>' | 'searchlight:<vox>' | 'custom'
    history          % cellstr of operations applied (chronological)
    additional_info  % struct, free-form
end

methods
    obj = rsm(dat, varargin)
    obj = rsm.from_design(design, varargin)        % static: model RDM from design columns
    obj = rsm.from_categorical(cat_vec, varargin)  % static: same-vs-different RDM
                                                   %   from a categorical/cellstr vector
                                                   %   (wraps rsa.rdm.categoricalRDM)
    obj = rsm.from_table(meta, spec, varargin)     % static: build multiple model RDMs
                                                   %   from a metadata table + spec
end
end
```

`replicate_table` is the key addition over the original plan: it makes the 3rd dimension self-describing. A 3D rsm with 45 slices (9 subjects × 5 runs) is unambiguous because `replicate_table` says which slice is which (subject, session, run). Phase 3's LME relies on this.

#### 1.2 — `rsm.from_categorical`, `rsm.from_metadata_distance`, `rsm.from_design`

Three static constructors for model RDMs.

```matlab
% (A) Categorical / same-vs-different (wraps rsa.rdm.categoricalRDM)
M = rsm.from_categorical(dat.metadata_table.subject_id);

% From multiple columns at once → array of model RDMs (one per column)
M = rsm.from_categorical(dat.metadata_table, ...
    {'subject_id','session_number','run_number','condition','bodysite'});
% M is a [1 x 5] array of rsm objects, each with is_dissimilarity=true

% (B) Continuous metadata distance (NEW per Phase 3 §6.4)
M = rsm.from_metadata_distance(dat.metadata_table, 'session_number', ...
    'metric', 'abs_diff');
% Entry (i,j) = |session_i - session_j|. Lets users model session-distance,
% time-since-baseline, run-distance, etc. as continuous predictors in rsa_lme.

% (C) Numeric/binary design matrix (existing rsa_regression entry point)
M = rsm.from_design(X, 'names', {'hot','warm','imagine'}, ...
                       'metric', 'seuclidean');
```

These replace both the `pdist(design(:,i),'seuclidean')` loop in `rsa_regression.m` lines 95–100 *and* the `DesignRSMs.{Bodysite,Condition,CrossCondition}` ad-hoc construction in the Sun et al. workflow.

#### 1.3 — `fmri_data.compute_rsm` (the omnibus constructor)

Single method that subsumes all four bespoke `generate_RSA*` variants:

```matlab
R = compute_rsm(dat, ...                       % omnibus form
    'group_by',         'condition', ...       % metadata column → k labels (row/col)
    'condition_collapse', 'mean', ...          % how to aggregate multiple rows
                                               %   per condition: 'mean'|'concat'|'none'
    'level',            'subject', ...         % 'subject'|'session'|'run'|'collapsed'
    'subject_var',      'subject_id', ...      % metadata col defining subject axis
    'session_var',      'session_number', ...  % required if level='session'
    'run_var',          'run_number', ...      % required if level='run'
    'fold_var',         'session_number', ...  % required if metric='crossnobis'
    'metric',           'correlation', ...     % see below
    'whiten',           'none', ...            % 'none'|'within_subject'|
                                               %   'session_difference'
    'whiten_method',    'covdiag', ...         % 'covdiag' (Ledoit-Wolf) | 'diag'
    'diagonal_correction', 'none', ...         % 'none'|'across_session_mean'|'nan'
                                               %   (across_session_mean reproduces
                                               %   generate_RSA3 diagonal patch)
    'mask',             [], ...                % image_vector/atlas/region for apply_mask
    'parcellation',     [], ...                % atlas → returns [nParcels x 1] array of rsm
    'smooth_mm',        [], ...                % optional smoothing
    'treat_zero_as_data', false, ...
    'use_parallel',     false);

% Distance metrics implemented:
%   'correlation' / 'spearman'     — corr(dat, 'type', 'Pearson'/'Spearman')
%   'cosine'                       — canlab_compute_similarity_matrix
%   'euclidean' / 'seuclidean'     — pdist
%   'mahalanobis'                  — pdist with regularized covdiag
%   'crossnobis'                   — see 1.4 below
```

`level='session'` returns a 3D rsm where dim 3 indexes (subject, session) pairs and `replicate_table` records which is which. `level='collapsed'` averages across sessions/runs first, then computes one RSM per subject. `level='subject'` is the default (one RSM per subject from all that subject's images).

#### 1.4 — Crossnobis support (`metric='crossnobis'`)

Generalizes `generate_RSA_accept_crossnobis.m`. Requires a `fold_var` (typically `session_number` or `run_number`) that defines the split.

Algorithm (mirroring the existing code):
1. For each fold `f`, build `X_f` of shape `[k × voxels]` by averaging images per condition within that fold.
2. Drop voxels that are non-finite or zero-variance jointly across folds.
3. Mean-center each `X_f` across conditions.
4. Optional whitening (`'whiten','session_difference'`): noise = `X1 - X2`, voxel_var = `var(noise, 0, 1)`, `X_f_white = X_f ./ sqrt(voxel_var)`.
5. Cross-validated dissimilarity: `D(a,b) = mean( (X1[a]-X1[b]) .* (X2[a]-X2[b]) )` (generalizes to >2 folds as the mean over all distinct fold pairs).
6. Return as `rsm` with `is_dissimilarity=true`, `metric='crossnobis'`, `whitened.level='session_difference'` etc.

The fold-builder pattern from `build_accept_runfold_matrices` becomes a private helper `build_fold_pattern_matrices(dat, group_by, fold_var)` that's reusable for any crossvalidated RSA.

#### 1.5 — Within-subject whitening (`'whiten','within_subject'`)

Generalizes the workflow in `10252024 RSA Whitening Walkthrough.mlx`. Three operating modes:

```matlab
% Build per-replicate RSMs first, then whiten the stack
R = compute_rsm(dat, ...
    'level',    'run', ...
    'metric',   'correlation', ...
    'whiten',   'within_subject');     % vectorize upper-tri of each replicate;
                                       % covdiag → svd → inverse-sqrt → re-square
```

Centralized in `@rsm/private/whiten_within_subject.m`. Wraps `rsa.stat.covdiag` when present (soft-dep), falls back to a stock implementation of Ledoit-Wolf shrinkage when not (we own the simple version; it's ~40 lines).

Per the walkthrough's recommendation, **no across-subject whitening by default** — that one is exposed but documented as "use with care; can suppress meaningful inter-subject variability."

#### 1.6 — Diagonal correction (`'diagonal_correction','across_session_mean'`)

Reproduces the `generate_RSA3.m` patch where same-bodysite diagonal entries are replaced with the mean of across-session same-bodysite cells. Generalized: given a metadata grouping column and a session column, replace each diagonal cell with the mean of across-session same-grouping off-diagonal cells.

```matlab
R = compute_rsm(dat, ...
    'level',                'subject', ...
    'diagonal_correction',  'across_session_mean', ...
    'diagonal_group_by',    'bodysite', ...
    'session_var',          'session_number');
```

#### 1.7 — `rsm.plot`

```matlab
plot(R);                   % rank-transformed heatmap (Kriegeskorte style if available)
plot(R, 'mode', 'raw');    % skip rank transform
plot(R, 'subject', 3);     % slice a 3D rsm
plot(R, 'mean');           % average across 3rd dim before plotting
plot(R, 'subjects', 'grid');  % subplot grid (one per subject/replicate)
plot(R, 'mds');            % 2D MDS scatter
plot(R, 'dendrogram');     % hierarchical clustering
plot(R, 'matched_pairs', matched_pairs);   % highlight (i,j) cells with white rectangles
                                           %   (from generate_RSA plot_RSA_matrix)
plot(R, 'block_borders_by', 'bodysite', ...
       'border_color', 'auto');  % color-coded borders around metadata-grouped blocks
                                  % (from 08192024 Representational Drift)
```

Soft-dep probe (`@rsm/private/probe_rsatoolbox.m`) returns a capability struct; `plot` uses `rsa.util.RDMcolormap` + `rsa.rdm.rankTransform` when present, else `parula` + tiedrank.

### What this phase reproduces from the corpus

- **Sun et al. Figs 7A–C** (Representational Similarity Analysis.mlx)
- **All four `generate_RSA*.m` variants** collapse to single `compute_rsm` calls
- **`10252024 RSA Whitening Walkthrough`** within-subject path (`'whiten','within_subject'`)
- **`generate_RSA_accept_crossnobis`** (`'metric','crossnobis','fold_var','session_number'`)
- **`generate_RSA3` diagonal correction** (`'diagonal_correction','across_session_mean'`)
- **`08192024 Representational Drift`** RSM heatmaps with bodysite block borders (`'block_borders_by'`)
- **`01282025 RSM Reliability`** Figs A/B (subject-level + subject-grid RSM plots)

### File layout (Phase 1)

```
CanlabCore/
  @rsm/
    rsm.m                       % classdef + constructor
    from_categorical.m          % static method  (NEW vs. v1 plan)
    from_metadata_distance.m    % static method  (NEW vs. v1 plan; per Phase 3 §6.4)
    from_design.m               % static method
    from_table.m                % static method
    plot.m
    mean.m
    fisher_z.m
    inv_fisher_z.m
    to_rdm.m
    to_rsm.m
    subset.m
    reorder.m
    display.m
    size.m
    isempty.m
    private/
      local_rank_transform.m         % fallback when rsatoolbox absent
      validate_metadata.m
      probe_rsatoolbox.m
      whiten_within_subject.m        % NEW
      whiten_session_difference.m    % NEW (crossnobis helper)
      build_fold_pattern_matrices.m  % NEW (crossnobis helper)
      ledoit_wolf_shrinkage.m        % NEW fallback (used if rsa.stat.covdiag absent)
      apply_diagonal_correction.m    % NEW
  @fmri_data/
    compute_rsm.m               % omnibus constructor — subsumes generate_RSA*
  RSA_tools/
    RSA_Pipeline_Phased_Plan.md
    README.md                   % added in Phase 5
```

### Verification (Phase 1 acceptance)

- Each `generate_RSA*.m` produces output that matches a single `compute_rsm` call to floating-point tolerance on a synthetic dataset with planted structure.
- Crossnobis distance: on a dataset where condition A and B differ by a known amplitude in signal but share noise, recovered `D(A,B)` is within 5% of the analytic crossnobis value (Walther et al. 2016).
- Within-subject whitening: reproduces `X_subject_avg` from `10252024 RSA Whitening Walkthrough.mlx` line 180 to floating-point tolerance.
- `plot(R, 'matched_pairs', ...)` renders white rectangles matching `generate_RSA.m`'s `plot_RSA_matrix`.
- `compute_rsm` runs without rsa.* on the path (soft-dep test).

---

## Phase 2 — Cell extraction + group inference + reliability

### Diagnosis

Cell-extraction logic (`process_*_correlations` in the workflow, the inline reshaping in `12111024 WASABI RSA Session Effects`, and the matched-bodysite logic across all WASABI notebooks) is the most-duplicated machinery in the corpus. **Phase 2 collapses it into three declarative methods.**

Reliability and drift were not in v1 of the plan; the corpus has them everywhere:

- `01282025 RSM Reliability Analysis and Visualization` — `ICC(2,k)` per subject across runs, per superordinate condition, per bodysite, per individual condition (row-level).
- `08072024 Run-Level RDM Analysis with RSA Toolbox` — `ICC(3,k)` per subject for each pain signature × condition (NPS, SIIPS).
- `08192024 Representational Drift in RSA` — within-condition pattern correlation vs. session distance.

### Proposed design

#### 2.1 — `rsm.cells` (unchanged from v1)

```matlab
vals = R.cells('hot', 'hot');           % within-condition (upper triangle, excluding diag)
vals = R.cells('hot', 'warm');          % between-condition
T    = R.cells_table({'hot','warm','imagine','hw','hi','iw', ...
                      'ipsi','contra','same_bs','diff_bs'});
```

#### 2.2 — `rsm.contrast` and `rsm.ttest_contrasts` (unchanged from v1)

```matlab
spec = {
    'hot_vs_warm',     {'hot','hot'},        {'warm','warm'};
    'hi_vs_hw',        {'hot','imag'},       {'hot','warm'};
    'ipsi_vs_contra',  {'ipsi','ipsi'},      {'contra','contra'};
    'same_vs_diff_bs', {'matching','matching'}, {'nonmatch','nonmatch'};
};
T = R.ttest_contrasts(spec, 'tail','right', 'within',true, 'fdr',true);
```

#### 2.3 — `rsm.reliability`  *(NEW vs. v1)*

Wraps the `ICC.m` File Exchange utility (already used throughout the corpus); soft-deps when present, falls back to a stock implementation we provide.

```matlab
% Whole-RSM reliability across replicates (the most common form)
[icc_value, ci] = R.reliability('icc_type', '2-k', 'across', 'replicate');

% Per-subject reliability across the subject's replicates (3D rsm with replicate axis)
icc_by_subject = R.reliability('across', 'replicate', 'by', 'subject_id');
% returns a table: subject_id | icc

% Per-grouping reliability — replaces the per-condition / per-bodysite loops in
% 01282025 RSM Reliability lines 60–195
icc_table = R.reliability_by_grouping('groupings', {'hot','warm','imagine', ...
                                                    'leftface','rightface',...});
% returns a table: grouping | icc | ci | mean | sd

% Row-level reliability (per-condition reliability of an RSM row)
icc_per_cond = R.reliability_per_condition('icc_type', '2-k');
% returns one ICC per row of the RSM
```

#### 2.4 — `rsm.drift`  *(NEW vs. v1)*

```matlab
% Within-condition similarity vs session distance (08192024 Representational Drift)
T = R.drift('reference', 'first', ...     % compare each replicate to first replicate
            'group_by',  'condition', ...
            'plot',      true);
% T: replicate_index | distance_from_ref | correlation | condition

% Or against a fitted slope
T = R.drift('reference', 'all', 'fit', 'linear');
% T: condition | slope | intercept | R^2 | p
```

#### 2.5 — `rsa_group_inference` (free function, unchanged from v1)

#### 2.6 — `plot_rsm_contrast_bars` (unchanged from v1)

Add: `plot_rsm_session_lineplot(merged_table, condition, rois, ...)` adapted from `12111024 WASABI RSA Session Effects`'s `plot_roi_lineplots` — generalize so it works on any per-session/per-replicate statistic.

### What this phase reproduces

- **Sun et al. Figs 7D–F** (within/between condition + ipsi/contra + same/diff bodysite bar plots)
- **`01282025 RSM Reliability` Figs 1–4** (ICC by subject / by condition / by bodysite / by individual condition)
- **`08192024 Representational Drift` Figs** (session-distance correlation plots)
- **`12111024 WASABI RSA Session Effects` line plots** (per-ROI per-contrast time series)

### File layout (Phase 2)

```
@rsm/
  cells.m
  cells_table.m
  contrast.m
  ttest_contrasts.m
  reliability.m                % NEW vs. v1
  reliability_by_grouping.m    % NEW vs. v1
  reliability_per_condition.m  % NEW vs. v1
  drift.m                      % NEW vs. v1
  private/
    icc_fallback.m             % NEW (stock ICC(2,k) / ICC(3,k) when File Exchange ICC absent)
RSA_tools/
  rsa_group_inference.m
  plot_rsm_contrast_bars.m
  plot_rsm_session_lineplot.m  % NEW vs. v1
```

### Verification (Phase 2 acceptance)

- `ttest_contrasts` reproduces the workflow's `pairwise_stats_table` to floating-point tolerance.
- `reliability` reproduces `01282025 RSM Reliability` `reliability_superordinate` and `reliability_bodysite` tables to floating-point tolerance.
- `drift` reproduces the bar plots in `08192024 Representational Drift` lines 380–400.
- `plot_rsm_session_lineplot` reproduces `plot_roi_lineplots` output visually.

---

## Phase 3 — Multi-level / LME modeling  *(NEW phase vs. v1)*

> **Detailed design lives in [RSA_Phase3_LME_Design.md](RSA_Phase3_LME_Design.md).** The summary below is enough to scope the phase; the linked doc is the source of truth for the API, table contract, nested-model spec, and validation plan.

### Diagnosis

`08072024 Run-Level RDM Analysis with RSA Toolbox.mlx` is the only place in the corpus that handles the full nested structure (subject ↑ session ↑ run ↑ condition ↑ bodysite). It does this by:

1. Concatenating all images into a single big object.
2. Computing one giant RSM `corr(dat.dat)`.
3. Building categorical model RDMs for each grouping variable (Subject, Session, Run, Condition, Bodysite) via `rsa.rdm.categoricalRDM`.
4. Adding pairwise (Cond×Bodysite, etc.) and 3-way interaction RDMs.
5. Vectorizing the upper triangle of the big RSM as `Y_vectorized`.
6. Vectorizing each model RDM the same way → columns of `all_regressors`.
7. Fitting `fitlm(all_regressors, Y_vectorized, 'VarNames', ...)` for fixed-effects regression.
8. Fitting `fitlme(mdl_tbl_sub, formula, 'FitMethod', 'REML')` for random-effects (random intercepts/slopes by Subject).
9. Comparing nested models (10 LMEs in the workflow) via likelihood ratio / AIC.

**This is the "multi-level modelling across scan runs, sessions, subjects, and tasks" pattern you flagged.** It's also the one piece of the corpus that is genuinely under-tooled in CanlabCore today — there's nothing parallel to this workflow anywhere in the toolbox.

### Proposed design

#### 3.1 — `rsm.build_design_rdms`

Promotes `rsa.rdm.categoricalRDM` calls + interaction-RDM construction into a single method.

```matlab
% Build a design RDM stack from metadata columns + optional interactions
[X, names] = R.build_design_rdms( ...
    'main_effects',     {'subject_id','session_number','run_number','condition','bodysite'}, ...
    'interactions',     {{'condition','bodysite'}, {'session_number','condition'}}, ...
    'three_way',        {{'subject_id','condition','bodysite'}});

% Returns:
%   X       p x m matrix of vectorized model RDMs (one per main effect or interaction)
%   names   {m x 1} regressor names: {'Subject','Session','Run','Condition','Bodysite',
%                                     'CondxBodysite', 'SesxCond', 'SubxCondxBodysite'}
```

Internally calls `rsa.rdm.categoricalRDM` if present, else our own implementation (~10 lines: `cat_vec(:) ~= cat_vec(:)'`).

#### 3.2 — `rsm.rsa_lm` (fixed effects)

```matlab
% Fit fitlm on the vectorized RSM with categorical/numeric model RDMs as predictors
mdl = R.rsa_lm( ...
    'predictors',       {'subject_id','session_number','run_number','condition','bodysite'}, ...
    'interactions',     {{'condition','bodysite'}}, ...
    'response_transform','fisherz', ...     % 'none'|'fisherz'|'rank'
    'standardize',      true);

% Returns a LinearModel — disp(mdl), coefTest, etc. all work.
```

#### 3.3 — `rsm.rsa_lme` (random effects)

The main multi-level entry point.

```matlab
% Random intercept by subject, fixed effects for everything else
mdl = R.rsa_lme( ...
    'fixed_effects',    {'session_number','run_number','condition','bodysite'}, ...
    'random_effects',   {'(1 | subject_id)'}, ...
    'fit_method',       'REML');

% Random slopes for condition by subject:
mdl = R.rsa_lme( ...
    'fixed_effects',    {'condition','bodysite','session_number'}, ...
    'random_effects',   {'(condition | subject_id)'});

% Or pass a Wilkinson formula directly (escape hatch):
mdl = R.rsa_lme('formula', ...
    'similarity ~ Condition + Bodysite + Session + (Condition | Subject)');
% Returns a LinearMixedModel.
```

Both methods wrap a shared internal pipeline:
1. `build_design_rdms` to get the predictor matrix and names.
2. Assemble a `mdl_tbl` table with columns `(similarity, Subject, Session, Run, Condition, Bodysite, ...)` — each row is one (i,j) pair from the vectorized RSM.
3. Hand to `fitlm` or `fitlme` (MATLAB Statistics Toolbox).

#### 3.4 — `rsm.compare_models` (nested LRT)

```matlab
% Compare nested LME models — replaces the manual 10-model loop in
% 08072024 Run-Level RDM Analysis lines 2633–2691
[T, best] = R.compare_models({
    'similarity ~ 1 + (1|Subject)';
    'similarity ~ Condition + (1|Subject)';
    'similarity ~ Condition + Bodysite + (1|Subject)';
    'similarity ~ Condition*Bodysite + (1|Subject)';
    'similarity ~ Condition*Bodysite + Session + (1|Subject)';
    'similarity ~ Condition*Bodysite + Session + (Condition|Subject)';
}, 'criteria', {'aic','bic','loglik','lrt'});

% T: model | aic | bic | loglik | df | lrt_chi2 | lrt_p
```

#### 3.5 — `rsm.rsa_lme_by_parcel`  *(bridges to Phase 4)*

Lifted to atlas level in Phase 4 as `fmri_data.rsa_parcelwise_lme`.

### What this phase reproduces

- **`08072024 Run-Level RDM Analysis` Sections "Mixed Effects Models Instead"** (all 10 nested LME comparisons)
- **`08072024` Sections "fitlm" per-subject + pooled** (rsa_lm and rsa_lme entry points)
- **Categorical RDM construction loops** at lines 1625–1629, 1807–1811 → single `build_design_rdms` call

### File layout (Phase 3)

```
@rsm/
  build_design_rdms.m
  rsa_lm.m
  rsa_lme.m
  compare_models.m
  private/
    vectorize_upper_tri.m       % helper used by all three
    assemble_lme_table.m
    fallback_categoricalRDM.m   % when rsa.rdm.categoricalRDM absent
RSA_tools/
  rsa_lme_formula_parser.m      % helps users assemble Wilkinson formulas from
                                % predictor + interaction specs
```

### Verification (Phase 3 acceptance)

- `rsa_lm` reproduces the `fitlm` outputs in `08072024` lines 1725–1737 (full model + each main effect) to floating-point tolerance.
- `rsa_lme` reproduces `mdl_random_full` from `08072024` line 2602 to within REML convergence tolerance.
- `compare_models` reproduces the nested LRT ordering of the 10 models from `08072024` lines 2633–2691.
- Runs in <5 minutes on a synthetic dataset of size (n_subjects=20, n_sessions=2, n_runs=5, n_conditions=24).

---

## Phase 4 — Parcelwise + searchlight + brain maps

### 4.1 — `fmri_data.rsa_parcelwise` (folds in Phase 3 LME)

```matlab
results = rsa_parcelwise(dat, ...
    'atlas',         canlab2024_atlas, ...
    'group_by',      'condition', ...
    'subject_var',   'subject_id', ...
    'metric',        'spearman', ...           % or 'crossnobis', 'mahalanobis', ...
    'whiten',        'within_subject', ...     % optional
    'contrasts',     spec, ...                 % Phase 2 contrast spec
    'lme',           'similarity ~ Condition + (Condition|Subject)', ...
                                               % optional: per-parcel LME from Phase 3
    'correction',    'fdr', ...
    'fdr_scope',     'contrast', ...           % 'contrast' (per contrast across parcels) | 'all'
    'tail',          'right', ...
    'use_parallel',  true);

% results.per_parcel_rsms    [nParcels x 1] array of rsm objects
% results.contrast_table     long table: parcel | contrast | t | p | FDR_P | sig
% results.lme_table          long table: parcel | term | beta | se | p | FDR_P  (if 'lme')
% results.maps               struct: maps.<contrast_name> = statistic_image
% results.atlas, .history
```

### 4.2 — `assign_vals_to_atlas`, `region2statistic_image`

Unchanged from v1.

### 4.3 — `image_vector.searchlight_rsa`

```matlab
results = searchlight_rsa(dat, model_rdms, ...
    'radius',      6, ...
    'group_by',    'condition', ...
    'subject_var', 'subject_id', ...
    'metric',      'spearman', ...
    'compare',     'spearman', ...        % RSM vs model
    'mask',        gray_matter_mask, ...
    'permutations', 1000);                % permutation test for inference
```

Reuses `@image_vector/searchlight.m`'s sphere infrastructure.

### What this phase reproduces

- **Sun et al. Figs 7G, S8, S9** (parcelwise brain maps)
- **`12111024 WASABI RSA Session Effects`** session-stratified parcelwise maps (via `'lme'` with `Session` as fixed effect)
- Future searchlight RSA (not in any current workflow but the natural progression)

### File layout (Phase 4)

```
@fmri_data/
  rsa_parcelwise.m
@image_vector/
  searchlight_rsa.m
RSA_tools/
  assign_vals_to_atlas.m
  region2statistic_image.m
```

### Verification (Phase 4 acceptance)

- Reproduces `condition_results_summary` and `condition_contrast_summary` tables from Sun et al. workflow.
- `results.maps.<contrast>` is a valid `statistic_image` — `.threshold/.montage/.table` work.
- `searchlight_rsa` permutation p-map has uniform distribution under randomized data.
- Optional `'lme'` argument reproduces `08072024 Run-Level RDM Analysis` per-subject LME results when applied parcelwise.

---

## Phase 5 — Kriegeskorte interop polish + docs + tests

### 5.1 — Centralized soft-dep probe

`@rsm/private/probe_rsatoolbox.m` returns a capability struct:

```matlab
caps.rdm_rankTransform    = ~isempty(which('rsa.rdm.rankTransform'));
caps.rdm_squareRDM        = ~isempty(which('rsa.rdm.squareRDM'));
caps.rdm_categoricalRDM   = ~isempty(which('rsa.rdm.categoricalRDM'));
caps.util_RDMcolormap     = ~isempty(which('rsa.util.RDMcolormap'));
caps.fig_MDSConditions    = ~isempty(which('rsa.fig.MDSConditions'));
caps.fig_dendrogramConditions = ~isempty(which('rsa.fig.dendrogramConditions'));
caps.stat_covdiag         = ~isempty(which('rsa.stat.covdiag'));
caps.compareRefRDM2candRDMs = ~isempty(which('rsa.compareRefRDM2candRDMs'));
caps.bootstrapRDMs        = ~isempty(which('rsa.bootstrapRDMs'));
caps.any                  = any(struct2array(caps));
```

Every method that uses `rsa.*` calls this and branches.

### 5.2 — `rsm.compare` (formal RDM-vs-RDM inference)

Wraps `rsa.compareRefRDM2candRDMs` when present (Spearman/Kendall, signed-rank inference, FDR across candidate models, noise ceilings); falls back to our own stock implementation when not.

```matlab
result = R.compare(model_rdms, ...
    'correlation_type',  'kendall_taua', ...
    'relatedness_test',  'subjectRFXsignedRank', ...
    'differences_test',  'subjectRFXsignedRank', ...
    'multiple_testing',  'fdr', ...
    'noise_ceiling',     true);
```

### 5.3 — README + examples + tests

```
RSA_tools/
  README.md
  examples/
    rsa_quickstart.m                            % 30-line example, public dataset
    rsa_sun2026_replication.m                   % full Sun et al. workflow
    rsa_kragel_generalization.m                 % wraps rsa_regression via @rsm
    rsa_crossnobis_acceptance.m                 % generalize_RSA_accept_crossnobis
    rsa_whitening_walkthrough.m                 % 10252024 walkthrough
    rsa_reliability_icc.m                       % 01282025 reliability
    rsa_drift_session.m                         % 08192024 + 12111024 session/drift
    rsa_multilevel_lme.m                        % 08072024 LME corpus
  tests/
    test_rsm_construction.m
    test_compute_rsm_metrics.m
    test_crossnobis.m
    test_whitening.m
    test_cells_and_contrasts.m
    test_reliability.m
    test_drift.m
    test_lme.m
    test_parcelwise.m
    test_searchlight.m
    test_softdep_fallback.m                     % rsa/ off the path
```

### Verification (Phase 5 acceptance)

- Every example script runs end-to-end on a freshly cloned CanlabCore with no extras.
- `test_softdep_fallback` passes with `rmpath(genpath(<rsatoolbox>))` in effect.
- README covers each metric, each level (subject/session/run), and links each example to its inspiring workflow.

---

## Risks & open issues

- **Open question — replicate axis vs subject axis on 3D `rsm`.** When `level='session'` we have one slice per (subject, session) pair. The replicate axis is now compound. `replicate_table` makes this self-describing, but cell-extraction and reliability methods need to know whether to aggregate within-subject across sessions, across-subject within-session, or pool. Default = follow the `level` field; expose `'aggregate_across', {'subject'|'session'|'run'}` for explicit control. Confirm at Phase 1 kickoff.
- **Crossnobis with >2 folds.** The `generate_RSA_accept_crossnobis` implementation assumes exactly 2 folds (one session pair). Generalize to k folds by averaging cross-products over all distinct pairs (Walther et al. 2016 eq. 5). Need a leave-one-fold-out variant too — confirm preference.
- **Whitening interaction with crossnobis.** Crossnobis already includes its own whitening (session-difference noise). Within-subject covdiag whitening on top of that is double-counting noise. The `compute_rsm` validator should warn if both are requested.
- **LME on full RSM scales poorly.** A k=200 condition RSM has 19,900 cells, and stacking subjects/sessions/runs grows that by 10–100×. `fitlme` will struggle past ~10⁶ rows. Phase 3 must document this and offer a `'subsample'` option (random or stratified) for prototyping.
- **ICC fallback.** The File Exchange `ICC.m` is the de-facto standard. Our `icc_fallback.m` must give numerically identical output for ICC(2,k) and ICC(3,k) — write unit tests against `ICC.m` output before declaring ready.
- **`fmri_data.compute_rsm` API breadth.** The omnibus form is verbose. We'll provide three short aliases (`compute_rsm_per_subject`, `compute_rsm_per_session`, `compute_rsm_crossnobis`) that pre-fill `level` and the relevant defaults. Each is a 5-line wrapper.
- **Atlas-vs-mask-vs-region inputs to ROI methods.** The `generate_RSA*.m` wrappers switch on `class(roi)`. Standardize: any method that takes a spatial restriction accepts `image_vector | atlas | region | char (atlas keyword)` via a single helper.
- **Replication coverage of the corpus.** I've inventoried 12 `.mlx` files plus four `.m` paradigm wrappers and the published Sun et al. workflow. Files I haven't yet read deeply (DistractMap, Novus, accept paradigm publication figures, DEMO1/DEMO2, `b_imensa_240116`, `12222023`) may reveal additional patterns. Phase 5 example coverage is the safety net — if an example can't be written with the current API, the API is incomplete.
- **MATLAB MCP integration for testing.** The MATLAB MCP server (`mcp__MATLAB__run_matlab_file`, `mcp__MATLAB__evaluate_matlab_code`, `mcp__MATLAB__check_matlab_code`) now lets me run example scripts and unit tests live during development. Use it for Phase 1 smoke tests before declaring `compute_rsm` ready.

---

## Sequencing recommendation

1. **Weeks 1–2:** Phase 1. Ship after `compute_rsm` reproduces all four `generate_RSA*` variants via single calls.
2. **Week 3:** Phase 2. Validate ICC/drift against `01282025` and `08192024` reference outputs.
3. **Weeks 4–5:** Phase 3 (the new multi-level phase). Validate LMEs against `08072024` reference outputs. This is the largest single phase by scope.
4. **Weeks 6–7:** Phase 4 (parcelwise + searchlight + maps).
5. **Week 8:** Phase 5 (polish + docs + tests).

Each phase remains independently shippable and mergeable. Phase 1 alone gives users a coherent RSM-construction API with crossnobis and whitening — that's already a substantial uplift over the status quo.
