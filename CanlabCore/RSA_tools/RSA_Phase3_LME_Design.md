# RSA Phase 3 — Multi-Level / LME Design Specification

**Author:** Michael Sun (drafted with Claude, May 2026)
**Audience:** Tor (review), Michael (execution)
**Status:** Design doc for review before any Phase 3 implementation begins
**Parent plan:** [RSA_Pipeline_Phased_Plan.md](RSA_Pipeline_Phased_Plan.md)

This document fleshes out the multi-level modeling phase. The other four phases are stable; this one was promoted from "implicit" to its own phase after reviewing the corpus and is large enough to warrant a dedicated spec.

The reference workflow is `C:\Temp\08072024 Run-Level RDM Analysis with RSA Toolbox.mlx`. The contract below mirrors what that workflow actually does, generalized to be paradigm-agnostic.

---

## 1. The contract — what an LME row represents

The single design decision that resolves everything else:

> **One LME row = one (i, j) upper-triangle pair from a single subject's image-level RSM.**

Cross-subject pairs (where image i is from subject A and image j is from subject B) are **excluded** from the LME table. This is what makes `Subject` a clean random-effects grouping variable: every row belongs unambiguously to one subject, so `(1 | Subject)` and `(Condition | Subject)` random-effects terms have well-defined semantics.

Cross-subject pairs are still available — they feed the *fixed-effects* `rsa_lm` model (which has `SameSubject` as one of its categorical predictors). They just don't participate in the LME.

### Schematic of the LME table

For a study with N_s subjects, where subject *s* has n_s images, the LME table has:

```
∑_s [n_s × (n_s − 1) / 2]   rows
```

For the WASABI corpus reference (9 subjects × ~135 images each = ~1215 images total): **per-subject pair count ≈ 9050, summed ≈ 81,000 rows**. `fitlme` handles this comfortably.

Columns:

| Column | Type | Meaning |
|---|---|---|
| `Y` | double | Similarity at pair (i, j); Fisher-z optionally applied |
| `Subject` | categorical | Subject ID — random-effects grouping |
| `<main effect 1>` | binary 0/1 | 1 iff i and j share that level (e.g. `SameSession`) |
| `<main effect 2>` | binary 0/1 | same |
| `...` | | |
| `<interaction>` | binary 0/1 | Element-wise AND of constituent main effects |

Conventionally we drop the `Same` prefix — `Session` is shorthand for "i and j share a session." Documented in the README.

---

## 2. Public API surface

### 2.1 — Top-level entry points

Two equivalent ways to fit an LME, depending on what the user already has:

```matlab
% Path A: one-shot from an fmri_data object
mdl = rsa_lme(dat, formula, varargin)         % @fmri_data method
% Builds per-subject image-level RSMs internally, then assembles the LME table.

% Path B: from a precomputed image-level rsm
mdl = rsa_lme(R, formula, varargin)           % @rsm method
% R must have been produced with level='image' (one image per row, no condition collapsing).
```

Both share a common internal helper `assemble_lme_table` (see §3).

### 2.2 — Common varargin (both paths)

```matlab
%   'predictors'        cellstr of metadata columns to include as main-effect RDMs.
%                       Default: all categorical columns in metadata_table except
%                       the subject_var.
%                       Example: {'session_number','run_number','condition','bodysite'}
%
%   'interactions'      cell of cellstr — each inner cell is a set of predictor names
%                       to combine via element-wise AND.
%                       Default: {} (no interactions).
%                       Example: {{'condition','bodysite'}, {'session_number','condition'}}
%
%   'three_way'         cell of cellstr — same as 'interactions' but length 3.
%                       Example: {{'subject_id','condition','bodysite'}}
%
%   'subject_var'       string. Metadata column identifying subjects.
%                       Default: 'subject_id'.
%
%   'response_transform' 'none' | 'fisherz' | 'rank'.
%                       Default: 'fisherz' (matches the corpus).
%
%   'fit_method'        'REML' | 'ML'.
%                       Default: 'REML'.
%
%   'standardize'       logical. Standardize predictors before fit.
%                       Default: false.
%
%   'use_parallel'      logical. Build the per-subject blocks in parallel (parfor).
%                       Default: false. Per-subject blocks are independent.
%
%   'subsample'         [] | positive integer N | fraction 0<f<1.
%                       Random subsample of within-subject pairs per subject before
%                       LME fit. Useful for very large studies where the full table
%                       is too big for fitlme (>~10^6 rows).
%                       Default: [] (no subsampling).
%
%   'stratify_subsample' cellstr. When 'subsample' is set, stratify the draw by these
%                       columns to preserve marginal balance.
%                       Default: {} (uniform).
%
%   'verbose'           logical. Default: true.
```

### 2.3 — Formula handling

Two forms accepted:

```matlab
% Form 1 — explicit Wilkinson formula (most flexible)
mdl = rsa_lme(dat, 'Y ~ Session + Run + Condition + Bodysite + (1 | Subject)');

% Form 2 — assembled from name-value args (no manual formula string)
mdl = rsa_lme(dat, ...
    'fixed_effects',  {'session_number','run_number','condition','bodysite'}, ...
    'random_effects', {'(1 | subject_id)'}, ...
    'interactions',   {{'condition','bodysite'}});
% Internally constructs: 'Y ~ Session + Run + Condition + Bodysite + CondxBodysite + (1|Subject)'
```

Form 2 is friendlier; Form 1 is the escape hatch. The internal `rsa_lme_formula_parser` (RSA_tools/) builds Form 1 from Form 2 — exposed as a public utility for users who want to preview the formula before fitting.

### 2.4 — Fixed-effects variant: `rsa_lm`

Same shape, but pools across all (i,j) pairs (including cross-subject) and uses `fitlm`:

```matlab
mdl = rsa_lm(dat, ...
    'predictors',     {'subject_id','session_number','run_number','condition','bodysite'}, ...
    'interactions',   {{'condition','bodysite'}, {'session_number','condition'}});
% Internally: builds the omnibus RSM (corr(dat.dat)), vectorizes its upper triangle,
% builds categorical RDMs for each predictor + interaction, calls fitlm.
% Returns LinearModel.
```

Note `subject_id` is *one of the predictors* here (the "SameSubject" categorical RDM), not the random-effects grouping. This is the workflow's lines 1730–1737 pattern.

### 2.5 — Model comparison: `compare_models`

```matlab
[T, best] = R.compare_models( ...
    {'Y ~ 1 + (1|Subject)';
     'Y ~ Condition + (1|Subject)';
     'Y ~ Condition + Bodysite + (1|Subject)';
     'Y ~ Condition*Bodysite + (1|Subject)';
     'Y ~ Condition*Bodysite + Session + (1|Subject)';
     'Y ~ Condition*Bodysite + Session + (Condition|Subject)';
     'Y ~ Condition*Bodysite + Session*(Condition+Bodysite) + (Condition|Subject)';
    }, ...
    'criteria', {'aic','bic','loglik','lrt'}, ...
    'lrt_pairs', 'sequential');     % 'sequential' compares each to previous; 'vs_first' compares all to first

% T columns: model | n_params | aic | bic | loglik | lrt_chi2 | lrt_df | lrt_p
% best: index of model with lowest AIC (or BIC; configurable via 'select_by')
```

`compare_models` uses MATLAB's `compare(LME1, LME2)` under the hood for nested LRTs.

### 2.6 — Per-subject LM fits: `rsa_lm_by_subject`

For when you want the per-subject fixed-effects fits and don't want to bother with LME:

```matlab
% Fits the same fitlm separately per subject, returns aggregated stats
T = rsa_lm_by_subject(dat, ...
    'predictors',  {'session_number','run_number','condition','bodysite'}, ...
    'interactions', {{'condition','bodysite'}});

% T: subject_id | term | beta | se | t | p | partial_R2 | full_R2
% This is the per-subject loop in 08072024 lines 1879–1900.
```

Group-level inference on this output: paired t-test of betas across subjects (uses `rsa_group_inference` from Phase 2).

### 2.7 — Convenience: ICC and BLUP extraction

The workflow extracts subject-specific ICCs and BLUPs from the fitted LME (lines 2184–2280). Expose as a method on the returned `LinearMixedModel`:

```matlab
% Not a new method — just a free helper since we can't extend LinearMixedModel:
icc = rsa_lme_icc(mdl);
%   icc.intercept_subject     — variance(Subject random intercept) / total variance
%   icc.slope_<term>_subject  — same for random slopes
%   icc.per_subject_table     — BLUP-based per-subject ICC

blups = rsa_lme_blups(mdl);
%   table: subject | term | estimate | se
```

---

## 3. Internal contract — `assemble_lme_table`

The single helper that produces the table. Concrete signature:

```matlab
function [tbl, info] = assemble_lme_table(rsm_3d, metadata_per_image, subject_var, ...
                                          predictors, interactions, three_way, ...
                                          response_transform)
% Inputs
%   rsm_3d              [k_s × k_s × N_subjects] cell array — one image-level RSM per subject
%                       (k_s = number of images for that subject; can vary across subjects)
%                       OR: [N_total × N_total] omnibus RSM with subject_id in metadata
%                       (auto-detected)
%   metadata_per_image  table with N_total rows (one per image); columns include
%                       subject_var + all predictors used
%   subject_var         string — name of the subject column
%   predictors          cellstr — main-effect predictor names
%   interactions        cell of cellstr — pairwise interactions to construct
%   three_way           cell of cellstr — three-way interactions to construct
%   response_transform  'none' | 'fisherz' | 'rank'
%
% Outputs
%   tbl                 LME-ready table. Columns: Y, Subject, <predictor 1>, …,
%                       <interaction 1>, …
%   info                struct with bookkeeping:
%                         .n_pairs_per_subject  [N_subjects x 1]
%                         .predictor_provenance struct(predictor -> metadata column)
%                         .formula_skeleton     suggested default formula
```

Per-subject block construction:

```matlab
for s = 1:N_subjects
    is_s = metadata_per_image.(subject_var) == subjects(s);
    meta_s = metadata_per_image(is_s, :);
    R_s    = rsm_3d{s};   % k_s x k_s

    y = upper_triangle(R_s);                          % [n_pairs_s × 1]
    if strcmp(response_transform, 'fisherz'), y = atanh(y); end

    block = table;
    block.Y = y;
    block.Subject = repmat(subjects(s), n_pairs_s, 1);

    % Main effects
    for p = 1:numel(predictors)
        rdm = categorical_rdm(meta_s.(predictors{p}));   % k_s x k_s binary
        block.(predictor_short_name(predictors{p})) = upper_triangle(rdm);
    end

    % Pairwise interactions = elementwise AND
    for k = 1:numel(interactions)
        terms = interactions{k};
        v = ones(n_pairs_s, 1);
        for t = 1:numel(terms), v = v .* block.(predictor_short_name(terms{t})); end
        block.(interaction_name(terms)) = v;
    end

    % Three-way interactions: same logic over 3 terms
    ...

    tbl = [tbl; block];   %#ok<AGROW>
end

% Subject must be categorical for fitlme to treat it as a grouping variable
tbl.Subject = categorical(tbl.Subject);
```

`categorical_rdm(v)` reproduces `rsa.rdm.categoricalRDM` semantics: returns binary `k × k` matrix where entry (i,j) is 1 iff `v(i) == v(j)` and i≠j. Implemented as `double(v(:) == v(:)')` with zeroed diagonal — no rsatoolbox dependency.

`upper_triangle(M)` returns `M(triu(true(size(M)), 1))` — the workflow's idiom verbatim.

`predictor_short_name` strips conventional suffixes: `session_number → Session`, `subject_id → Subject` (rejected as conflict in LME mode; auto-renamed in lm mode), `condition → Condition`, etc. Override via `'short_names'` name-value pair.

`interaction_name({'condition','bodysite'}) → 'CondxBodysite'` to match workflow naming. Configurable.

---

## 4. Nested-model spec format

Three input shapes for `compare_models`:

### 4.1 — Cellstr of Wilkinson formulas (verbose, explicit)

```matlab
{'Y ~ 1 + (1|Subject)';
 'Y ~ Condition + (1|Subject)';
 'Y ~ Condition + Bodysite + (1|Subject)'}
```

### 4.2 — Cell of structs (programmatic)

```matlab
{ struct('fixed', {{}}, 'random', {{'(1|Subject)'}});
  struct('fixed', {{'condition'}}, 'random', {{'(1|Subject)'}});
  struct('fixed', {{'condition','bodysite'}}, 'random', {{'(1|Subject)'}});
  struct('fixed', {{'condition','bodysite','condition:bodysite'}}, ...
         'random', {{'(condition|Subject)'}});
}
```

Each struct is converted to a Wilkinson string by the parser. Use `:` for interactions, `*` for main + interaction.

### 4.3 — Build sequence: `add_term`

```matlab
seq = rsa_model_sequence('Y ~ 1 + (1|Subject)');
seq = seq.add_term('Condition');
seq = seq.add_term('Bodysite');
seq = seq.add_term('Condition:Bodysite');
seq = seq.add_term('Session');
seq = seq.add_term('Session:Condition');
seq = seq.add_term('Session:Bodysite');
seq = seq.upgrade_random({'(Condition|Subject)'});

formulas = seq.formulas;    % cellstr ready for compare_models
[T, best] = R.compare_models(formulas);
```

This is sugar around 4.1; users who don't want to write 10 formulas by hand can build a sequence imperatively. Mirrors what the workflow does in spirit but makes the additivity explicit.

---

## 5. Validation strategy

Every Phase 3 method validated against the `08072024` workflow on a fixed reference dataset (extract once, save to `tests/fixtures/lme_reference.mat`):

### 5.1 — Numerical reproduction

| Workflow line | Test | Acceptance |
|---|---|---|
| Lines 1804–1832 (per-subject categorical RDMs) | `assemble_lme_table` on reference data | Each main-effect column matches `binRDM_*_sub{i}` vectorized to bit |
| Lines 1880–1884 (per-subject regressor stack) | `assemble_lme_table` interaction columns | `CondxBodysite` etc. match to bit |
| Lines 1886–1900 (per-subject fitlm fits) | `rsa_lm_by_subject` | β, SE, p match within MATLAB's default tolerances |
| Lines 1944–1957 (LME table assembly) | `assemble_lme_table` table | Same rows, same columns, same values |
| Line 2129 (basic LME REML fit) | `rsa_lme(R, 'Y ~ Session + Run + Condition + Bodysite + (1|Subject)')` | LME coefficients, log-likelihood, REML criterion match |
| Lines 2226–2249 (BLUPs and subject-specific ICCs) | `rsa_lme_blups`, `rsa_lme_icc` | Match to MATLAB's default tolerances |

### 5.2 — Subsampling correctness

For a small synthetic dataset:
- Fit LME on full table → record coefficients.
- Fit LME on 50% stratified subsample → coefficients should be within 1 SE of full-data coefficients on average over 100 random seeds.

### 5.3 — Soft-dep correctness

Run all tests twice: once with rsatoolbox on path (uses `rsa.rdm.categoricalRDM`), once with it removed (uses our `fallback_categoricalRDM`). Outputs must agree to bit.

### 5.4 — MATLAB MCP live verification

The `mcp__MATLAB__run_matlab_file` tool lets us execute the reference workflow's exact LME calls and our `rsa_lme` calls side-by-side and diff their `disp(mdl)` outputs. Use this for the Phase 3 acceptance gate.

---

## 6. Resolved decisions  *(closed 2026-05-29)*

All six open questions resolved. Lock these in for Phase 3 implementation; no further sign-off needed before code starts.

1. **Per-subject row counts vary** → **No weighting by default; expose `'weight_by_subject', true`.** Matches `fitlme`'s native behavior and the `08072024` workflow. LME's random effects already absorb subject-level variance; explicit weighting is an opt-in.

2. **`Subject` conflict in `rsa_lme`** → **Error with clear message.** If `subject_id` (or whatever `subject_var` is set to) appears in `predictors` for `rsa_lme`, throw:
   > `subject_id is the random-effects grouping in rsa_lme; remove it from 'predictors' or use rsa_lm if you want Subject as a fixed effect.`

   No silent coercion. Same-vs-different "SameSubject" is reachable only via `rsa_lm`.

3. **Default `response_transform`** → **`'fisherz'`.** LME residuals are unbounded under Fisher-z; bounded under raw correlations. Users who want raw correlations (to match the `08072024` numbers exactly) pass `'response_transform','none'`.

4. **Continuous-distance predictors** → **Add `rsm.from_metadata_distance` in Phase 1.** Single distance metric (`'abs_diff'`) to start; the broader `{abs_diff, squared_diff, log_diff, custom_fcn}` set is a future extension. Lets users model session-distance, run-distance, time-since-baseline, etc. as continuous predictors.
   ```matlab
   M = rsm.from_metadata_distance(meta, 'session_number', 'metric', 'abs_diff');
   % Then in rsa_lme:
   mdl = R.rsa_lme('predictors', {'condition','bodysite'}, ...
                   'numeric_predictors', {M}, ...
                   'random_effects', {'(1|Subject)'});
   ```

5. **Parcelwise LME output granularity** → **Long table + `statistic_image` per term.**
   ```matlab
   results.lme_table   % long: parcel | term | beta | se | t | p | FDR_P | sig
   results.maps.<term> % statistic_image (beta or t per parcel — configurable)
   ```
   Per-parcel `LinearMixedModel` objects are not retained by default (memory); enable with `'keep_models', true` for deep inspection.

6. **Cross-subject pairs** → **Top-level `rsa_lm` uses all pairs (workflow default); `rsa_parcelwise(dat, 'lm', ...)` uses within-subject only (matches LME).** Both expose `'pair_scope', 'within_subject' | 'all'` for override.

---

## 7. Status

**Signed off 2026-05-29.** All §6 questions resolved with recommended defaults; design is complete. Phase 3 implementation can proceed independently of (and in parallel with) Phase 1 — `assemble_lme_table` only needs an `fmri_data` object plus the metadata table; it doesn't depend on the rest of `@rsm`.
