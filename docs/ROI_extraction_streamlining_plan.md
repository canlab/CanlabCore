# ROI / Image Data Extraction — Streamlining & Overhaul Plan

Status: PROPOSAL (draft)
Scope: All data-extraction code in CanlabCore — object methods (`@fmri_data`, `@image_vector`, `@region`, `@atlas`) and the standalone `Data_extraction/` (and `peak_coordinates/`) functions.
Audience: CanlabCore maintainers.

This document is grounded in a full read of the actual `.m` sources and a repo-wide call-graph/usage analysis (`grep -rn --include='*.m'`). Function signatures and call sites cited below were extracted from source, not guessed.

---

## 1. Problem statement: redundancy and drift among extraction methods

The toolbox currently exposes **at least 15 distinct functions** that all do some variant of "pull voxel values out of images, restrict to a mask/atlas/region, and (optionally) reduce to per-region means or pattern expression." They have accreted over ~20 years across three eras: (a) the legacy SPM **`clusters`-struct** era, (b) the standalone procedural era, and (c) the current **object-oriented** era. The result is redundancy, silent behavioral drift, and dead code.

Concrete symptoms found in the source:

- **Two near-identical standalone twins.** `Data_extraction/extract_image_data.m` and `Data_extraction/extract_from_rois.m` share essentially the same body and H1 text. They differ only in default averaging mode (`'contiguous_regions'` vs `'unique_mask_values'`) and output argument order. **Neither has any live caller** — the only `extract_image_data(...)` call in the repo resolves to a *local subfunction* inside `Parcellation_tools/parcel_images.m`, not the standalone.

- **Two overlapping OO `extract_roi_averages` implementations that have drifted.** `@fmri_data/extract_roi_averages.m` (`[cl, cl_roisum, cl_demeanedpattern] = extract_roi_averages(obj, mask_image, varargin)`) is the rich version with pattern expression, `statistic_image`/`fmri_mask_image` special-casing, and atlas-aware resampling. `@image_vector/extract_roi_averages.m` (`cl = extract_roi_averages(obj, mask, varargin)`) is a thinner copy that **explicitly `error`s on `'pattern_expression'`** and whose actual default (`'unique_mask_values'`) disagrees with its own help text (which says contiguous). Its own header admits: *"Better to have only one function of record in the future."*

- **Inconsistent resampling direction across methods doing the same job.** `extract_roi_averages` (both), `apply_mask`, and `apply_parcellation` resample the **mask/atlas → data** space. `@atlas/extract_data` resamples the **data → atlas** space (then `apply_parcellation` re-matches internally). `@region/extract_data` does **no resampling** at all (mm-coordinate intersection). Same conceptual operation, three different space contracts — a correctness and reproducibility hazard.

- **A documented-but-unimplemented method.** `apply_atlas` is named in the help headers of `@image_vector/image_vector.m` (lines 23, 143) and `@fmri_data/fmri_data.m` (lines 23, 141) as "Computes the mean value or pattern expression for each reference region specified in an atlas object," but **no `apply_atlas.m` file exists anywhere**. Users following the docs hit a missing method. The real implementations are `@image_vector/apply_parcellation.m` and `@atlas/extract_data.m`.

- **An unresolved merge-conflict file shipped in the repo.** `Data_extraction/sphere_roi_tool_2008_conflict.m` literally contains `<<<<<<< .mine` / `=======` / `>>>>>>> .r1197` markers (lines 1-5, 43-52, 59-65) and a duplicate `function sphere_roi_tool_2008` definition. It cannot parse and is a name-collision hazard.

- **A fully dead `Grandfathered/` folder.** `Data_extraction/Grandfathered/timeseries4.m` and `check_timeseries_vals.m` have no external callers; every live call resolves to local subfunctions inside `cluster_orthviews.m`.

- **A legacy clusters-struct chain with an orphaned head.** `extract_raw_data` → `tor_extract_rois` / `extract_indiv_peak_data`. `extract_raw_data` has no live external caller, and `extract_indiv_peak_data` is reachable *only* through it.

- **Doc-quality debt.** `poorly_documented_functions.txt` flags `apply_mask` (line 27, 0/3), `check_timeseries_vals` (200), `timeseries4` (201), `sphere_roi_tool_2008` + its conflict twin (204-205), and `extract_raw_data` (521).

The cost: contributors cannot tell which function is canonical; bug fixes land in one copy and not its twin (e.g. the `@fmri_data` version got the `'nearest'`-interp fix for integer atlas codes and the `atlas2region` fix; the `@image_vector` version did not); and users are pointed (by help text) at dead standalones and a non-existent `apply_atlas`.

---

## 2. Target architecture: one core engine, thin wrappers

### 2.1 Principle

There should be exactly **one computational core** for "restrict-to-mask + reduce." Every public method becomes a thin, well-documented wrapper that normalizes its inputs to that core and shapes the outputs to its return type. No two functions should independently re-implement voxel intersection, weight normalization, or pattern expression.

### 2.2 What the core should be

The strongest existing candidate for the core is **`@image_vector/apply_parcellation.m`** plus **`canlab_pattern_similarity`** (the pattern-expression primitive that almost everything already funnels through). Reasons it should be the engine of record:

- It already produces the full superset of useful outputs: `[parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg, voxel_count, parcel_ste]`.
- It handles the integer-label-preserving resample correctly (`resample_space(parcels, dat, 'nearest')`), full-width output via `condf2indic(parcels.dat, 'integers', n_orig_parcels)`, and the `probability_maps = []` speedup for atlases.
- It computes true weighted means via a scaled matrix product and routes pattern expression through `canlab_pattern_similarity`, which is the same primitive `apply_mask`, `@region/extract_data`, and `@fmri_data/extract_roi_averages` already call.
- `@atlas/extract_data` is *already* a thin wrapper over it. `@fmri_data/extract_measures_batch` already calls it directly.

Recommended concrete structure:

```
canlab_extract_core(data_obj, mask_or_atlas, 'reduce', <mean|pattern|both>, ...)   % NEW or = refactored apply_parcellation
        |
        |-- handles: input normalization (char/struct/stat_image/atlas -> canonical),
        |            ONE documented space contract (resample mask/atlas -> data, 'nearest' for integer labels),
        |            voxel intersection, weight L1/L2 normalization, canlab_pattern_similarity dispatch
        v
  returns the full matrix superset (means, pattern, valence, rmsv, counts, ste)

Wrappers (public, stable names, shape outputs only):
  @fmri_data/extract_roi_averages   -> region-object output (cl, cl_roisum, cl_demeanedpattern)
  @image_vector/extract_roi_averages-> region-object output (mean only) — OR deprecate-merge (see table)
  @atlas/extract_data               -> matrix/table output (already thin)
  @region/extract_data              -> region-object output, KEEPS its mm-coordinate path (distinct contract — see below)
  apply_parcellation                -> public alias / direct core entry (keep name; it IS the core)
  apply_mask                        -> KEEP as the global (non-region-aware) masking primitive
```

### 2.3 Deliberate exception: `@region/extract_data`

`@region/extract_data` must **not** be naively folded into the resampling core. Its mm-coordinate intersection (`mm2voxel` + `intersect(...,'rows')`, no interpolation) is a *different and intentional* contract — its header explicitly notes results differ from a resample-then-extract approach. It should keep its coordinate path but **share the reduction step** (mean + `canlab_pattern_similarity` weighting) with the core via a small shared helper, so the "what does pattern expression mean" logic lives in one place.

### 2.4 Space-contract unification

Pick **one** documented default: *resample the mask/atlas/parcels to the data space, using `'nearest'` for integer-labeled parcellations and the atlas-aware path for atlases.* `@atlas/extract_data`'s data→atlas direction should be retained only as an explicit, documented performance option (`'resample','data_to_atlas'`) for large-N, not as a silent divergence. Every wrapper documents which contract it uses and why.

### 2.5 Fix the docs-vs-reality gaps

- Either implement `apply_atlas` as a one-line alias to `apply_parcellation`/`extract_data`, **or** remove it from the class help headers. (Recommend: add a thin `apply_atlas` alias so the documented name resolves — lowest-surprise for users.)
- Correct the `@image_vector/extract_roi_averages` help/default mismatch.
- Repoint the "non-object-oriented alternative" help in `@fmri_data/extract_roi_averages.m:96` and `@atlas/extract_data.m:114` away from the dead `extract_image_data` standalone.

---

## 3. Disposition table

Verdicts: **KEEP-CORE** (canonical, becomes/anchors the engine) · **KEEP** (distinct, still needed) · **CONSOLIDATE-INTO-X** (merge behavior into X, leave thin wrapper) · **DEPRECATE** (keep callable, emit warning, schedule removal) · **REMOVE** (delete now or end of plan).

| Function | Path (relative to CanlabCore/) | Live callers | Disposition | Rationale & risk |
|---|---|---|---|---|
| `@image_vector/apply_parcellation` | `@image_vector/apply_parcellation.m` | 5 (atlas/extract_data, fmri_data/extract_measures_batch, image_vector/wedge_plot_by_atlas) | **KEEP-CORE** | Becomes the consolidation target / engine. Produces the full output superset; already the computational backbone. Risk: low; central, so any refactor must be test-guarded hard. |
| `@image_vector/apply_mask` | `@image_vector/apply_mask.m` | 57 (13 distinct @class methods incl. `@statistic_image/threshold`) | **KEEP-CORE** | Most-used masking primitive in the toolbox; global (one-mask) masking, not region-aware. Keep as the masking layer the core uses. Doc debt: flagged 0/3 (poorly_documented_functions.txt:27) — improve header. Risk: very high blast radius; do NOT change semantics, only document. |
| `@fmri_data/extract_roi_averages` | `@fmri_data/extract_roi_averages.m` | ~5 (region/region.m, fmri_data/canlab_connectivity_preproc) | **KEEP-CORE (wrapper)** | The primary OO ROI extractor and the version with the bug fixes (nearest-interp for unique values, atlas2region). Re-home its compute onto the core; keep its region-object output contract and the `cl_roisum`/`cl_demeanedpattern` outputs. Risk: medium — output struct/region shape must be byte-stable for callers. |
| `@image_vector/extract_roi_averages` | `@image_vector/extract_roi_averages.m` | tests + inheritance | **CONSOLIDATE-INTO @fmri_data version / core** | Thinner drifted copy; errors on pattern expression; help/default mismatch; its own header asks for a single function of record. Make it delegate to the shared core (mean-only path) so the two cannot drift again. Risk: medium — it is the inherited superclass method; subclasses other than fmri_data rely on it. Must verify class dispatch unchanged. |
| `@atlas/extract_data` | `@atlas/extract_data.m` | 1 (region/ttest_table_by_condition) + public | **KEEP (thin wrapper)** | Already delegates to `apply_parcellation`. Keep; align its space contract documentation (data→atlas) as an explicit option. Risk: low. |
| `@region/extract_data` | `@region/extract_data.m` | 1 (region/ttest_table_by_condition) + public | **KEEP (distinct contract)** | mm-coordinate, no-resample path is intentional and unique. Share only the reduction helper with the core. Risk: low if coordinate path untouched. |
| `@image_vector/extract_gray_white_csf` | `@image_vector/extract_gray_white_csf.m` | callers via batch | **KEEP** | Purpose-built tissue-compartment summaries; thin over `apply_mask`. No redundancy. Risk: none. |
| `@fmri_data/extract_measures_batch` | `@fmri_data/extract_measures_batch.m` | orchestrator | **KEEP** | Orchestration only; already calls the core (`apply_parcellation`) directly. Risk: none. |
| `@region/check_extracted_data` | `@region/check_extracted_data.m` | QC | **KEEP** | Standalone QC sanity check. Minor latent bug (`isok` overwritten each loop) — fix opportunistically. Risk: none. |
| `extract_image_data` | `Data_extraction/extract_image_data.m` | **0 live** | **DEPRECATE → REMOVE** | Orphaned standalone twin of `extract_from_rois`; superseded by OO methods. The OO help text falsely advertises it. Deprecate with a pointer to `extract_roi_averages`, remove in Phase 3. Risk: low — but external/downstream scripts may call it; hence deprecate-first, not delete-now. |
| `extract_from_rois` | `Data_extraction/extract_from_rois.m` | **0 live** | **DEPRECATE → REMOVE** | Older duplicate of the above (different default + arg order). Same path. Risk: low; deprecate-first for downstream safety. |
| `extract_raw_data` | `Data_extraction/extract_raw_data.m` | 0 live external | **DEPRECATE** (KEEP-LEGACY-WRAPPER) | Head of a legacy EXPT/clusters chain; no live external caller but it drives `tor_extract_rois` + `extract_indiv_peak_data`. Deprecate; do not remove until the chain is retired. Doc-flagged (521). Risk: medium — old user scripts/EXPT pipelines. |
| `extract_indiv_peak_data` | `Data_extraction/extract_indiv_peak_data.m` | 1 (only `extract_raw_data`) | **DEPRECATE → REMOVE (with its parent)** | Dead once `extract_raw_data` goes. Remove together. Risk: low. |
| `tor_extract_rois` | `Data_extraction/tor_extract_rois.m` | 34 across ~14 legacy files | **KEEP-LEGACY-WRAPPER** | High usage but only from legacy procedural/cluster tooling (`mask2clusters`, `cluster_orthviews`, `mask_*`, `parcel_*`, `hewma_*`, `classify_naive_bayes*`); never from an @class method. Long pole — keep until those tools migrate. Note local-subfunction shadow in `cluster_orthviews.m:503`. Risk: high to remove; safe to keep. |
| `extract_contrast_data` | `Data_extraction/extract_contrast_data.m` | 2 (`cluster_tool_getbetas`) | **KEEP-LEGACY-WRAPPER** | Used only by the GUI `cluster_tool`; chains to `tor_extract_rois` + `extract_ind_peak`. Keep until GUI migrates. Risk: medium. |
| `extract_ind_peak` | `peak_coordinates/extract_ind_peak.m` | 2 (`cluster_barplot`, `extract_contrast_data`) | **KEEP-LEGACY-WRAPPER** | Legacy individual-peak helper. Name-collision hazard with `extract_indiv_peak_data`. Keep until its two callers retire. Risk: medium. |
| `canlab_maskstats` | `Data_extraction/canlab_maskstats.m` | 2 (`canlab_glm_maskstats`) | **KEEP-LEGACY-WRAPPER** | Overlaps `apply_mask 'pattern_expression'` + `extract_roi_averages` but still wired into the GLM_Batch_tools pipeline. Keep until that pipeline is refactored, then DEPRECATE. Risk: medium. |
| `canlab_load_ROI` | `Data_extraction/canlab_load_ROI.m` | active | **KEEP** | Current, object-aware (returns `region`/`atlas`). Not an extractor of data — a region loader. No change. Risk: none. |
| `canlab_extract_ventricle_wm_timeseries` | `Data_processing_tools/canlab_extract_ventricle_wm_timeseries.m` | active | **KEEP** | Current; aCompCor-style nuisance extraction with PCA comps; not replaced by `extract_gray_white_csf`. Risk: none. |
| `timeseries_extract_slice` | `Data_extraction/timeseries_extract_slice.m` | active | **KEEP** | Low-level single-slice I/O; no object equivalent. Risk: none. |
| `sphere_roi_tool_2008` | `Data_extraction/sphere_roi_tool_2008.m` | 4 (ROI-building scripts) | **KEEP** | Active, clean, has 2023/2025 examples bridging to `region` via `cluster2region`. Risk: none. |
| `sphere_roi_tool_2008_conflict.m` | `Data_extraction/sphere_roi_tool_2008_conflict.m` | 0 | **REMOVE NOW** | Unresolved merge-conflict file; will not parse; duplicate function name. No callers. Zero risk. Doc-flagged (205). |
| `Grandfathered/timeseries4.m` | `Data_extraction/Grandfathered/timeseries4.m` | 0 (live calls hit local subfn in `cluster_orthviews.m`) | **REMOVE (end of Phase 1)** | Fully dead in `Grandfathered/`. Doc-flagged (201). Risk: very low — verify the `cluster_orthviews.m:765` local subfunction shadow is the actual resolver before deleting. |
| `Grandfathered/check_timeseries_vals.m` | `Data_extraction/Grandfathered/check_timeseries_vals.m` | 0 (used only by sibling Grandfathered file) | **REMOVE (with timeseries4)** | Dead validation helper. Doc-flagged (200). Risk: very low. |
| `apply_atlas` (documented, no file) | (none) | — | **CREATE alias OR REMOVE from docs** | No file exists; named only in class help headers. Add a thin alias to `apply_parcellation`, or strike from help. Risk: low; user-facing doc correctness. |

---

## 4. Phased migration plan

Each phase is gated by `canlab_run_all_tests` passing. Tests are function-based (`functiontests(localfunctions)`), discovered by the custom glob in `Unit_tests/canlab_run_all_tests.m` (files must be named `canlab_test_*.m`). Fixtures come from `Unit_tests/helpers/canlab_get_sample_fmri_data.m` (emotionreg, 30 images) and `canlab_get_sample_thresholded_t.m`.

### Phase 0 — Zero-risk hygiene (immediate)
1. `git rm Data_extraction/sphere_roi_tool_2008_conflict.m` (broken, 0 callers).
2. Confirm via `grep -rn 'timeseries4\|check_timeseries_vals'` that the only resolvers for the live calls are the local subfunctions in `cluster_orthviews.m`; then `git rm` the two `Grandfathered/` files (or, if any doubt remains, defer their deletion to end of Phase 1 after the characterization tests below exist).
3. Commit. Run `canlab_run_all_tests`; expect green (these files have no callers, so nothing should change).

### Phase 1 — Characterization tests + deprecation warnings (no behavior change)
Goal: lock in current behavior of the KEEP-CORE functions *before* refactoring, and start warning on the doomed ones.

Steps:
1. **Add value-correctness tests** (extend, do not duplicate, the existing `Unit_tests/image_vector/canlab_test_extract_roi.m`). The existing file already covers `(data, mask_filename)` return type and `.dat`-row-count==n-images. Add new local test functions covering the **currently uncovered** ground that the refactor will touch:
   - Synthetic known-value test: build an `fmri_data` whose `.dat` is a known constant (e.g. all 7s) over `brainmask_canlab.nii`, run `apply_parcellation` against `load_atlas('canlab2024')`, and `verifyEqual` that every parcel mean == 7 (within tolerance). Pin the analytic mean so the core refactor is guarded.
   - Atlas-input / `'unique_mask_values'` multi-region path (the case the old `old_to_integrate/check_roi_extraction.m` exercised) — verify number of regions and per-region means.
   - On-the-fly resampling: extract with a mask in a different space; verify it runs warning-free and returns sane shapes.
   - Pattern-expression equivalence: verify `@fmri_data/extract_roi_averages(..., 'pattern_expression')` and the corresponding `apply_parcellation(..., 'pattern_expression', w)` agree within tolerance on a shared ROI — this is the cross-check that lets us merge the two later.
   - Cross-method space-contract test: document/verify current resampling direction of each method so changes are visible.
   - Assertion style per house convention: `tc.verifyEqual`, `tc.verifyClass`, `tc.assumeNotEmpty(which('brainmask_canlab.nii'), ...)` to skip cleanly when data is off-path.
2. **Add deprecation warnings** (callable, no behavior change) to: `extract_image_data`, `extract_from_rois`, `extract_raw_data`, `extract_indiv_peak_data`. Use a standard one-liner at the top of each:
   ```matlab
   warning('CanlabCore:deprecated', ...
       '%s is deprecated and will be removed. Use @fmri_data/extract_roi_averages or @atlas/extract_data instead.', mfilename);
   ```
   Warnings must NOT change return values (downstream scripts keep working).
3. **Fix docs-vs-reality now** (cheap, high value): correct the `@image_vector/extract_roi_averages` help/default mismatch; repoint the "non-object-oriented alternative" pointers (`@fmri_data/extract_roi_averages.m:96`, `@atlas/extract_data.m:114`); add the `apply_atlas` alias (or strike it from the class headers).
4. Run `canlab_run_all_tests`. Green gate. The new characterization tests now define "correct."

How tests guard: the synthetic known-value and pattern-expression-equivalence tests pin the exact numerical contract of the core *before* any code moves, so Phase 2 cannot silently alter results.

### Phase 2 — Consolidate onto the core (behavior-preserving refactor)
Goal: make every wrapper delegate to one engine; eliminate drift.

Steps:
1. Factor the shared reduction logic (voxel intersection + weight L1/L2 normalization + `canlab_pattern_similarity` dispatch + the full output superset) into the core (`apply_parcellation`, possibly fronted by a `canlab_extract_core` helper). Keep `apply_parcellation`'s public signature unchanged.
2. Re-point `@fmri_data/extract_roi_averages` to compute via the core, then shape into `region` objects (`cl`, `cl_roisum`, `cl_demeanedpattern`). Preserve the exact output struct/region fields.
3. Make `@image_vector/extract_roi_averages` delegate to the same core (mean-only path). Remove the divergent inline resampling. This kills the fmri_data-vs-image_vector drift.
4. Unify the space contract: single documented default (mask/atlas → data, `'nearest'` for integer labels), with `@atlas/extract_data`'s data→atlas as an explicit documented `'resample'` option.
5. Share the reduction helper into `@region/extract_data` while preserving its mm-coordinate intersection path.
6. After **each** step, run `canlab_run_all_tests`. The Phase 1 tests (value-correctness, pattern-expression equivalence, shape, cross-method agreement) are the guardrail: any numeric or shape change fails the gate.

How tests guard: because Phase 1 pinned analytic means and cross-method equivalence, a refactor that changes resampling, normalization, or output shape fails immediately. The pattern-expression-equivalence test specifically guards the riskiest merge (fmri_data ↔ image_vector ↔ apply_parcellation).

### Phase 3 — Remove deprecated code
Goal: delete what Phase 1 warned about, once a deprecation window (recommend: one tagged release / ~6 months) has elapsed.

Steps:
1. Remove `extract_image_data.m` and `extract_from_rois.m` (0 callers, deprecated since Phase 1).
2. Remove `extract_raw_data.m` together with `extract_indiv_peak_data.m` (its only caller) — *only if* the EXPT/clusters legacy pipeline is confirmed retired. Otherwise keep them as KEEP-LEGACY-WRAPPER and re-evaluate next cycle.
3. Remove the `Grandfathered/` directory if not already done in Phase 0.
4. Re-run the full suite + a spot-check of `CANlab_help_examples` walkthroughs (see Risks).
5. For functions that must keep their public name for backward compat but whose body is gone, leave a **thin wrapper** that calls the core and emits the deprecation warning (see §5).

How tests guard: removal of zero-caller functions cannot break the suite by definition; the walkthrough spot-check (Phase 3 step 4) catches example scripts that the unit suite does not exercise.

---

## 5. Backward-compatibility strategy

1. **Deprecate, never silently delete, anything with a non-zero (or uncertain) external surface.** Functions named in user-facing help, in `CANlab_help_examples`, or plausibly called by downstream repos get a deprecation warning first and a removal only after a full release cycle.
2. **Keep function names as thin wrappers.** When `extract_image_data` / `extract_from_rois` behavior is fully subsumed by the OO methods, the *name* can survive as a one-screen wrapper that (a) emits `warning('CanlabCore:deprecated', ...)` once, and (b) forwards to `extract_roi_averages`/`apply_parcellation`, mapping args and return order. This keeps old scripts running while steering new code to the core.
3. **Stable public signatures.** `apply_parcellation`, `apply_mask`, `@fmri_data/extract_roi_averages`, `@atlas/extract_data`, `@region/extract_data` keep their exact input/output signatures throughout. Internal re-homing onto the core must not change argument names, option strings, or output order/fields.
4. **Warning hygiene.** Use a single warning identifier namespace (`CanlabCore:deprecated`) so users can `warning('off', 'CanlabCore:deprecated')` if needed, and so the unit tests can assert deprecation warnings deliberately (`verifyWarning`) without `verifyWarningFree` checks elsewhere tripping.
5. **`apply_atlas` resolution.** Add the documented-but-missing `apply_atlas` as a real (thin) alias rather than removing the promise from the docs — least surprising for users who copied the help.
6. **Changelog + deprecation table.** Maintain a short "deprecated / removed" list in the docs and release notes so downstream maintainers can grep their code ahead of removal.

---

## 6. Risks and what could break

- **`CANlab_help_examples` walkthroughs and tutorials** are not exercised by the unit suite (the runner skips `walkthroughs/` by default). They are the most likely place to call `extract_roi_averages`/`extract_data`/`apply_parcellation` with subtle expectations about output shape. **Mitigation:** before any Phase 3 removal and after Phase 2 consolidation, run a grep of `CANlab_help_examples` for every touched function name and execute the affected walkthroughs headless.

- **Downstream CANlab repos** (e.g. `canlab_single_trials`, `CanlabPrivate`, project analysis scripts) may call the standalone `extract_image_data`/`extract_from_rois`/`extract_raw_data` even though there are zero callers *inside CanlabCore*. The repo-internal call-graph cannot see them. **Mitigation:** this is exactly why these are DEPRECATE-then-REMOVE, not REMOVE-now; the warning gives downstream a release cycle to migrate.

- **Output-shape sensitivity.** `@fmri_data/extract_roi_averages` returns `region` objects with specific fields (`.dat`, `.all_data`, `.val`, `.val_descrip`) and the secondary `cl_roisum`/`cl_demeanedpattern` outputs populated only under `'pattern_expression'` + `nargout>1`. Callers like `@region/region.m:286-288` and `canlab_connectivity_preproc` depend on this. **Mitigation:** Phase 1 characterization tests pin these fields; refactor must keep them byte-stable.

- **Resampling-direction change is a correctness risk.** Unifying the space contract could shift values for code that implicitly relied on `@atlas/extract_data`'s data→atlas direction (different interpolation than mask→data). **Mitigation:** keep data→atlas as an explicit option; add the cross-method agreement test in Phase 1; document the chosen default loudly.

- **Local-subfunction shadowing.** `tor_extract_rois`, `timeseries4`, `check_timeseries_vals`, and `extract_image_data` each have like-named *local subfunctions* inside `cluster_orthviews.m` / `parcel_images.m` that currently satisfy the live calls. Deleting the standalone is safe **only** if the local subfunction is the true resolver. **Mitigation:** confirm with `which -all` / a path-order check in MATLAB before each deletion (Phase 0 step 2, Phase 3 step 3).

- **Name-collision hazards** between `extract_ind_peak` (peak_coordinates) and `extract_indiv_peak_data` (Data_extraction) — similar purpose, different signatures. Removing one must not change which the other's callers resolve to. **Mitigation:** verify callers explicitly; these are KEEP/DEPRECATE, not REMOVE-now.

- **The custom test runner only discovers `canlab_test_*.m`.** New guard tests must follow that exact naming or they will silently not run. **Mitigation:** add new local functions into the existing `canlab_test_extract_roi.m` (and a new `canlab_test_apply_parcellation.m` if warranted), and confirm they appear in the `canlab_run_all_tests` suite count.

---

### Bottom line
Delete the obviously dead/broken items now (conflict file, `Grandfathered/`), pin current behavior with synthetic value-correctness and cross-method-equivalence tests, converge the OO `extract_*`/`apply_parcellation` methods onto a single space-correct core (`apply_parcellation` + `canlab_pattern_similarity`) while keeping `@region/extract_data`'s mm-coordinate path distinct, and deprecate-then-remove the orphaned standalones over one release cycle behind `CanlabCore:deprecated` warnings.
