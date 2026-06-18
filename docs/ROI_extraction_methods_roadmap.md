# CANlabCore ROI / Atlas / Parcel Data-Extraction Roadmap

**Audience:** CANlab developers/maintainers (Wager lab).
**Scope:** Every method/function in CanlabCore that extracts data (means, pattern responses, timeseries) from images using ROIs, masks, atlases, or parcellations.
**Status:** Internal developer reference. Redundancy/removal notes here feed forward to a separate streamlining plan; this file does not itself delete anything.

> All paths absolute. Signatures quoted verbatim from source. "Resampling direction" is the single most important cross-method gotcha — see §4.

---

## 1. Conceptual model

CanlabCore is an object-oriented neuroimaging toolbox. Data extraction always combines **a data object** with **a region-definition object**.

- **Images = data objects.** `@image_vector` is the superclass; `@fmri_data` (subject/trial maps, timeseries, with `.Y`, `.X`, `metadata_table`), `@statistic_image` (adds `.sig`, `.ste`, `.p`), and `@fmri_mask_image` are subclasses. Internally, data live in `obj.dat` shaped **`[voxels x images]`**, with `obj.volInfo` carrying the spatial mapping (affine `.mat`, `.xyzlist`). The `replace_empty` / `remove_empty` invariant governs which voxels are present — reason about voxel positions accordingly.
- **Masks / ROIs.** A mask can be a char filename (→ built into `fmri_mask_image`), an `fmri_mask_image`, any `image_vector`/`fmri_data` whose `.dat` holds weights, or a `@statistic_image` (thresholded by `.sig > 0`). A **`@region`** object is the modern multi-region container (coordinate-based: `XYZ`, `XYZmm`, `.val` for optional local pattern weights, `.dat`/`.all_data` after extraction).
- **Atlases = `@atlas` objects.** Integer-coded `.dat` (one code per parcel), plus `.probability_maps`, `.labels`, `.label_descriptions`, `.references`. Loaded by keyword via `load_atlas`.
- **Parcellations.** Any integer-coded image (`atlas` preferred over a plain `fmri_data` parcel image, because atlas-aware resampling preserves integer labels). Extraction produces an `[images x parcels]` matrix.

Two extraction paradigms coexist:
1. **Resample-and-index** (object methods): bring mask/atlas and data into a common space, intersect in-mask voxels, average or apply weights.
2. **Coordinate intersection** (`@region/extract_data`): no resampling; match region mm-coordinates directly against the data's voxel list.

---

## 2. Categorized catalog

### (a) Primary object methods — RECOMMENDED

#### `@fmri_data/extract_roi_averages.m`
- **Path:** `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/@fmri_data/extract_roi_averages.m`
- **Purpose:** Region means (or per-region pattern expression) from an `fmri_data` object, returned as a `region` object array.
- **Signature:** `function [cl, cl_roisum, cl_demeanedpattern] = extract_roi_averages(obj, mask_image, varargin)`
- **Inputs:** `obj` (fmri_data); `mask_image` = char / `fmri_mask_image` / `image_vector` / `fmri_data` / `atlas` / `statistic_image` (empty → `obj.mask`). Options: `'unique_mask_values'` (DEFAULT — average over integer codes; atlas-style), `'contiguous_regions'` (average over contiguous blobs), `'pattern_expression'` (weighted similarity using mask values as weights), `'nonorm'` (skip L1 weight norm), `'cosine_similarity'`, `'correlation'`, `'noverbose'`/`'verbose'`, `'notables'`.
- **Outputs:** `cl` region array (`.dat` = mean or pattern expr; `.all_data` = `[images x voxels]`; `.val` = weights). `cl_roisum`, `cl_demeanedpattern` only populated under `'pattern_expression'` + `nargout>1`.
- **When to use:** The default ROI-mean / per-region pattern extractor for subject-level fmri_data. Richest version (statistic_image special-casing, pattern expression).
- **Calls:** `region(mask, average_over)` or `atlas2region(...)`; `canlab_pattern_similarity` for pattern expression; inline `resample_space`.
- **Caveat:** Header warns it can LOSE removed image data; mask is reloaded from disk if remapped (manual thresholding lost).

#### `@image_vector/extract_roi_averages.m`
- **Path:** `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/@image_vector/extract_roi_averages.m`
- **Purpose:** Simpler region-mean extractor inherited by all `image_vector` subclasses (fmri_data overrides it with the richer version above).
- **Signature:** `function cl = extract_roi_averages(obj, mask, varargin)`
- **Inputs:** `obj` (image_vector); `mask` (image_vector or filename). Options: `'unique_mask_values'` (DEFAULT), `'contiguous_regions'`, `'noverbose'`. `'pattern_expression'` → explicit `error` (only defined for fmri_data).
- **Outputs:** `cl` region array; `.dat` = per-image region mean; `.all_data` = `[images x voxels]`.
- **When to use:** Generic image_vector objects where the fmri_data override doesn't apply. Mean only; no weighting.
- **Calls:** `region()` / `atlas2region()`; inline `resample_space`.

#### `@atlas/extract_data.m`
- **Path:** `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/@atlas/extract_data.m`
- **Purpose:** Extract parcel means / per-parcel pattern responses for every region of an atlas.
- **Signature:** `function [parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg, voxel_counts, parcel_stes] = extract_data(obj, data_obj, varargin)`
- **Inputs:** `obj` (atlas); `data_obj` (fmri_data). Options forwarded to `apply_parcellation`: `'pattern_expression'` + weight_map, `'correlation'`, `'cosine_similarity'`, `'norm_mask'`, `'ignore_missing'`.
- **Outputs:** Same set as `apply_parcellation` (means `[images x parcels]`, pattern expression, valence, signed RMS pos/neg, voxel counts, per-parcel SE). Header usage line wrongly says `data_table`; it returns matrices.
- **When to use:** Whole-atlas parcel extraction with the atlas as the primary object. Thin wrapper.
- **Calls:** `apply_parcellation(data_obj, obj, varargin{:})` (delegates the whole computation).
- **Caveat:** Resamples **DATA → ATLAS** space (opposite of the others). Pre-resampling the atlas to functional space recommended for large N.

#### `@region/extract_data.m`
- **Path:** `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/@region/extract_data.m`
- **Purpose:** Attach voxelwise + averaged data (and optional local pattern response) to a region object array, by coordinate matching.
- **Signature:** `function [r, local_pattern_response] = extract_data(r, data_obj, varargin)`
- **Inputs:** `r` (region array with `XYZmm`; `.val` optionally local weights); `data_obj` (image_vector/fmri_data). Options: `'cosine_similarity'`, `'correlation'`.
- **Outputs:** `r` with `.dat` (region mean per image), `.all_data` (`[images x voxels]`, NaN where unmatched), `.source_images`. `local_pattern_response` only if `nargout>1` and `.val` nonempty.
- **When to use:** You already have `region` objects (e.g. from thresholded results) and want data without interpolation artifacts.
- **Calls:** `mm2voxel`; local `match_coordinates` (`intersect(...,'rows')`); `canlab_pattern_similarity`.
- **Caveat:** Coordinate-based — NO resampling. Slightly different from resample-first approaches by design.

#### `@image_vector/apply_parcellation.m`
- **Path:** `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/@image_vector/apply_parcellation.m`
- **Purpose:** Computational engine for parcel means + per-parcel pattern expression + valence + signed RMS.
- **Signature:** `function [parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg, voxel_count, parcel_ste] = apply_parcellation(dat, parcels, varargin)`
- **Inputs:** `dat` (image_vector/fmri_data); `parcels` (atlas preferred, or fmri_data). Options: `'pattern_expression'` FOLLOWED BY an fmri_data pattern object, `'correlation'`, `'cosine_similarity'`, `'norm_mask'`, `'ignore_missing'`.
- **Outputs:** `parcel_means` `[images x parcels]` (always, full-width, NaN for lost parcels); `parcel_pattern_expression`; `parcel_valence` (cosine vs unit vector); `rmsv_pos`/`rmsv_neg` (weights/cm³, vol-regularized by +1 cm³); `voxel_count`; `parcel_ste` (only for statistic_image with `.ste`).
- **When to use:** Direct, fastest path for atlas/parcel means and per-parcel pattern application; the engine behind `@atlas/extract_data` and `extract_measures_batch`.
- **Calls:** local `match_spaces`; `condf2indic`; `get_region_volumes`; `canlab_pattern_similarity`. Sets `parcels.probability_maps = []` first for speed.
- **Caveat:** One parcel image / one pattern at a time. Space order = pattern → (data → atlas).

#### `@image_vector/apply_mask.m`
- **Path:** `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/@image_vector/apply_mask.m`
- **Purpose:** Subset voxels to a single (global) mask, OR compute whole-image pattern expression.
- **Signature:** `function [dat, mask] = apply_mask(dat, mask, varargin)`
- **Inputs:** `dat` (image_vector & subclasses); `mask` (filename / fmri_mask_image / weighted image_vector / statistic_image). Options: `'pattern_expression'`, `'correlation'`, `'cosine_similarity'`, `'norm_mask'` (L2), `'ignore_missing'`, `'invert'`.
- **Outputs:** `dat` = masked object OR similarity output; `mask` = resampled, empty-removed mask.
- **When to use:** Masking to one region, or a single global dot-product/similarity (not region-aware). The building block under `extract_gray_white_csf`.
- **Calls:** local `get_nonempty_voxels`; `zeroinsert`; `canlab_pattern_similarity`; `check_properties(..., 'compress_index')` for atlas inputs.
- **Note:** Most heavily used masking primitive (~57 call sites across 13 method files). Doc-quality flagged (`poorly_documented_functions.txt:27`) — keep, improve header.

#### `@image_vector/extract_gray_white_csf.m`
- **Path:** `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/@image_vector/extract_gray_white_csf.m`
- **Purpose:** Per-tissue (gray/white/CSF) summary statistic across canonical tissue masks, with optional PCA components and L2 norms.
- **Signature:** `function [values, components, full_data_objects, l2norms] = extract_gray_white_csf(obj, varargin)`
- **Inputs:** `obj`. Options: `'eval'` + function handle (default `@(x1)(nanmean(x1,1))`), `'masks'` + 3-cell {gray, white, ventricle} (defaults `gray_matter_mask_sparse.img`, `canonical_white_matter.img`, `canonical_ventricles.img`).
- **Outputs:** `values` `[n_images x 3]`; `components` (1x3 cell of 5 PCs each); `full_data_objects` (1x3 cell); `l2norms` `[n_images x 3]`. Lazy by `nargout`.
- **When to use:** QC / nuisance tissue summaries in standard MNI space.
- **Calls:** `apply_mask` once per tissue; local `getnorms`.

#### `@fmri_data/extract_measures_batch.m`
- **Path:** `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/@fmri_data/extract_measures_batch.m`
- **Purpose:** Orchestrator that bundles outlier stats, tissue compartments, signature responses, and per-parcel means/patterns into one `DAT` struct.
- **Signature:** `function DAT = extract_measures_batch(data_obj)` (no varargin).
- **Inputs:** `data_obj` (fmri_data).
- **Outputs:** `DAT` struct (mahalanobis, rmssd, gray_white_csf_table, signature responses, `PARCELS` per parcellation with group t/p/FDR).
- **When to use:** One-shot QC + measures dump for a dataset.
- **Calls:** `mahal`, `preprocess`, `extract_gray_white_csf`, `apply_all_signatures`, `resample_space`, `load_atlas`, and **`apply_parcellation` directly** (NOT `atlas.extract_data`).

#### `@region/check_extracted_data.m`
- **Path:** `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/@region/check_extracted_data.m`
- **Purpose:** QC sanity check — re-read original images and confirm stored region averages still correlate > 0.999.
- **Signature:** `function isok = check_extracted_data(cl)`
- **Inputs:** `cl` (region with `.dat`, `.XYZ`, `.source_images` on path).
- **Outputs:** `isok` logical (NOTE: overwritten each loop, reflects only the last sampled region).
- **When to use:** After extraction, to verify integrity.
- **Calls:** `spm_get_data`, `check_valid_imagename`, `corrcoef`.

> **`apply_atlas` does NOT exist.** It is named only in the help headers of `@image_vector/image_vector.m` (lines 23, 143) and `@fmri_data/fmri_data.m` (lines 23, 141). The functional equivalent is `apply_parcellation` (+ `@atlas/extract_data`). Do not document it as callable.

---

### (b) Standalone current functions

| Function | Path | Purpose / signature | When to use |
|---|---|---|---|
| `extract_image_data` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/extract_image_data.m` | Generic mask/atlas extraction → voxel data + region averages. `function [imgdat, volInfo, cl] = extract_image_data(imgs_to_extract_from, mask_image, varargin)`. Default `'contiguous_regions'`. | Prefer the OO methods; header itself points to `@fmri_data/extract_roi_averages`. **Orphaned (0 live callers)** — removal candidate. |
| `canlab_maskstats` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/canlab_maskstats.m` | Pattern/weight-mask stats vs image sets. `function MASKSTATS = canlab_maskstats(msks, imgs, varargin)`. Measures: mean/std/dot_product/cosine_similarity/correlation/centered_dot_product/all. | Still wired into `canlab_glm_maskstats`. Otherwise prefer `apply_mask 'pattern_expression'` / `apply_nps` / `canlab_pattern_similarity`. |
| `canlab_load_ROI` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/canlab_load_ROI.m` | Load named published ROIs as `region`/`atlas` objects. `function [r, atlas_obj, default_color, region_file, image_file] = canlab_load_ROI(region_name, varargin)`. Option `'noatlas'`. | **KEEP.** Object-aware; interoperates with `addbrain`. |
| `canlab_extract_ventricle_wm_timeseries` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_processing_tools/canlab_extract_ventricle_wm_timeseries.m` | aCompCor-style ventricle+WM nuisance timeseries. `function [vw_nuisance, vw_nuisance_comps] = canlab_extract_ventricle_wm_timeseries(mask_image_dir, imgs, varargin)`. Option `'noplot'`. | **KEEP.** Purpose-built for nuisance regression (not replaced by `extract_gray_white_csf`). |
| `timeseries_extract_slice` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/timeseries_extract_slice.m` | Single slice across all volumes → `X×Y×time`. `function sl = timeseries_extract_slice(V, sliceno, orientation)`. | **KEEP.** Low-level; no object equivalent. |
| `sphere_roi_tool_2008` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/sphere_roi_tool_2008.m` | Extract data from a sphere around a coordinate (default r=10mm). `function [cl, all_data] = sphere_roi_tool_2008(imgs, varargin)`. | **KEEP** (active, bridges to `region` via `cluster2region`). |
| `extract_from_rois` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/extract_from_rois.m` | Near-duplicate of `extract_image_data`; default `'unique_mask_values'`, output order `[cl, imgdat]`. | **Orphaned duplicate (0 callers).** Removal candidate. |

---

### (c) Legacy / grandfathered / duplicate

Operate on the pre-object **`clusters` struct** (from `tor_extract_rois`) unless noted. Replaced by `region` objects + `extract_roi_averages` / `extract_data`.

| Function | Path | Status | Notes |
|---|---|---|---|
| `tor_extract_rois` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/tor_extract_rois.m` | **KEEP-LEGACY-WRAPPER** | Root of old clusters pipeline. 34 live call sites across legacy procedural tools (`mask2clusters`, `cluster_orthviews`, `mask_*`, `parcel_*`, `hewma_*`). Never called from an `@class` method. The long pole for migration. Works on raw XYZ voxel coords, no space transform. |
| `extract_raw_data` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/extract_raw_data.m` | KEEP-LEGACY → removal candidate | Needs `EXPT`+`clusters`. No live external caller. Head of a legacy chain (calls `tor_extract_rois`, `extract_indiv_peak_data`). |
| `extract_contrast_data` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/extract_contrast_data.m` | KEEP-LEGACY-WRAPPER | Used only by GUI `cluster_tool` (2 sites). Chains to `tor_extract_rois` + `extract_ind_peak`. |
| `extract_indiv_peak_data` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/extract_indiv_peak_data.m` | Removal candidate | Sole caller is the orphaned `extract_raw_data`. |
| `cluster_tmask` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/cluster_tmask.m` | Legacy helper | For `extract_indiv_peak_data`. |
| `extract_ind_peak` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/peak_coordinates/extract_ind_peak.m` | KEEP-LEGACY-WRAPPER | Used by `cluster_barplot`, `extract_contrast_data`. Name-collision hazard with `extract_indiv_peak_data`. |
| `timeseries4` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/Grandfathered/timeseries4.m` | Removal candidate | In `Grandfathered/`. No external caller — live `timeseries4(...)` calls resolve to a local subfunction inside `cluster_orthviews.m`. |
| `check_timeseries_vals` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/Grandfathered/check_timeseries_vals.m` | Removal candidate | In `Grandfathered/`. Only used by sibling `timeseries4`. |
| `sphere_roi_tool_2008_conflict.m` | `/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/Data_extraction/sphere_roi_tool_2008_conflict.m` | **DELETE IMMEDIATELY** | Contains unresolved merge-conflict markers (`<<<<<<< .mine` / `=======` / `>>>>>>> .r1197`). Uncallable, duplicate function name. Zero risk. |

---

## 3. "Which method should I use?" decision guide

| Task | Recommended method |
|---|---|
| Region means from subject-level `fmri_data` (filename/mask/atlas ROIs) | `@fmri_data/extract_roi_averages` |
| Region means from a generic `image_vector` subclass | `@image_vector/extract_roi_averages` |
| Per-region **pattern expression** (weighted similarity) | `@fmri_data/extract_roi_averages` with `'pattern_expression'`, or `apply_parcellation` with `'pattern_expression', pattern_obj` |
| Whole-atlas parcel means / valence / RMS (atlas is primary object) | `@atlas/extract_data` |
| Parcel means fastest, data is primary object | `apply_parcellation(dat, atlas)` |
| Attach data to existing `region` objects without interpolation | `@region/extract_data` |
| Single global mask subset, or one whole-image dot product | `@image_vector/apply_mask` |
| Gray/white/CSF compartment summaries (QC) | `@image_vector/extract_gray_white_csf` |
| aCompCor ventricle+WM nuisance regressors (timeseries) | `canlab_extract_ventricle_wm_timeseries` |
| Sphere ROI around a coordinate | `sphere_roi_tool_2008` |
| Load a named published ROI as an object | `canlab_load_ROI` |
| One-shot QC + measures dump | `@fmri_data/extract_measures_batch` |
| Verify a prior extraction's integrity | `@region/check_extracted_data` |
| Single slice across volumes | `timeseries_extract_slice` |
| Pattern/weight-mask stats inside the GLM batch pipeline | `canlab_maskstats` (legacy; else `apply_mask`/`canlab_pattern_similarity`) |

---

## 4. Space-handling notes (the key gotcha)

**Resampling direction differs by method, with numerical consequences.**

| Method | Direction | Interpolation | Consequence |
|---|---|---|---|
| `@fmri_data/extract_roi_averages` | mask/atlas → DATA | `'nearest'` for `'unique_mask_values'` & atlas-aware for atlas masks; else linear | Linear interp would corrupt integer codes (Wani's fix), so codes use nearest. Plain-mean masks use linear. |
| `@image_vector/extract_roi_averages` | mask → DATA | default (linear) | `isdiff==2` (missing volInfo) → error. |
| `apply_parcellation` | parcels → DATA (pattern → data → atlas order) | `'nearest'` (preserves integer labels) | Full-width output via `condf2indic` with `n_orig_parcels`; NaN for lost parcels. |
| `@atlas/extract_data` | **DATA → ATLAS** (then `apply_parcellation` re-matches) | nearest in the wrapped `apply_parcellation` | Opposite of all others. Resampling the *data* changes the data grid; pre-resample atlas to functional space for large N. |
| `apply_mask` | mask → DATA | (resample_space default) | `isdiff==3` (different non-empty voxels) is OK; fixes illegal `removed_voxels` length post-resample. |
| `@region/extract_data` | **none** — mm-coordinate intersection | n/a | No interpolation artifacts; results differ slightly from resample-first methods by design. Missing voxels → NaN. |
| `extract_gray_white_csf` | delegates to `apply_mask` (mask → data) | per `apply_mask` | After masking, `.dat==0` set to NaN before stats. |

**Rules of thumb:**
- Integer-labeled atlases/parcellations must use **nearest-neighbor** resampling. The OO methods handle this; if you roll your own, do not linear-interp integer codes.
- `@region/extract_data` is the choice when you want to avoid any interpolation (coordinate truth), at the cost of being slightly stricter about exact voxel matches.
- `@atlas/extract_data` moving DATA→ATLAS means downstream stats are computed on the resampled data grid — be deliberate when comparing its outputs to mask→data methods.

---

## 5. Call-graph summary & redundancies

### Internal call graph (who calls whom)
```
@atlas/extract_data ─────────────► apply_parcellation
@fmri_data/extract_measures_batch ► apply_parcellation (directly, NOT atlas.extract_data)
                                  ► extract_gray_white_csf ► apply_mask
                                  ► mahal, preprocess, apply_all_signatures, resample_space, load_atlas
apply_parcellation ──────────────► match_spaces (local), canlab_pattern_similarity, condf2indic, get_region_volumes
@fmri_data/extract_roi_averages ─► region()/atlas2region(), canlab_pattern_similarity
@image_vector/extract_roi_averages► region()/atlas2region()   [errors on 'pattern_expression']
apply_mask ──────────────────────► get_nonempty_voxels (local), zeroinsert, canlab_pattern_similarity
@region/extract_data ────────────► mm2voxel, match_coordinates (local), canlab_pattern_similarity
@region/check_extracted_data ────► spm_get_data, check_valid_imagename
@region/region (constructor) ────► extract_roi_averages

Legacy: extract_raw_data ► tor_extract_rois, extract_indiv_peak_data
        extract_contrast_data ► tor_extract_rois, extract_ind_peak
        extract_indiv_peak_data ► cluster_tmask
```

KEEP-CORE trio/quartet wired into the OO API: `extract_roi_averages`, `extract_data`, `apply_parcellation`, `apply_mask`. `@region/region.m` builds region objects via `extract_roi_averages`; `@atlas/extract_data.m` is the canonical caller of `apply_parcellation`.

### Redundancies (feed-forward to the streamlining plan — handled in a separate file)
1. **`sphere_roi_tool_2008_conflict.m`** — broken merge-conflict file, zero callers. Delete now.
2. **Twin orphaned ROI extractors** — `extract_image_data.m` and `extract_from_rois.m` share H1 text and job; neither has a live caller (the one `extract_image_data` call resolves to a same-file local subfunction in `parcel_images.m`). Collapse to one or remove both; update the dead "non-object-oriented alternative" pointers in `@fmri_data/extract_roi_averages.m:96` and `@atlas/extract_data.m:114`.
3. **`Data_extraction/Grandfathered/` is fully dead** — `timeseries4.m`, `check_timeseries_vals.m` have no external callers (live calls resolve to local subfunctions in `cluster_orthviews.m`). Directory removable.
4. **Legacy clusters chain** — `extract_raw_data` → (`tor_extract_rois`, `extract_indiv_peak_data` → `cluster_tmask`); `extract_contrast_data` → (`tor_extract_rois`, `extract_ind_peak`). Group as KEEP-LEGACY, schedule migration to `extract_roi_averages`. `tor_extract_rois` (34 sites) stays until `mask2clusters`/`cluster_orthviews`/`mask_*`/`parcel_*`/`hewma_*` are migrated.
5. **`canlab_maskstats`** overlaps with `apply_mask 'pattern_expression'` + `extract_roi_averages`; keep only until `canlab_glm_maskstats` is refactored.
6. **Name-collision hazards** — `extract_ind_peak` vs `extract_indiv_peak_data`; local subfunction shadows for `tor_extract_rois`, `timeseries4`, `check_timeseries_vals`, `extract_image_data` inside `cluster_orthviews.m` / `parcel_images.m`. Verify shadowing is intentional before deleting any standalone.

**Doc-quality flags** (`/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore/poorly_documented_functions.txt`): `apply_mask` (27), `check_timeseries_vals` (200), `timeseries4` (201), `sphere_roi_tool_2008` + conflict twin (204–205), `extract_raw_data` (521).
