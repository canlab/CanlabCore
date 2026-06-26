# CANlab Surface / Grayordinate Data Object — Design & Implementation Plan

## 1. Purpose and scope

This document specifies the design and a phased implementation plan for a new CANlab
object, **`fmri_surface_data`**, that represents cortical-surface (vertex) and
grayordinate (CIFTI: surface vertices + subcortical voxels) brain data inside the
existing CANlab object model. The goal is a first-class sibling of `fmri_data` that
(a) losslessly holds and round-trips CIFTI `.dscalar/.dtseries/.dlabel.nii` and GIFTI
`.surf/.func/.shape/.label.gii` files, (b) reuses the entire geometry-agnostic CANlab
analysis stack (`predict`, `ica`, `mean`, `descriptives`, parcellation,
etc.) by keeping data in the canonical flat `[units × images]` `.dat` matrix, (c) renders
natively on bundled fs_LR / fsaverage meshes, and (d) maps to/from volumetric `fmri_data`
using warps **already vendored in the repo** — all with **no external MATLAB toolbox,
FreeSurfer, or Connectome Workbench binary required at runtime**. Scope for v1 is the
object, native I/O, native rendering, group-template volume↔surface mapping, and surface
parcellation; per-subject ribbon-constrained mapping and fsaverage↔fs_LR deformation are
explicitly deferred.

This plan reconciles three independent architecture proposals. Where they diverged, the
conflicts and resolutions are stated explicitly in §3.

---

## 2. Background: CIFTI grayordinates and GIFTI for CANlab developers

**The grayordinate model.** HCP "grayordinates" represent the cortex as *surface
vertices* and the subcortex/cerebellum as *voxels*. The standard "91k" space has
**91,282 grayordinates** = 29,696 left-cortex vertices + 29,716 right-cortex vertices
(the non-medial-wall subset of 32,492 fs_LR-32k vertices per hemisphere) + 31,870
subcortical 2 mm MNI voxels spread across ~19 `CIFTI_STRUCTURE_*` volume structures
(accumbens, amygdala, caudate, cerebellum, diencephalon, hippocampus, pallidum, putamen,
thalamus L/R, brainstem). The key insight for CANlab: **grayordinate data is already a
flat `[grayordinates × maps]` matrix** — exactly the layout CANlab's `.dat` expects.

**CIFTI-2 file format.** A CIFTI-2 file is a single uncompressed `.nii` (NIfTI-2)
container: a 540-byte little-endian NIfTI-2 binary header, a 4-byte extension flag at
byte 540, one or more header extensions (the CIFTI XML lives in the extension whose
`ecode == 32`), then the raw data matrix beginning at `vox_offset`. NIfTI `dim[1..4]=1`;
the real CIFTI lengths live in `dim[5]` (values per row, fastest) and `dim[6]` (rows);
`dim[0]` is 6 or 7. **The data matrix is stored row-major**, so a MATLAB reader `fread`s
then permutes `[2 1]` to column-major. The XML (`<CIFTI Version="2"><Matrix>...`) carries
one `MatrixIndicesMap` per dimension; a `BRAIN_MODELS` map holds `BrainModel` elements
each with `IndexOffset`, `IndexCount`, `ModelType` (`SURFACE`/`VOXELS`), `BrainStructure`,
`SurfaceNumberOfVertices`, and either `<VertexIndices>` (0-based vertices that carry data,
i.e. medial wall excluded) or `<VoxelIndicesIJK>` + a `<Volume>` 4×4 affine.
Intent codes / extensions: `.dscalar.nii` = 3006 `ConnDenseScalar`, `.dtseries.nii` =
3002 `ConnDenseSeries`, `.dlabel.nii` = 3007 `ConnDenseLabel` (integer keys indexing a
per-map `LabelTable`). **No external toolbox is required to read/write CIFTI-2** — only
`fopen/fread/fwrite` + XML parsing (`xmlread` or regexp). `wb_command` is needed only to
convert legacy CIFTI-1, which v1 does not support.

**GIFTI (.gii) format.** GIFTI is a UTF-8 XML container. Root
`<GIFTI Version="1.0" NumberOfDataArrays="N">` holds N `<DataArray>` elements. A geometry
`.surf.gii` has two arrays: `NIFTI_INTENT_POINTSET` (N×3 FLOAT32 vertex mm coords) and
`NIFTI_INTENT_TRIANGLE` (M×3 INT32 0-based faces — add 1 for MATLAB `patch`). Per-vertex
data files use `NIFTI_INTENT_NONE/SHAPE/TIME_SERIES/LABEL`. Each `<Data>` payload is
`ASCII`, `Base64Binary`, or `GZipBase64Binary` (base64 of a **raw zlib** stream, header
`0x78 0x9C` — **not gzip**). A fully native MATLAB reader/writer was implemented and
**validated bit-exact** on the repo's S1200 32k surf.gii (32,492 verts / 64,980 faces;
verts/faces `maxdiff = 0` on round-trip). Three gotchas were solved and must be carried
over verbatim: (1) decompress with `java.util.zip.InflaterInputStream` (not
`GZIPInputStream`); (2) base64 with `matlab.net.base64decode/encode` (not apache-commons,
which corrupts bytes); (3) inflate entirely inside the JVM via `IOUtils.copy` (in-place
`Inflater.inflate(buf)` returns garbage to MATLAB because MATLAB passes a copy).

**Why CANlab can do this with almost no new format code:** CANlab already ships fs_LR-32k
and fsaverage meshes as plain `.mat` files (faces + vertices), the MNI→fs_LR/fsavg
barycentric resample matrices, and the CBIG MNI↔fsaverage registration-fusion warps — so
the only genuinely new I/O is a self-contained CIFTI/GIFTI parser.

---

## 3. Design decisions & rationale

### D1. Class name — `fmri_surface_data`
All three proposals independently converged on `fmri_surface_data`. **Decision:
`fmri_surface_data`.** It is unambiguous, prefixed to avoid collisions, and signals
"surface/grayordinate analogue of `fmri_data`." *(Open question Q1: confirm the
`canlab_` prefix vs. a bare `surface_data`/`fmri_surface_data`; the existing siblings
`fmri_data`/`image_vector`/`atlas`/`region` are unprefixed, so a bare name like
`fmri_surface_data` would match house style. Recommend confirming with the user.)*

### D2. Superclass — subclass `image_vector`, **not** `fmri_data`, **not** from scratch
All three proposals agree: **`classdef fmri_surface_data < image_vector`.**

- Subclassing `image_vector` inherits the load-bearing, purely-`.dat`/`removed_*`
  methods verbatim by dispatch: `remove_empty`, `replace_empty`, `get_wh_image`,
  `descriptives`, `ica`, `mahal`, `pca`, `image_math`, `history`, `isempty`,
  `enforce_variable_types`.
- **Not** `fmri_data`: `fmri_data` hard-codes a mandatory `volInfo` and an
  `fmri_mask_image`-typed `mask` in `run_checks_and_fixes` (errors if `volInfo` empty),
  and its `cat/horzcat/plot/predict/ttest/regress` assume `interp3`/orthviews resampling.
  Inheriting those would fight the surface design at every turn.
- **Not** from scratch: that would forfeit dozens of inherited methods and the
  constructor polymorphism contract.

Trade-off accepted: `cat`, `horzcat`, `plot`, `predict`, `ttest`, `regress` live **only**
in `@fmri_data`, so an `image_vector` subclass does **not** inherit them. We provide our
own (see §5), copying the geometry-agnostic cores and stripping the volume rebuild/resample
steps. This is the deliberate price of not inheriting `fmri_data`'s volume baggage.

### D3. How `volInfo` is replaced — a `brain_model` descriptor + repurposed `volInfo` slot

This is where the proposals **conflicted** and is resolved explicitly:

- **Proposal A** repurposes the inherited `volInfo` *property slot* to hold the whole
  surface+volume `brain_model` descriptor (keeping the literal name `volInfo`).
- **Proposal B/C** add a **new** `brain_model`/`cifti`/`geom` property as the single
  source of truth, and (B) additionally keep `volInfo` populated to describe **only the
  volumetric sub-block** so the subcortex still round-trips through existing voxel code.

**Resolution (graft the best of both):**
1. The surface/grayordinate geometry truth lives in a **new property `brain_model`** (a
   CIFTI-`BrainModel`-mirroring struct; the single source of truth for round-trip and for
   the flat-row↔space map). This is cleaner and less leak-prone than overloading `volInfo`
   semantically (Proposal A's chief risk was inherited methods dereferencing
   `volInfo.mat`).
2. The inherited **`volInfo` slot is kept populated, but describes ONLY the volumetric
   sub-block** of the grayordinates (its `.mat`/`.dim`/`.xyzlist`/`.wh_inmask` come from
   the subcortical `VoxelIndicesIJK` + affine). This is Proposal B's key idea: it means
   the subcortex is a *valid voxel space*, inherited methods that read `volInfo` find a
   real (if partial) struct, and `to_fmri_data` on the subcortex is almost free. For
   **surface-only** objects (`.func.gii`, cortex-only CIFTI), `volInfo` is empty and the
   relaxed run-checks tolerate it.
3. Mesh geometry (faces/vertices per hemisphere, medial wall) lives in a **new property
   `geom`**, loaded lazily from bundled `.mat` assets (this is what render/area/contiguity
   patch and compute on).

`brain_model` fulfills `volInfo`'s three load-bearing roles: **(a)** row↔space map =
per-structure index lists (`vertlist` 0-based per hemisphere; `voxlist` IJK per
subcortical structure) — the `wh_inmask`/`xyzlist` analogue; **(b)** reconstruct target =
scatter rows into dense per-hemisphere vertex arrays (`reconstruct_image` override) +
subcortical volume; **(c)** contiguity source = mesh-graph connected components (cortex) /
3-D connectivity (subcortex), replacing `volInfo.cluster`.

### D4. Re-declare `fmri_data` per-image annotation properties **by name**
Because we subclass `image_vector` (where these don't exist), we re-declare `X`, `Y`,
`Y_names`, `covariates`, `covariate_names`, `images_per_session`, `metadata_table`,
`image_metadata`, `additional_info` with the **identical `fmri_data` names** so analysis
methods bind by name. These are per-*image* (column) annotations, independent of
voxel-vs-vertex.

### D5. Preserve method contracts exactly
`compare_space` must keep its **4-valued return code** (`0` same / `1` diff / `2` missing /
`3` same-space-diff-grayordinates), not a boolean — callers (`cat`, `apply_mask`) branch
on it. A single overridable **`rebuild_like(obj, newdat)`** helper replaces the hard-coded
`image_vector('dat',…,'volInfo',…)`+`fmri_data` re-wrap in `mean`/`ica`/etc. so they emit
`fmri_surface_data` carrying `brain_model`+`geom`.

### D5b. NO empty-squeezing — `.dat` is always the full grayordinate set (user directive, 2026-06-25)
Unlike `fmri_data` (which stores sparse 3-D volumes and must drop out-of-brain / empty
voxels to save space), CIFTI grayordinate data is **already compact**: the medial wall is
excluded by construction (it is simply absent from each surface model's `vertlist`), and
there are no out-of-brain voxels. So this object does **not** reproduce the
`remove_empty`/`replace_empty` space-squeezing contract:

- `.dat` is **always** `[nGrayordinates × nMaps]`, in exact 1:1 row correspondence with
  `brain_model`. It is never shrunk to drop zero/NaN rows. There is therefore no "full vs
  reduced" duality and no `replace_empty`-before-reconstruct dance.
- `remove_empty` and `replace_empty` are **overridden as identity no-ops** (return the
  object unchanged) so any inherited method that calls them internally still composes, and
  a user who calls them by habit gets a valid object back. (Zeroing a grayordinate's value
  is allowed; it just never removes the row.)
- `removed_voxels` is kept (inherited) but is **vestigial**: always `false(nGray,1)`. It
  exists only so inherited methods that read it find a valid all-false vector. `removed_images`
  likewise stays `false(nImg,1)`; *selecting* a subset of maps is `get_wh_image`'s job, not a
  space-saving squeeze.
- **Masking is intrinsic, not space-squeezing.** Because all objects share a standardized
  grayordinate space (fs_LR-32k / 91k), a "mask" is just a same-length logical over
  grayordinates (or another `fmri_surface_data` on the same space). `apply_mask` therefore
  zeros/selects matching rows directly — no `fmri_mask_image`, no resample, no
  empty-removal. This is a major simplification over `@fmri_data/apply_mask`.

This eliminates the plan's largest risk class (the "forgot to `replace_empty`" family) and
several override headaches; `reconstruct_image`/`write`/`cat` no longer need a
`replace_empty` pre-step.

### D6. Statistic / atlas variants deferred to a later phase
`ttest`/`regress` results and `.dlabel` parcellations want parallel per-grayordinate
fields (`.p/.ste/.sig`; integer label keys + RGBA). For v1, `ttest`/`regress` return a
`fmri_surface_data` carrying those as extra fields; dedicated subclasses
(`fmri_surface_statistic_image`, `fmri_surface_atlas`) that keep those parallel fields in
sync under `get_wh_image` (column subsetting) are a later phase. *(No `remove_empty`/
`replace_empty` sync needed — they are no-ops, D5b.)*

### Open questions — RESOLVED (user, 2026-06-25)
- **Q1. Class name → `fmri_surface_data`** (unprefixed, matching house style of
  `fmri_data`/`atlas`/`region`). Statistic/atlas subclasses → `fmri_surface_statistic_image`
  / `fmri_surface_atlas`.
- **Q2. Canonical default space → fs_LR-32k (91k)** for native CIFTI; fsaverage-164k for the
  CBIG mapper path. Confirmed.
- **Q3. Mesh assets → ship the in-repo HCP S1200 meshes for v1, add CC0 TemplateFlow
  `tpl-fsLR`/`tpl-fsaverage` as the license-clean default in a later milestone (M8).**
- **Q4. CIFTI I/O → hand-roll natively** (`canlab_read_cifti`/`canlab_write_cifti`), modeled
  on the BSD-2 cifti-matlab core, so CanlabCore stays 100% permissive and self-contained.
  The GIFTI codec is already hand-rolled and validated bit-exact.

---

## 4. Object model

`.dat` is `single [n_grayordinates × n_maps]`; rows in fixed CIFTI order: left-cortex
non-medial vertices, then right-cortex non-medial vertices, then subcortical voxels.

| Property | Type | Description |
|---|---|---|
| `dat` | `single [nGray × nImg]` | **Inherited, unchanged.** The CIFTI cdata matrix (cortex L, cortex R, subcortical voxels). All generic `.dat` math operates here. For `.dlabel` holds integer label keys. |
| `brain_model` | struct (**new**, the `volInfo` replacement) | **Single source of truth** for surface/grayordinate geometry, mirroring CIFTI `BrainModel` 1:1. Fields: `.models` (struct array per BrainModel: `.struct_name` e.g. `CIFTI_STRUCTURE_CORTEX_LEFT`, `.type` `'surf'|'vox'`, `.index_offset`, `.index_count`, `.vertlist` 0-based vertex indices, `.surf_nverts`, `.voxlist` 3×N IJK, `.mat` 4×4 affine for vox models); `.grayordinate_type` e.g. `'91k'`; `.cluster` (computed connected-component labels). Replaces `wh_inmask/xyzlist` with `vertlist/voxlist`; per-vox-model affines replace one global `mat`. |
| `geom` | struct (**new**, mesh cache) | Per-hemisphere meshes (loaded lazily from bundled `.mat`): `.faces_lh/.faces_rh` [M×3 1-based], `.vertices_lh/.vertices_rh` [N×3 mm] for one or more surface types (`midthickness/inflated/pial/white/sphere`, faces shared across types), `.space_name` (`fsLR_32k`/`fsavg_164k`), `.medialwall_lh/.medialwall_rh` logical masks. What render/area/contiguity use. |
| `volInfo` | struct (**inherited, repurposed for the volume sub-block only**) | Populated to describe ONLY the subcortical/cerebellar voxel models (`.mat`=their affine, `.dim`, `.xyzlist/.wh_inmask` from `voxlist`) so the subcortex is a valid voxel space and `to_fmri_data` is near-free. **Empty for surface-only objects** (run-checks relaxed). Surface part is described by `brain_model`+`geom`, never `volInfo`. |
| `removed_voxels` | logical [nGray × 1] | **Inherited but vestigial (see D5b).** Always `false` — grayordinate data is already compact, so rows are never squeezed. Kept only so inherited methods that read it find a valid all-false vector. |
| `removed_images` | logical [nImg × 1] | **Inherited but vestigial (see D5b).** Always `false`. Map subsetting is `get_wh_image`'s job, not a space-saving squeeze. |
| `history` / `source_notes` / `image_names` / `fullpath` / `files_exist` / `dat_descrip` | cell / char / logical | **Inherited verbatim.** Provenance + file metadata. `fullpath` reinterpreted by `write()` as `.dscalar/.dtseries/.dlabel.nii` or `.gii` path; extension drives intent. |
| `X` | double [nImg × p] | **New (fmri_data name).** Design/predictor matrix, per-map, so `regress`/`predict` bind. |
| `Y` / `Y_names` | double [nImg × q] / cell | **New (fmri_data names).** Outcomes + names, per-map. |
| `covariates` / `covariate_names` | double / cell | **New (fmri_data names).** Per-map nuisance covariates. |
| `images_per_session` | double vector | **New (fmri_data name).** Run lengths for `.dtseries` data. |
| `metadata_table` | table (nImg rows) | **New (fmri_data name).** Per-map metadata; `write()` exports a companion CSV like `fmri_data.write`. |
| `image_metadata` / `additional_info` | struct | **New (fmri_data names).** Acquisition flags + free-form info; preserves the `fmri_data` API surface. |
| `mask` | logical [nGray × 1] or `fmri_surface_data` | **New, lightweight (see D5b).** Optional same-length logical over grayordinates (or another object on the same space). No `fmri_mask_image`, no resampling, no empty-squeezing — masking just zeros/selects rows. Medial wall is already excluded by `brain_model`, so a `mask` is only for further sub-selection (e.g. cortex-only, an ROI). |
| `intent` | char | **New.** CIFTI/GIFTI intent (`dscalar`/`dtseries`/`dlabel`/`func`/`shape`/`label`) so `write()` round-trips the correct `intent_code` + extension. |
| `series_info` | struct | **New.** For `.dtseries`: `SeriesStart/Step/Unit/Exponent` from the CIFTI SERIES map. |
| `label_table` | struct/table | **New (only `intent=='dlabel'/label`).** Per-key Name + RGBA, mirroring atlas label tables; lets a surface `.dlabel` convert cleanly to/from a surface atlas. |
| `surface_space` | char | **New.** Canonical space tag (`'fsLR_32k_91k'` default, `'fsLR_32k'`, `'fsaverage_164k'`). Gatekeeps `compare_space`; drives which meshes/warps load. |

**How `brain_model` plays `volInfo`'s role across surface + volume grayordinates:** for
**cortical** grayordinates there is no affine — the row↔space map is
`models(i).vertlist` into `models(i).surf_nverts` (per hemisphere), reconstruction
scatters rows into dense vertex arrays (medial wall → NaN), and contiguity is connected
components on the `geom.faces` edge graph. For **subcortical** grayordinates the
`models(i).voxlist` + `models(i).mat` *is* a normal voxel space — these are mirrored into
the inherited `volInfo` slot so the subcortex reconstructs into a 3-D volume and exports
to `fmri_data` via existing voxel code.

---

## 5. Method surface

Classification: **Inherited** (works verbatim by dispatch), **Override** (same name,
surface-specific body), **New** (not present in `@image_vector`; mirrors `@fmri_data`).
"Mirrors fmri_data" = same name/signature/role as the `fmri_data` method.

| Method | Class | Behavior (one line) |
|---|---|---|
| `fmri_surface_data` (constructor) | New | Mirrors `image_vector`/`fmri_data` polymorphism: empty / struct→fieldname-copy / `'key',value` / filename-autoload (`.dscalar/.dtseries/.dlabel.nii`/`.gii` extension test replaces `'.nii'`+exist) / `isa(image_vector)` recast (copy matching fields, cast `.dat` to single). Builds `brain_model`+`volInfo`(vol sub-block); lazy `geom`. Relaxes mandatory-`volInfo` check. |
| `remove_empty` | Override → **no-op** (D5b) | Returns the object unchanged. Grayordinate data is already compact; rows are never squeezed. Defined so inherited callers compose and habitual user calls are safe. |
| `replace_empty` | Override → **no-op** (D5b) | Returns the object unchanged (`.dat` is always full-size already). Removes the "forgot to replace_empty" bug class. |
| `get_wh_image` | Inherited | Map (column) selection; subsets `.dat`-sized + per-image fields (X/Y/metadata_table). |
| `descriptives` | Inherited | Numeric `.dat` summaries; only the "n voxels" label reads cosmetically as grayordinates. |
| `ica` / `mahal` / `pca` | Inherited | Pure `.dat` decomposition/stats; geometry-agnostic. Re-wrap (if any) via `rebuild_like`. |
| `image_math` / `history` / `isempty` / `enforce_variable_types` | Inherited | Pass `.dat`/metadata around; bind unchanged. |
| `mean` | Override (math reused) | Reuse `.dat` averaging + `group_by`; replace the hard-coded `image_vector(...,'volInfo',...)`+`fmri_data` re-wrap with `rebuild_like` → `fmri_surface_data` carrying `brain_model`+`geom`. |
| `threshold` | Override (raw branch reused) | Raw value thresholding on `.dat` delegates to inherited logic; override only the cluster-extent `k` / `trim_mask` paths to mesh connected components (cortex) + 3-D (subcortex), medial-wall aware. |
| `apply_mask` | Override (simplified, D5b) | Same-space logical/object mask → zero or select matching `.dat` rows directly. No `fmri_mask_image`, no resample, no empty-removal. Pattern-expression / dot-product math on `.dat` reused verbatim where applicable. |
| `reconstruct_image` | Override | Scatter `.dat` rows directly (no `replace_empty` pre-step; `.dat` is always full) into dense per-hemisphere vertex arrays via `brain_model.models.vertlist` (medial wall→NaN) and subcortical models into a 3-D volume via their `.mat`. Returns `{lh_vertdata, rh_vertdata, subvol}`. |
| `resample_space` | Override | `interp3`-in-mm is meaningless on meshes. Surface: barycentric mesh↔mesh via bundled `resample_from_*_to_*.mat` structs; subcortex: existing voxel resample. Same signature. |
| `compare_space` | Override | Compare `surface_space` + per-model `struct_name`/`index_count`/`surf_nverts` + vox-model `dim`/`mat`. **Preserves the 0/1/2/3 return contract.** |
| `reparse_contiguous` | Override | Connected components on the mesh edge graph (cortex, via `External/matlab_bgl` or `graph/conncomp`) + 6/18/26 voxel connectivity (subcortex). Writes `brain_model.cluster`. |
| `write` | Override | Native `canlab_write_cifti`/`canlab_write_gifti` (M1), dispatched on `intent`/extension; keeps `fullpath/fname/overwrite` convention; exports `metadata_table` CSV. No `replace_empty` pre-step needed (`.dat` already full). |
| `apply_parcellation` | Override (core reused) | `condf2indic`→column-normalize→`parcel_means = dat.dat' * parcels.dat` reused verbatim; override only space-matching (shared grayordinate index / mesh resample) and `get_region_volumes`→per-parcel surface **area** (face geometry) for `rmsv`; exclude medial wall from normalization. |
| `region` / `surface_region` | Override / New | Per-region vertex-index lists + per-vertex `.val/.Z` + centroid (mean vertex coord) via mesh components; subcortical models route through the existing volumetric `region` path. Mirrors `cifti_struct_2_region_obj` field conventions for the vox part (which currently drops `'surf'` — the gap to fill). |
| `surface` / `render_on_surface` | Override | Object already holds per-vertex values → skip `render_on_surface`'s `interp3`; feed `.dat` columns straight into `FaceVertexCData` + the split gray/color colormap builders. Reuse `surface_outlines.m`; reuse `make_surface_figure.m` body (swap `gifti()`→`load()`). Register into `@fmridisplay/surface.m` + `addbrain`. |
| `montage` / `orthviews` / `slices` | Override | Route cortex through the surface renderer; subcortex through the existing slice montage; warn where a 2-D slice view is meaningless for cortex. |
| `cat` | New (mirrors fmri_data) | Establish common grayordinate space via surface `compare_space` (resample if needed), hcat `.dat` + per-map fields (X/Y/covariates/metadata_table/series_info). No `replace_empty` pre-step (`.dat` always full). Reuses `@fmri_data/cat.m` logic minus volume resample. |
| `horzcat` | New (mirrors fmri_data) | Thin wrapper over `cat`. |
| `predict` | New (mirrors fmri_data) | Copy `@fmri_data/predict.m` algorithm core (operates on `.dat`+X/Y, geometry-agnostic); skip volume weight-map rebuild, wrap weights via `rebuild_like`. |
| `ttest` / `regress` | New (mirrors fmri_data) | Per-row stats on `.dat` columns (geometry-agnostic core); return a `fmri_surface_data` carrying parallel `.p/.ste/.sig` per-grayordinate fields (statistic subclass deferred). |
| `plot` | New (mirrors fmri_data) | Surface QC: reuse the geometry-agnostic panels (covariance, global signal, per-map histograms); replace the volume montage/orthviews panel with a 4-view surface render of the mean map. |
| `to_fmri_data` | New | Export the subcortical/volumetric models to a `fmri_data` (+ writeable `.nii`) via the repurposed `volInfo` sub-block, reusing `extract_vol_from_cifti.m`. Surface models dropped here (use `surf2vol`). |
| `vol2surf` (method on `fmri_data`/`image_vector`) | New | Volume → surface via CBIG RF (returns `fmri_surface_data`). See §7. |
| `surf2vol` | New | Surface → volume via CBIG RF (returns `fmri_data`; `.nii` via `fmri_data.write`). See §7. |
| `rebuild_like` (helper) | New | Overridable rebuild of a `fmri_surface_data` from a new `.dat`, carrying `brain_model`+`geom`. Used by `mean`/`ica`/`predict`/etc. to avoid the volume re-wrap. |

---

## 6. Import / Export — native, no external toolbox

**Hard requirement:** no `gifti`, FieldTrip `ft_read_cifti`, `wb_command`, or
cifti-matlab `@xmltree` (LGPL) / `ft_cifti` (GPL) at runtime. Verify at each milestone via
`mcp__matlab__detect_matlab_toolboxes` on a clean path. Recommend hand-rolling the
reader/writer (modeled on the BSD-2 cifti-matlab core) so CanlabCore stays 100%
permissive and self-contained.

### CIFTI-2 reader (`read_cifti_native`)
1. `fopen(fname,'r','l')`; `fread` `sizeof_hdr` int32 @0 (=540; else byte-swap).
2. Read the 540-byte NIfTI-2 header at exact offsets: `dim` int64@16, `datatype`
   int16@12, `vox_offset` int64@168, `scl_slope`@176/`scl_inter`@184, `intent_code`
   int32@504, `intent_name` char16@508.
3. `fseek(540)`; read int32 extension flag; loop esize/ecode blocks until `vox_offset`;
   keep the `ecode==32` bytes as the XML char array.
4. Parse XML with a **regexp pull-parser** (preferred — avoids even the `xmlread`/JRE
   dependency and namespace/entity edge cases) or `xmlread`. Walk
   `CIFTI>Matrix>MatrixIndicesMap`; for `BRAIN_MODELS` read each `BrainModel`'s
   `IndexOffset/Count/ModelType/BrainStructure/SurfaceNumberOfVertices` and `sscanf`
   `VertexIndices`/`VoxelIndicesIJK`; read `Volume`+16-float affine (`reshape(...,4,4)'`,
   row-major). For SCALARS/LABELS iterate `NamedMap>MapName` (+`LabelTable>Label`); SERIES
   `Start/Step/Exponent/Unit`.
5. `fseek(vox_offset)`; `fread([rowLen nRows], precision)`; permute `[2 1]` to
   column-major; apply `scl_slope/inter` when nonzero. Populate `.dat`, `brain_model`,
   `volInfo`(vol sub-block via the split logic of `get_cifti_data.m`), `intent`,
   `series_info`, `label_table`, `image_names`.

### CIFTI-2 writer (`write_cifti_native`)
Build XML (`sprintf`, escape `&<>`), pad to /16, `esize=len+8` rounded, `ecode=32`;
`vox_offset=round-up(540+4+esize)`; write 540-byte header (magic `n+2` + `0D 0A 1A 0A`,
`dim[0]=6/7`, `dim[1..4]=1`, `dim[5]=rowLen`, `dim[6]=nCols`, datatype/bitpix per class,
`intent_code`+`intent_name` per the 3006/3002/3007 table), ext flag=1, esize/ecode/XML,
pad, then data transposed back to row-major. Target: byte/value-level round-trip
verifiable by `cifti_diff`-style comparison.

### GIFTI reader/writer (`read_gifti_native` / `write_gifti_native`)
**Validated bit-exact** (reference at `/tmp/test_gii_native.m`). Read: `fileread`,
`regexp(txt,'(?s)<DataArray .*?</DataArray>','match')` (the `(?s)` dotall flag is
required; avoid `\b`), pull attrs + `<Data>` payload, `matlab.net.base64decode`,
GZipBase64 inflate via `java.io.ByteArrayInputStream` →
`java.util.zip.InflaterInputStream` → `org.apache.commons.io.IOUtils.copy` (IOUtils ships
with MATLAB; **not** `GZIPInputStream`), `typecast`, `swapbytes` if endian mismatch,
RowMajor `reshape(vals,[Dim1 Dim0])'`. `POINTSET`→vertices, `TRIANGLE`→faces(+1),
`LABEL`→keys+LabelTable, `NONE/SHAPE/TIME_SERIES`→per-vertex data. Write: RowMajor byte
flatten of the transpose, zlib-compress JVM-side via `DeflaterOutputStream`,
`matlab.net.base64encode`, emit GIFTI XML. A no-Java path (ASCII/Base64Binary) is fully
native for writing; reading external GZip files still needs the JVM `Inflater`.

### Vendorable, permissively-licensed code
| Tool | License | Use |
|---|---|---|
| **cifti-matlab** core (`cifti_read`, `cifti_write`, `cifti_parse_xml`, `read_nifti2_hdr`) | **BSD-2-Clause** (take `read_nifti2_hdr` under its BSD option) | Model the hand-rolled CIFTI reader/writer on it. **Avoid** the `@xmltree` (LGPL) and `ft_cifti` (GPL) subdirs. |
| **gllmflndn/gifti** | **MIT** | Optional *fallback* only; ships per-platform MEX with no pure-MATLAB fallback, so prefer the validated native codec. |
| **CBIG RF warps** (already in repo) | **MIT** | vol↔surf mapping assets (§7). Keep `CBIG_LICENSE` alongside. |

`get_cifti_data.m` (diminfo model walk → cortex L/R + volumes) and
`extract_vol_from_cifti.m` (volumetric grayordinates → `fmri_data`) from
`canlabSurface/cifti_utils/` are ported in for the split and the `to_fmri_data` bridge.

---

## 7. Volume ↔ Surface mapping (native v1)

**Zero external binaries.** Reuse the **CBIG Registration-Fusion (RF-ANTs) warps already
vendored** at `CanlabCore/Cifti_plotting/CBIG_registration_fusion_surf2vol_vol2surf/`
(MIT; the `.mat` tables verified on disk). Replace FreeSurfer `MRIread/MRIwrite` with SPM12
`spm_vol/spm_read_vols/spm_write_vol` (already a hard CanlabCore dependency) and the MARS
kd-tree with the precomputed vertex-index grids. v1 standardizes the RF path on
**fsaverage** and the bundled barycentric path on **fs_LR-32k / fsavg-164k**.

### `vol2surf` (method on `fmri_data`/`image_vector` → `fmri_surface_data`, fsaverage)
Port `CBIG_RF_projectMNI2fsaverage.m`: (1) get the object's 3-D volume via
`reconstruct_image` and its 4×4 `mat` from `volInfo.mat`; (2) load
`lh/rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat` → `ras` (3×163842 MNI-RAS
per fsaverage7 vertex); (3) `vox = inv(mat(1:3,1:3))*(ras - mat(1:3,4))`,
`matcoord=[vox(2,:)+1; vox(1,:)+1; vox(3,:)+1]`; (4) per map column
`vertex_vals = interpn(vol3d, matcoord..., interp)` — `'linear'` for continuous,
`'nearest'` for label maps; OOB→0. Result populates `.dat` with an `fsavg_164k`
`brain_model`.

### `surf2vol` (method on `fmri_surface_data` → `fmri_data` + `.nii`)
Port `CBIG_RF_projectfsaverage2Vol_single.m`, **fast nearest gather path** (no kd-tree):
(1) load `..._avgMapping.vertex.mat` → `lh_vertex/rh_vertex` (256³ nearest-vertex-index
grids) and the liberal cortex mask NIfTI (defines grid + vox2ras); (2) vectorized gather
`out = zeros(numel(mask)); idx = lh_vertex>0; out(idx) = lh_vals(lh_vertex(idx));` (same
rh) — no loop/search; (3) reshape to the mask volume, wrap as `fmri_data` with the mask's
vox2ras (lands in MNI152). **The subcortical grayordinates need no mapping** — they are
already voxels+affine in the `volInfo` sub-block, so reconstruct them directly and add to
the cortical projection for a full-brain `fmri_data`. `fmri_data.write()` emits the
`.nii`. Optional linear path uses `prop.mat` + one `knnsearch` (Statistics Toolbox).

### Second native path (preferred when data is already fs_LR grayordinate)
`fsLR_32k ↔ MNI` via the repo's `resample_from_<space>_to_fsLR_32k(_nearestneighbor).mat`
barycentric structs through `render_on_surface`'s existing projection code — targets
fs_LR directly, no fsaverage↔fs_LR deformation needed.

### Deferred (documented limitation)
RF is a **fixed group-template MNI152↔fsaverage correspondence** — correct for group MNI
maps, **not** a per-subject native-space ribbon mapper. Best-in-class
ribbon-constrained per-subject mapping (needs subject white+pial) and fsaverage↔fs_LR HCP
deformation are deferred.

---

## 8. Reuse map

All paths under `/Users/f003vz1/Documents/GitHub/`.

| Existing file / method | How reused |
|---|---|
| `CanlabCore/@image_vector/{get_wh_image,descriptives,ica,mahal,pca,image_math,history,isempty,enforce_variable_types}.m` | **Inherit verbatim** by dispatch (touch only `.dat`). |
| `CanlabCore/@image_vector/{remove_empty,replace_empty}.m` | **Override as no-ops** (D5b) — not inherited; grayordinate `.dat` is never squeezed. |
| `CanlabCore/@image_vector/mean.m` (rebuild ~L141-143), `apply_mask.m` (pattern math), `threshold.m` (raw branch), `apply_parcellation.m` (aggregation core) | **Copy body, swap rebuild/space step** via `rebuild_like` / surface `resample_space`. |
| `CanlabCore/@image_vector/image_vector.m` (L252-279 property block; L358 `.nii`+exist branch); `CanlabCore/@fmri_data/fmri_data.m` (L313-351 analysis fields; L441-471 recast-from-sibling; L713 mandatory-volInfo) | **Mirror constructor polymorphism**, swap reader + extension test, relax volInfo check. |
| `CanlabCore/@fmri_data/{cat,horzcat,predict,ttest,regress,plot}.m` | **New by copying body** (geometry-agnostic cores; strip volume rebuild/resample) — NOT inherited by an `image_vector` subclass. |
| `canlabSurface/cifti_utils/get_cifti_data.m` | **Port** as the CIFTI→grayordinate splitter (cortex L/R + volumes). |
| `canlabSurface/cifti_utils/extract_vol_from_cifti.m` | **Port near-verbatim** for `to_fmri_data` (builds volInfo/voxlist/affine, wraps `fmri_data`). |
| `CanlabCore/Cifti_plotting/CBIG_registration_fusion_surf2vol_vol2surf/CBIG_RF_projectMNI2fsaverage.m`, `CBIG_RF_projectfsaverage2Vol_single.m` + `.mat` warps | **Port** (drop FreeSurfer→`spm_vol`, kd-tree→`vertex.mat`) for `vol2surf`/`surf2vol`. |
| `CanlabCore/@image_vector/render_on_surface.m` | **Reuse colormap/split-colormap/indexmap subfunctions**; feed per-vertex `.dat` straight to `FaceVertexCData` (skip `interp3`). |
| `CanlabCore/Cifti_plotting/surface_outlines.m` | **Reuse verbatim** (pure MATLAB). |
| `CanlabCore/Cifti_plotting/make_surface_figure.m` | **Reuse body** (swap `gifti()`→`load()`; fix hardcoded `/Users/sgeuter` addpath). |
| `CanlabCore/@fmridisplay/surface.m`, `CanlabCore/Visualization_functions/addbrain.m` | **Register** the object so `o2.surface{}` workflows keep working (`'hcp inflated'`, `'foursurfaces_hcp'`). |
| `CanlabCore/Cifti_plotting/cifti_struct_2_region_obj.m` | **Mirror field conventions** (XYZ/XYZmm/M/Z/voxlist+1) for the subcortical region path; it explicitly drops `'surf'` — the gap to fill. |
| `CanlabCore/@atlas/{num_regions,get_region_volumes,atlas2region}.m`, `condf2indic` | **Reuse / mirror** for surface atlas + parcellation (area replaces voxel volume). |
| `External/matlab_bgl` (or `graph/conncomp`) | **Reuse** for mesh connected components in `reparse_contiguous`/`region`. |
| `Neuroimaging_Pattern_Masks/.../2016_Glasser_*`, `2018_Schaefer_Yeo_*`, `2024_CANLab_atlas/` (`.label.gii`, `.dlabel.nii`) | **Test/atlas content** for surface parcellation import. |

**Known bugs to fix on reuse:** `make_surface_figure.m` hardcoded `/Users/sgeuter` addpath;
`render_cifti_on_brain.m` L188 extracts `CORTEX_LEFT` for the right hemisphere.

---

## 9. Standard mesh assets

CANlab already ships fs_LR-32k and fsaverage meshes as `.mat` (faces + vertices, loadable
straight into `patch()` with no `gifti` dependency) under
`CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/`.

| Asset | Vertex count / hemi | Source / license | Use |
|---|---|---|---|
| `S1200.{L,R}.midthickness_MSMAll.32k_fs_LR.surf.gii`, `S1200.{L,R}.sphere.32k_fs_LR.mat`, inflated `.mat` | 32,492 (29,696 L / 29,716 R non-medial) | **HCP S1200** Open-Access Data Use Terms (redistributable, not relicensable) — **already in repo** | Default fs_LR-32k display/geometry. |
| `fsavg_{inflated,white,sphere}_{lh,rh}.mat` | 163,842 | Already in repo | fsaverage-164k display + CBIG RF target. |
| `resample_from_{MNI152NLin2009cAsym,MNI152NLin6Asym,colin27}_to_{fsLR_32k,fsavg_164k}[_nearestneighbor].mat` | weights/vertices [Ntarget × 9] | Already in repo | MNI↔surface barycentric resampling (vol2surf path B, mesh↔mesh). |
| CBIG RF warps + `FSL_MNI152_FS4.5.0_cortex_estimate` mask | fsaverage7 = 163,842 | **MIT** (CBIG) — already in repo | `vol2surf`/`surf2vol` (§7). |
| **TemplateFlow `tpl-fsLR`** (32k/59k/164k midthickness/inflated/veryinflated/sphere + nomedialwall) | 32,492 / 59,292 / 163,842 | **CC0-1.0** | *Recommended add* for a clean-license default (pending Q3). |
| **TemplateFlow `tpl-fsaverage`** (white/pial/sphere + curv/sulc, 10k/41k/164k) | 2,562 / 40,962 / 163,842 | **CC0-1.0** | Optional fsaverage support with permissive license. |

Hemisphere convention (from HCP/CANlab code): **left=1, right=2**; surface handle `.Tag`
must contain `'left'`/`'right'` for projection; medial-wall mask value 1 = valid cortex.

---

## 10. Implementation phasing

Each milestone has a verification step using sample data + a unit test in
`CanlabCore/Unit_tests/`. Run `mcp__matlab__detect_matlab_toolboxes` on a clean path at
M1, M2 and M5 to prove zero external-toolbox runtime.

**M1 — Native I/O foundation (START HERE; highest value, no class yet). ✅ DONE 2026-06-25.**
Built (in `CanlabCore/Surface_tools/`): `canlab_read_cifti` / `canlab_write_cifti` (NIfTI-2
+ ecode-32 XML, hand-rolled regexp pull-parser; supports dscalar/dtseries/dlabel, surface +
voxel BrainModels, subcortical affine, scalar/series/label maps; writer does faithful
re-emit *and* full XML regeneration), `canlab_read_gifti` / `canlab_write_gifti` (.surf/
.func/.shape/.label; ASCII/Base64/GZipBase64), and shared `canlab_zlib_inflate` /
`canlab_zlib_deflate` JVM helpers. *Deliverable met:* self-contained, **proven to round-trip
exactly with cifti_read AND gifti removed from the path** (cifti cdata maxdiff=0, gifti
verts/faces maxdiff=0). *Verified:* `Unit_tests/surface_data/canlab_test_surface_io.m`
(5/5 pass) — shipped S1200 surf.gii (32,492 verts / 64,980 faces, bit-exact across all 3
encodings), GIFTI + CIFTI label tables, synthetic dscalar/dlabel grayordinate round-trips
(cdata + BrainModel start/count/vertlist/voxlist + subcortical affine preserved), and an
opportunistic real `.dscalar.nii` round-trip. Dev-time cross-validation against the
`cifti_read`/`gifti` oracles: bit-exact on dscalar/dtseries/dlabel/surf/label files.
Key gotchas solved: MATLAB regexp `\b` does NOT word-boundary here (use `\s`); CIFTI maps
can be self-closing `<MatrixIndicesMap .../>`; the affine element is
`TransformationMatrixVoxelIndicesIJKtoXYZ`; strip whitespace only for base64 (not ASCII)
payloads; inflate must run fully JVM-side.

*Naming note:* the design's `read_cifti_native` etc. were implemented under CANlab's
`canlab_<name>` standalone-function convention.

**M2 — Class skeleton + inheritance proof + import. ✅ DONE 2026-06-25.**
Built `@fmri_surface_data/fmri_surface_data.m` (`classdef fmri_surface_data < image_vector`)
with the full property set (`brain_model`/`geom`/`intent`/`series_info`/`label_table`/
`surface_space`/`mask`/X/Y/covariates/images_per_session/metadata_table/...). Constructor
polymorphism implemented: empty / cifti-struct / gifti-struct / plain-struct / `'key',value`
/ filename-autoload (dispatches `.gii`→`canlab_read_gifti`, `.nii`→`canlab_read_cifti`) /
`isa(image_vector)`-recast. `brain_model` + `.dat` are built directly from M1's reader output
(the reader already returns the `diminfo{1}.models` split, so the `get_cifti_data` port was
unnecessary). `surface_space`/`grayordinate_type` inferred from vertex counts. No-op
`remove_empty`/`replace_empty` overrides added (`@fmri_surface_data/{remove_empty,
replace_empty}.m`); `removed_voxels`/`removed_images` initialised all-false (D5b). *Verified:*
`Unit_tests/surface_data/canlab_test_surface_object_basic.m` (7/7 pass) — all 7 construction
paths (real dscalar/dtseries/dlabel, GIFTI label, `.surf.gii` geometry, struct, empty),
`.dat` stays full under `remove_empty`/`replace_empty` (no squeeze), `get_wh_image` subsets
maps only (rows preserved), `removed_*` vestigial, inherited method surface present.
`volInfo` sub-block population was **deferred to M3** (not needed for the M2 inheritance proof).

*M3 note (volInfo leak, Risk #1 confirmed):* inherited `descriptives` dereferences
`dat.volInfo.n_inmask` (its only volInfo use). Inherited methods that read `volInfo.*`
(`descriptives`, `montage`, `orthviews`, `flip`, `interpolate`, `get_xyzmm_coordinates`,
`rebuild_volinfo_from_dat`) need surface-aware guards/overrides in M3 — audit every
`volInfo.*` dereference and guard with `isfield` or override.

**M3 — Spatial overrides + interop. ✅ DONE 2026-06-26.**
Built `@fmri_surface_data/`: `compare_space.m` (0/1/2/3 contract, brain_model-based),
`reconstruct_image.m` (returns a struct: dense per-hemisphere vertex arrays with medial
wall = NaN, + a [X×Y×Z×maps] subcortical volume + its volInfo), `to_fmri_data.m` (exports
the subcortical voxel sub-block to an `fmri_data` in MNI), `write.m` (dispatches to
`canlab_write_cifti`/`canlab_write_gifti`; rebuilds the maps dimension from
intent+image_names/label_table/series_info; faithful re-emit of stashed source XML when the
layout is unchanged, else regenerate), and `rebuild_like.m` (re-wrap new data carrying the
geometry). Plus a private `build_volinfo_subblock.m` (0-based CIFTI → 1-based SPM affine
conversion; wh_inmask aligned with grayordinate rows) used by the constructor (populates the
inherited `volInfo` slot for the subcortical sub-block; empty for surface-only) and
`to_fmri_data`. *Verified:* `Unit_tests/surface_data/canlab_test_surface_space_recon.m`
(7/7) — volInfo sub-block + affine, `to_fmri_data` values vs subcortical rows,
`reconstruct_image` cortex NaN-medial-wall + subcortex voxel values, the full compare_space
contract (0/1/2/3), `rebuild_like` (geometry preserved, row-count enforced), and native
CIFTI dscalar/dlabel write→read round-trips (maxdiff 0, affine + map names + label table
preserved). Full surface suite (M1+M2+M3) = **19/19**.

*Deferred from M3:* (a) `reparse_contiguous` — needs the cortical mesh edge graph, which
depends on a bundled-mesh `geom` loader; moved to M5 (rendering) where that loader is built.
(b) Surface-aware overrides of the volInfo-dependent inherited QC/display methods
(`descriptives`, `montage`, `orthviews`, `flip`) — moved to M5/M6. `to_fmri_data` already
gives a clean volumetric escape hatch for those in the meantime.

**M4 — Volume ↔ surface mapping. ✅ DONE 2026-06-26.**
Built `@image_vector/vol2surf.m` (volume → fsaverage_164k `fmri_surface_data`: reconstruct
volume, sample the CBIG RF-ANTs `ras` MNI→fsaverage per-vertex coords with `interpn`;
`'interp'` linear/nearest) and `@fmri_surface_data/surf2vol.m` (fsaverage_164k → `fmri_data`
in MNI 2 mm: scatter the same `ras` coords with `accumarray` mean; `'reference'` to set the
target grid; writes `.nii` via the returned fmri_data). Self-consistent inverse pair, fully
native — no FreeSurfer/Workbench. Plus `Surface_tools/canlab_cbig_warp_path.m` (resolves the
vendored MIT warps) and the overrides `@fmri_surface_data/{mean,apply_mask,threshold}.m`
(mean across maps via `rebuild_like`; apply_mask = zero out-of-mask rows, no fmri_mask_image/
resample per D5b; threshold = raw-value only). *Verified:*
`Unit_tests/surface_data/canlab_test_surface_vol_map.m` (6/6) — vol2surf size/space, vertex
value == `interpn` sample (exact), **vol→surf→vol cortical correlation r = 0.9999**, surf2vol
rejects non-fsaverage objects, nearest-interp preserves integer labels, and mean/apply_mask/
threshold. Full surface suite (M1–M4) = **25/25**.

*Deferred from M4:* (a) the fs_LR-32k barycentric path (the `resample_from_*_to_fsLR_32k.mat`
structs are surface→surface fsnative→fs_LR deformations needing a vol→fsnative step first; the
fsaverage CBIG path is the v1, fsaverage↔fs_LR remains a later enhancement); (b)
cluster-extent `threshold` (needs the cortical mesh graph — M5); (c) combining the
subcortical sub-block into `surf2vol` output (use `to_fmri_data` for subcortex meanwhile).

**M5 — Rendering. ✅ DONE 2026-06-26.**
Built `@fmri_surface_data/surface.m` (native 4-panel L/R×lateral/medial render on the bundled
inflated mesh; `'surftype'`/`'which_image'`/`'clim'`/`'pos_colormap'`/`'neg_colormap'`;
`'existingsurface'` and `'mni_surface'` modes), `@fmri_surface_data/render_on_surface.m`
(override: colors given patch handles — DIRECT per-vertex truecolor when the patch matches a
hemisphere's vertex count, e.g. any addbrain fs_LR/fsaverage mesh; else projects to a volume
via `obj_to_volume` and reuses `@image_vector/render_on_surface`), `@fmri_surface_data/plot.m`
(QC: histogram, per-map mean±sd, coverage, mean-map render), and the deferred
`@fmri_surface_data/reparse_contiguous.m` (cortex = mesh edge-graph `conncomp`; subcortex =
`spm_clusters`; writes `brain_model.cluster`). Helpers: `Surface_tools/canlab_surface_vertexcolors.m`
(value→truecolor split colormap, medial wall/zero = gray) and the private
`load_surface_geom.m` (bundled S12000 inflated / S1200 midthickness/sphere / fsaverage
inflated meshes) and `obj_to_volume.m`. **Native-space rendering needs no resampling** —
addbrain's `hcp inflated` is the fs_LR-32k template (32492) and `inflated`/`fsavg` is fsaverage
(163842), exact matches. **Rendering on any other (e.g. MNI pial) surface** works via the
volume projection. *Verified:* `Unit_tests/surface_data/canlab_test_surface_render.m` (7/7) —
native 4-panel render (truecolor per vertex, correct vertex counts for both spaces), medial
wall renders gray, direct coloring of an addbrain native patch, via-volume render onto an
addbrain MNI surface, `reparse_contiguous` (every active grayordinate labeled, inactive = 0),
and `plot`. Visually confirmed (transcriptomic gradient on inflated fs_LR; smooth map on
fsaverage and on an MNI surface). Full surface suite (M1–M5) = **32/32**.

*Deferred from M5:* deep `@fmridisplay`/`addbrain` registration (the standalone `surface`/
`render_on_surface` cover the rendering need; fmridisplay montage integration is optional
polish); `surface_outlines` overlay and a colorbar/legend on the native path; surface-aware
overrides of `montage`/`orthviews`/`flip`/`descriptives` (still routed via `to_fmri_data`).

**M6 — Analysis parity. ✅ DONE 2026-06-26.**
Built `@fmri_surface_data/`: `cat.m` + `horzcat.m` (concatenate along maps; `compare_space==0`
required; per-map fields X/Y/covariates/image_names/metadata_table concatenated),
`ttest.m` and `predict.m` and `ica.m` (**delegated** to the corresponding `@fmri_data`/
`@image_vector` methods via a private `as_fmri_data_proxy.m` — wraps the object as an
`fmri_data` with a dummy 1-D volInfo treating each grayordinate as a voxel — then remaps
geometry-bearing results back with `rebuild_like`; this reuses the entire battle-tested
algorithm/CV machinery), and a native lightweight `regress.m` (per-grayordinate OLS with
`obj.X`; betas in `.dat`, t/p/se/dfe in `additional_info.statistic`). `ttest`/`regress`
return an `fmri_surface_data` carrying the stat in `.dat` (t / betas) plus parallel fields in
`additional_info.statistic` (full surface statistic_image subclass deferred to M7). `predict`
returns `[cverr, stats, optout]` with `stats.weight_obj` remapped to a surface object.
*Verified:* `Unit_tests/surface_data/canlab_test_surface_analysis.m` (4 pass + ica filtered) —
cat/horzcat (+ space-mismatch error), ttest vs MATLAB `ttest` (matches to single precision),
regress OLS vs `\` (and recovered slope), predict CV round-trip with remapped weight map.
Full surface suite (M1–M6) = **36/36** (+ ica skipped where its toolbox is absent).

*Upstream fix made:* `@image_vector/ica.m` had a stray, never-functional debug line
(`figure; plot(B(:), W(:), …)` referencing undefined `B`/`W`) that made `ica` error on *every*
call for all `image_vector` subclasses — removed it. `ica` still requires `icatb_fastICA`
(GIFT/icatb toolbox); it is the one `fmri_surface_data` method that is not fully self-contained.

*Deferred from M6:* full `glm_map` integration for `regress` (contrasts/diagnostics — use
`to_fmri_data` + `fmri_data.regress`/`glm_map` for now); a dedicated
`fmri_surface_statistic_image` subclass with native cluster-extent `threshold` (M7).

**M7 — Parcellation / region + surface atlas. ✅ DONE 2026-06-26.**
Built `@fmri_surface_data/`: `apply_parcellation.m` (parcel means `[nMaps × nParcels]` from a
`.dlabel` surface object or an integer key vector; key 0/NaN = background/medial wall excluded;
optional `'area'` weighting via per-vertex mesh surface area; returns labels + summary table),
`surface_region.m` (contiguous clusters → struct array with `.struct`/`.type`/`.grayord_rows`/
`.vertex_indices`/`.XYZmm` centroid/`.numVox`/`.val`), and cluster-extent thresholding added to
`threshold.m` (`'k', N` via `reparse_contiguous`). Surface atlases are just `.dlabel`
`fmri_surface_data` objects (label keys + label_table), so they load through the standard
constructor — no separate import path needed. *Verified:*
`Unit_tests/surface_data/canlab_test_surface_parcellation.m` (5/5) — parcel means == keys on a
known synthetic parcellation and the real Gordon333+Tian atlas, label_table names, area
weighting finite, `compare_space` enforcement, cluster-extent threshold (no surviving cluster
< k), and `surface_region` partitioning the active grayordinates.

*Deferred from M7:* the dedicated `fmri_surface_statistic_image` / `fmri_surface_atlas`
subclasses (the single class + `additional_info.statistic` + `.dlabel` objects cover the
functionality; subclasses are polish for `get_wh_image`-synced stat fields).

**M8 — Docs + deprecation. ✅ DONE 2026-06-26.**
Authored `docs/fmri_surface_data_methods.md` (full method/option reference, updated through M7)
and `docs/fmri_surface_data_walkthrough.m` (runnable cell-mode walkthrough incl. parcellation +
group analysis; verified end-to-end). Added deprecation pointers in the external-dependent
`Cifti_plotting/render_cifti_on_brain.m` and `plot_surface_map.m` headers steering new code to
the native `fmri_surface_data`. *Deferred (optional, asset fetch — not code):* bundling CC0
TemplateFlow `tpl-fsLR`/`tpl-fsaverage` meshes for a relicensable default (the in-repo HCP/
FreeSurfer meshes work today; this is license hygiene); a `validate_object`-style suite (the
7-file, 41-test surface suite covers M1–M7).

---

## Status: v1 complete (M1–M8)

All eight milestones delivered. `fmri_surface_data` reads/writes CIFTI-2 + GIFTI natively,
holds grayordinate data with full `fmri_data`-style method parity (construct, reconstruct,
to_fmri_data, write, compare_space, cat/horzcat, mean/apply_mask/threshold, ttest/regress/
predict/ica, vol2surf/surf2vol, surface/render_on_surface/plot, reparse_contiguous/
apply_parcellation/surface_region), and runs with **no external toolbox** (sole exception:
`ica`, which needs the GIFT/icatb toolbox like the base class). Surface suite: **41 tests
passing** across 7 files. Best-in-class per-subject nonlinear mapping and fsaverage↔fs_LR
deformation remain the main future enhancements.

---

## 11. Risks & open questions

**Risks**
1. **`volInfo` dual-role leak.** Inherited methods read `volInfo.mat`/`.dim`/`.wh_inmask`
   (counts: mat 25, wh_inmask 63, xyzlist 41, image_indx 37, cluster 22). Mitigation:
   the surface truth lives in `brain_model`; `volInfo` describes *only* the volume
   sub-block, and is empty for surface-only objects. Audit every `volInfo.*` dereference
   (highest-risk overlooked inheritors: `rebuild_volinfo_from_dat`,
   `get_xyzmm_coordinates`, `extract_roi_averages`, `montage`, `orthviews`, `flip`,
   `interpolate`) and override or guard with `isfield`; do NOT rely on inheritance for any
   method touching those fields.
2. ~~**remove/replace_empty contract**~~ **ELIMINATED (D5b).** `.dat` is always the full
   grayordinate set; `remove_empty`/`replace_empty` are no-ops and `removed_voxels` is a
   vestigial all-false vector. This whole risk class (the "forgot to `replace_empty`" bugs)
   is removed by design. Remaining care: keep `.dat` rows in 1:1 order with `brain_model`.
3. **Native parser correctness.** CIFTI XML edge cases (CDATA, BrainModels tiling with no
   gaps, row-major permute, `scl_slope`); GIFTI GZip = **raw zlib not gzip**, apache base64
   corruption, in-place JVM inflate garbage. Mitigation: validate against real HCP files;
   keep `/tmp/test_gii_native.m` as the oracle; carry the three GIFTI gotchas verbatim.
4. **`compare_space` 4-valued contract** must be preserved (not boolean) or `cat`/
   `apply_mask` break subtly.
5. **`cat`/`horzcat`/`plot`/`predict`/`ttest`/`regress` are NOT inherited** (only in
   `@fmri_data`) — forgetting to provide them yields silent "method not found".
6. **`mean`/rebuilders hardcode** the `image_vector`+`fmri_data` re-wrap — must route
   through `rebuild_like` or they silently emit volume objects from surface data.
7. **Space proliferation.** Native CIFTI yields fs_LR-32k; CBIG RF yields fsaverage —
   different topologies. fsaverage↔fs_LR deformation is deferred, so the two cannot be
   `cat`'d without resampling. Gate `cat`/`compare_space` on `surface_space`.
8. **RF is a fixed group correspondence** — wrong for per-subject native space; document
   clearly. Exclude medial wall from `apply_parcellation` column-normalization.
9. **Statistic/atlas parallel fields** (`.p/.ste/.sig`; integer keys) need to be carried
   through `get_wh_image` (column subsetting) or they desync from `.dat` (handled in M7
   subclasses). Simpler than `fmri_data`, since `remove_empty`/`replace_empty` are no-ops.
10. **Surface area** for `get_region_volumes`/`rmsv` cannot be a no-op stub — compute from
    face geometry early (M7).
11. **Known reuse bugs:** `make_surface_figure.m` hardcoded addpath;
    `render_cifti_on_brain.m` L188 CORTEX_LEFT-for-right copy-paste — fix on reuse.
12. **License hygiene (hard requirement):** must NOT vendor `@xmltree` (LGPL) / `ft_cifti`
    (GPL) / `gifti` MEX. Hand-roll permissively; ship CBIG MIT + HCP/TemplateFlow terms;
    verify before bundling.

**Open questions (need user confirmation — see §3):**
- **Q1.** Class name: `fmri_surface_data` (prefixed) vs. house-style `fmri_surface_data`?
- **Q2.** Confirm fs_LR-32k (91k) as the canonical default space.
- **Q3.** Add CC0 TemplateFlow `tpl-fsLR`/`tpl-fsaverage` meshes for a clean-license
  default, or rely on the HCP S1200 surfaces already shipped?
- **Q4.** Hand-roll the CIFTI reader/writer (recommended, fully self-contained) vs. vendor
  the BSD-2 cifti-matlab core?
