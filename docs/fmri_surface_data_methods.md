# `fmri_surface_data` methods, organized by area

This is a functional index of methods available on an `fmri_surface_data`
object — the CANlab object for **cortical-surface and grayordinate (HCP CIFTI)
data**. It is the surface analogue of `fmri_data`: data are stored flat as a
`[grayordinates × maps]` matrix in `.dat`, with a `brain_model` describing the
cortical-surface vertices and subcortical voxels, so the familiar CANlab method
names (`mean`, `threshold`, `apply_mask`, `surface`, `write`, …) work on surface
data and interoperate with the rest of the toolbox.

`fmri_surface_data` is a subclass of `image_vector`, so it inherits the
purely-`.dat` methods listed in
[image_vector_methods.md](image_vector_methods.md) (`get_wh_image`, `mahal`,
`pca`, `image_math`, …). Only methods owned by `fmri_surface_data` (defined in
`@fmri_surface_data/`) are listed below. Type `methods(my_obj)` for the live list.

**Everything runs natively in MATLAB — no gifti toolbox, FieldTrip, Connectome
Workbench, or FreeSurfer is required at runtime** (sole exception: `ica`, which
needs the GIFT/`icatb` toolbox, like the base `image_vector` method).

- **Tutorial:** [Surface & grayordinate data how-to](workflows/fmri_surface_data_howto.md)
- **Runnable script:** [`CanlabCore/docs/fmri_surface_data_walkthrough.m`](../CanlabCore/docs/fmri_surface_data_walkthrough.m)
- **Design rationale / roadmap:** [`CanlabCore/docs/fmri_surface_data_design_plan.md`](../CanlabCore/docs/fmri_surface_data_design_plan.md)

## Concept

The HCP "grayordinate" model represents the cortex as *surface vertices* and the
subcortex/cerebellum as *voxels*. Because that is already a flat
`[grayordinates × maps]` matrix, it maps directly onto CANlab's `.dat`. Two
differences from `fmri_data`:

- **`brain_model` replaces `volInfo`** as the geometry source of truth (cortical
  vertices have no voxel coordinates). The inherited `volInfo` is populated to
  describe only the subcortical voxel sub-block (empty for surface-only objects).
- **No empty-squeezing:** grayordinate data is already compact, so `.dat` is
  *always* the full set; `remove_empty` / `replace_empty` are no-ops.

## Properties

`fmri_surface_data` inherits the `image_vector` properties (`dat`, `volInfo`,
`removed_voxels`, `removed_images`, `image_names`, `fullpath`, `files_exist`,
`history`, `dat_descrip`, `source_notes`) and adds:

| Property | Description |
|---|---|
| `brain_model` | Geometry source of truth (mirrors CIFTI BrainModels). `.models{i}`: `.struct`, `.type` (`surf`/`vox`), `.start`, `.count`, `.numvert`, `.vertlist` (0-based), `.voxlist`; plus `.vol` (`.dims`, `.sform`), `.grayordinate_type`, `.cluster`. |
| `geom` | Cortical mesh cache (faces/vertices) used for rendering. |
| `intent` | `dscalar` / `dtseries` / `dlabel` / `func` / `shape` / `label`. |
| `series_info` | For `.dtseries`: `.start`/`.step`/`.unit`/`.exponent`. |
| `label_table` | For `.dlabel`/`.label`: struct array `.key`/`.name`/`.rgba`. |
| `surface_space` | e.g. `fsLR_32k`, `fsaverage_164k`. Gatekeeps `compare_space`; drives mesh/warp choice. |
| `mask` | Optional `[nGray × 1]` logical (or same-space object) for `apply_mask`. |
| `X` / `Y` / `covariates` / `images_per_session` / `metadata_table` / … | Per-map annotations, same names/roles as `fmri_data`. |
| `additional_info` | Free-form struct; also holds `.statistic` (from `ttest`/`regress`) and the source CIFTI xml/hdr. |

## Construction and import/export

| Method | From | One-liner |
|---|---|---|
| `fmri_surface_data` | `@fmri_surface_data` | Constructor: from a CIFTI/GIFTI file, a reader struct, `'key',value` pairs, or an `image_vector` recast |
| `write` | `@fmri_surface_data` | Write natively to CIFTI-2 (`.nii`) or GIFTI (`.gii`); rebuilds map metadata; faithful re-emit or regenerate |
| `to_fmri_data` | `@fmri_surface_data` | Export the subcortical voxel grayordinates as a volumetric `fmri_data` (writeable to NIfTI) |

## Geometry and space

| Method | From | One-liner |
|---|---|---|
| `reconstruct_image` | `@fmri_surface_data` | Unpack `.dat` into dense per-hemisphere vertex arrays (medial wall NaN) + a subcortical volume |
| `compare_space` | `@fmri_surface_data` | Compare grayordinate spaces (0 same / 1 diff space / 2 missing / 3 diff in-data), preserving the contract |
| `reparse_contiguous` | `@fmri_surface_data` | Label contiguous clusters (cortex = mesh edge graph; subcortex = 26-connectivity) into `brain_model.cluster` |
| `rebuild_like` | `@fmri_surface_data` | Wrap new `[nGray × K]` data into an object carrying this geometry (used internally) |

## Volume ↔ surface mapping

| Method | From | One-liner |
|---|---|---|
| `vol2surf` | `@image_vector` | Project a volumetric image (MNI152) onto the fsaverage-164k surface via the CBIG RF mapping (`interpn`) |
| `surf2vol` | `@fmri_surface_data` | Project an fsaverage-164k object back to an MNI152 `fmri_data` volume (self-consistent inverse of `vol2surf`) |

## Data operations

| Method | From | One-liner |
|---|---|---|
| `mean` | `@fmri_surface_data` | Average across maps → single-map object |
| `apply_mask` | `@fmri_surface_data` | Zero out-of-mask grayordinates (no `fmri_mask_image`, no resample; `.dat` stays full) |
| `threshold` | `@fmri_surface_data` | Raw-value threshold; optional cluster-extent (`'k', N`) via `reparse_contiguous` |
| `cat` / `horzcat` | `@fmri_surface_data` | Concatenate objects along maps (same space required); also `[a, b, …]` |
| `remove_empty` / `replace_empty` | `@fmri_surface_data` | No-ops (grayordinate data is already compact; overrides so inherited callers compose) |

## Analysis

`ttest`, `predict`, and `ica` delegate to the corresponding `fmri_data` /
`image_vector` methods — treating each grayordinate as a feature — and remap
geometry-bearing results back to the surface, so those algorithms are identical
to their volumetric counterparts. `regress` is a native grayordinate-wise OLS
implementation.

| Method | From | One-liner |
|---|---|---|
| `ttest` | `@fmri_surface_data` | Grayordinate-wise one-sample t-test across maps (t-map in `.dat`; p/ste/sig in `additional_info.statistic`) |
| `regress` | `@fmri_surface_data` | Grayordinate-wise OLS onto `obj.X` (betas in `.dat`; t/p/se in `additional_info.statistic`) |
| `predict` | `@fmri_surface_data` | Cross-validated multivariate prediction (`obj.Y`); weight map remapped to a surface object |
| `ica` | `@fmri_surface_data` | Spatial ICA (requires the GIFT/`icatb` toolbox) |

## Rendering and QC

| Method | From | One-liner |
|---|---|---|
| `surface` | `@fmri_surface_data` | Render on the native mesh (4-panel), on `'existingsurface'` handles, or an `'mni_surface'` addbrain surface |
| `render_on_surface` | `@fmri_surface_data` | Color existing patch handles: directly if the mesh matches, else via a volume projection |
| `plot` | `@fmri_surface_data` | QC panel: value histogram, per-map mean±sd, coverage, mean-map render |

## Parcellation and regions

| Method | From | One-liner |
|---|---|---|
| `apply_parcellation` | `@fmri_surface_data` | Parcel means `[nMaps × nParcels]` from a `.dlabel` object or key vector (optional `'area'` weighting) |
| `surface_region` | `@fmri_surface_data` | Summarize contiguous clusters as region-like structs (`.struct`, `.XYZmm`, `.numVox`, `.val`, …) |

## Volume-only methods (redirected or masked)

`fmri_surface_data` is a full subclass of `image_vector`, but many `image_vector`
methods assume a single 3-D volume and are not meaningful for surface/grayordinate
data. These are handled so you get useful behavior or a clear message instead of a
cryptic error, and so `methods(obj)` stays clean:

- **Redirected to the subcortical volume:** `orthviews`, `montage`, `slices` route
  the subcortical grayordinates through `to_fmri_data` and call the `fmri_data`
  method (so you see the subcortex on slices). For a cortex-only object they give a
  clear error pointing to `surface(obj)`.
- **Masked (hidden + informative error):** volume-only methods with no surface
  meaning — e.g. `flip`, `isosurface`, `interpolate`, `resample_space`,
  `extract_gray_white_csf`, `searchlight`, `slice_movie`, `trim_mask`,
  `read_from_file`, `extract_roi_averages`, … — are overridden in a `methods
  (Hidden)` block. They no longer appear in `methods(obj)` / tab-completion; calling
  one raises `fmri_surface_data:unsupportedMethod` with a pointer to the surface
  equivalent (e.g. `apply_parcellation` instead of `extract_roi_averages`). The
  object remains an `image_vector` subclass — only these specific methods are masked.

## Surface data in Neuroimaging_Pattern_Masks

CanlabCore ships one small continuous sample map
(`Sample_datasets/CIFTI_surface_examples/S1200.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii`).
**Surface atlases and additional CIFTI maps live in the companion
`Neuroimaging_Pattern_Masks` (NPM) repository**, a standard CANlab dependency
(add it with `canlab_toolbox_setup`). Once on the path, load any of them by
filename, e.g. `fmri_surface_data(which('transcriptomic_gradients.dscalar.nii'))`.
Summary of what's available (representative files per family):

| Atlas / map family | Type | Space | Notes |
|---|---|---|---|
| **Gordon-333 + Tian S1–S4** (`Gordon333.32k_fs_LR_Tian_Subcortex_S*.dlabel.nii`) | `.dlabel` | fs_LR-32k, 91k | Gordon cortical parcels + Tian subcortex (4 subcortical scales) |
| **Schaefer-400 + Tian S1–S4** (`Schaefer2018_400Parcels_{7,17}Networks_..._S*.dlabel.nii`) | `.dlabel` | fs_LR-32k, 91k | Schaefer 400 (7 or 17 networks) + Tian subcortex |
| **HCP-MMP (Glasser) + Tian S1–S4** (`Q1-Q6_...CorticalAreas..._Tian_Subcortex_S*.dlabel.nii`) | `.dlabel` | fs_LR-32k, 91k | 360 HCP-MMP areas + Tian subcortex |
| **Glasser-2016 (HCP-MMP 360)** (`Glasser_2016.32k.{L,R}.label.gii`) | `.label.gii` | fs_LR-32k | Cortex only, per hemisphere |
| **JulichBrain 3.0.3** (`{lh,rh}.JulichBrainAtlas_3.0.3[.fslr_32k].label.gii`) | `.label.gii` | fsaverage & fs_LR-32k | Cytoarchitectonic cortex, per hemisphere |
| **CANLab2024 (coarse)** (`openCANLab2024_..._coarse_2mm.dlabel.nii`, `CANLab2024_fsaverage-164k_hemi-{L,R}_coarse.label.gii`) | `.dlabel` / `.label.gii` | MNI152NLin6Asym 91k / fsaverage-164k | CANLab whole-brain atlas |
| **CIT168 amygdala nuclei** (`CIT168_AmyNuc.dlabel.nii`) | `.dlabel` | subcortical | Amygdala subnuclei |
| **Transcriptomic gradients** (`transcriptomic_gradients.dscalar.nii`) | `.dscalar` | fs_LR-32k, 91k | Gene-expression gradient maps (continuous) |
| **HCP group ICAs** (`hcp_d{15,25,50}_ICs.dscalar.nii`) | `.dscalar` | fs_LR-32k, 91k | Resting-state group ICA component maps |
| **HCP spectral bases** (`spectral_bases_200.dscalar.nii`) | `.dscalar` | fs_LR-32k, 91k | Spatial-basis functions |

(Grayordinate parcellations = cortex surface + subcortical volume; `.label.gii`
files are cortex-only, per hemisphere. Paths are under
`Neuroimaging_Pattern_Masks/Atlases_and_parcellations/` and
`.../spatial_basis_functions/`.)

## Native I/O functions (`CanlabCore/Surface_tools/`)

Standalone functions backing the object; usable directly, no external toolbox.

| Function | Purpose |
|---|---|
| `canlab_read_cifti` / `canlab_write_cifti` | Native CIFTI-2 reader / writer |
| `canlab_read_gifti` / `canlab_write_gifti` | Native GIFTI reader / writer |
| `canlab_surface_vertexcolors` | Map values → truecolor (split hot/cool; NaN/zero = gray) |
| `canlab_cbig_warp_path` | Resolve the vendored CBIG registration-fusion warp paths |
