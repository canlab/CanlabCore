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

## Surface spaces

Every object carries a `surface_space` tag naming the cortical mesh its vertices
live on. The supported spaces, their resolutions, and their sources:

| `surface_space` | Verts/hemi | What it is | Source |
|---|---|---|---|
| `fsaverage_164k` | 163,842 | **FreeSurfer fsaverage** average-subject spherical surface, 7th-order icosahedron. The output of `vol2surf`. | Fischl et al. 1999 |
| `fsaverage6` | 40,962 | fsaverage 6th-order icosahedron — a **nested subset** (first N vertices) of fsaverage-164k. | Fischl et al. 1999 |
| `fsaverage5` | 10,242 | fsaverage 5th-order icosahedron (nested subset). | Fischl et al. 1999 |
| `fsaverage4` | 2,562 | fsaverage 4th-order icosahedron (nested subset). | Fischl et al. 1999 |
| `fsLR_32k` | 32,492 | **HCP fs_LR-32k** — the left/right-symmetric grayordinate cortical mesh used by CIFTI `.dscalar`/`.dtseries`/`.dlabel` files. | Van Essen et al. 2012 |
| `onavg_41k` | 40,962 | **onavg** (OpenNeuro Average) equal-area template, den-41k — vertices placed so cortical regions are sampled at uniform resolution. Resample-only (via its `space-fsLR` registration sphere). | Feilong et al. 2024 (CC0) |
| `onavg_10k` | 10,242 | onavg equal-area template, den-10k. | Feilong et al. 2024 (CC0) |

**Grayordinate objects (the "91k" CIFTI model).** A native CIFTI file combines
the `fsLR_32k` cortex (59,412 vertices = 32,492 × 2 minus the medial wall) with
**31,870 subcortical/cerebellar voxels**, for 91,282 grayordinates — cortex on the
surface, subcortex in a volume. The cortex uses `fsLR_32k`; the subcortical block
is a standard MNI152 volume (extract it with `to_fmri_data`). See Glasser et al.
2013 for the grayordinate model.

**Moving between spaces.** `resample_surface(obj, target_space)` resamples the
cortical data between any of these spaces natively (see the [Surface resampling](#surface-resampling-and-volume--surface-mapping)
section); `resample_surface(obj, 'list')` prints the keyword list. `vol2surf` /
`surf2vol` convert between an MNI152 **volume** and `fsaverage_164k`.

*References.* Fischl B, Sereno MI, Tootell RBH, Dale AM (1999), *Hum Brain Mapp*
8(4):272–284. Van Essen DC, Glasser MF, Dierker DL, Harwell J, Coalson T (2012),
*Cereb Cortex* 22(10):2241–2262. Glasser MF et al. (2013), *NeuroImage*
80:105–124. Feilong M, Guo J, Gobbini MI, Haxby JV (2024), *Nat Methods*
21:2069–2078 (onavg; TemplateFlow `tpl-onavg`, CC0).

## Properties

`fmri_surface_data` inherits the `image_vector` properties (`dat`, `volInfo`,
`removed_voxels`, `removed_images`, `image_names`, `fullpath`, `files_exist`,
`history`, `dat_descrip`, `source_notes`) and adds:

| Property | Description |
|---|---|
| `brain_model` | Geometry source of truth (mirrors CIFTI BrainModels). `.models{i}`: `.struct`, `.type` (`surf`/`vox`), `.start`, `.count`, `.numvert`, `.vertlist` (0-based), `.voxlist`; plus `.vol` (`.dims`, `.sform`), `.grayordinate_type`, `.cluster`. |
| `geom` | **Attached CUSTOM mesh geometry only** (`.vertices`/`.faces`), set when the object is built from a `.surf.gii` that carries geometry. It is **empty for a data-only CIFTI** (e.g. a `.dscalar`) — that is expected, not a bug. Standard-space meshes are **not** stored here: for a known `surface_space` (fs_LR-32k, fsaverage-164k, …) `surface()` and resampling load the mesh on demand from the bundled assets (`private/load_surface_geom`, `canlab_canonical_brains/Canonical_brains_surfaces`). Use `.geom` only to carry a non-standard mesh no `surface_space` keyword can reproduce. |
| `imagetype` | `dscalar` / `dtseries` / `dlabel` / `func` / `shape` / `label`. |
| `series_info` | For `.dtseries`: `.start`/`.step`/`.unit`/`.exponent`. |
| `label_table` | For `.dlabel`/`.label`: a MATLAB table with variables `key`, `name`, `rgba` (Nx4). |
| `surface_space` | The cortical mesh the vertices live on — one of `fsaverage_164k`, `fsaverage6`, `fsaverage5`, `fsaverage4`, `fsLR_32k` (see [Surface spaces](#surface-spaces)). Gatekeeps `compare_space`; drives mesh/warp choice and `resample_surface`. |
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

## Surface resampling and volume ↔ surface mapping

| Method | From | One-liner |
|---|---|---|
| `resample_surface` | `@fmri_surface_data` | Resample the cortex to another **surface** space — `fsaverage_164k` ↔ `fsLR_32k` (HCP CIFTI) and nested `fsaverage6`/`5`/`4`. Barycentric by default (weights built once, reused across maps), nearest for binary/label maps; subcortex carried through. `resample_surface(obj,'list')` prints the spaces |
| `vol2surf` | `@image_vector` | Project a volumetric image (MNI152) onto the fsaverage-164k surface. **Is the CBIG RF-ANTs mapper natively** (Wu et al. 2018) — a reimplementation of `CBIG_RF_projectMNI2fsaverage` (`interpn`), no FreeSurfer needed |
| `surf2vol` | `@fmri_surface_data` | Project an fsaverage-164k object back to an MNI152 `fmri_data` volume — native inverse using the same CBIG RF-ANTs warp (`accumarray` scatter) |

**`resample_surface(obj, target_space[, 'interp', 'nearest'])`.** Target spaces:
`fsaverage_164k` (aliases `fsaverage`, `fsavg`, `164k`), `fsaverage6`/`5`/`4`,
`fs_LR_32k` (aliases `fsLR`, `hcp`, `32k`). fsaverage↔fs_LR uses the vendored HCP
registration ("deformed") sphere so both meshes share one spherical frame and are
resampled with barycentric (or nearest) interpolation; fsaverage down-sampling is
the **exact nested icosahedral subset**. Barycentric weights depend only on
geometry, so they are computed once and applied to every map as a sparse
matrix-multiply (a 50-map object costs about the same as one map). Continuous data
defaults to barycentric (`'linear'`); binary masks and `.dlabel` images default to
`'nearest'`.

```matlab
s   = vol2surf(ttest(load_image_set('emotionreg')));  % fsaverage_164k
s32 = resample_surface(s, 'fsLR_32k');                % -> HCP fs_LR-32k
s6  = resample_surface(s, 'fsaverage6');              % -> nested fsaverage6
lab = resample_surface(atlas_surf, 'fsLR_32k', 'interp', 'nearest');  % labels
```

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
| `surface` | `@fmri_surface_data` | Native mode returns a **managed `fmridisplay`** (4-panel, matching meshes painted directly); also `'existingsurface'` handles or an `'mni_surface'` addbrain surface (both return a struct) |
| `render_on_surface` | `@fmri_surface_data` | Color existing patch handles: directly if the mesh matches, else via a volume projection |
| `plot` | `@fmri_surface_data` | QC panel: value histogram, per-map mean±sd, coverage, mean-map render |

`surface` / `render_on_surface` accept the **same color vocabulary as the volume
visualization pipeline** (`addblobs` / `set_colormap`): `clim` / `cmaprange`,
`colormap` / `colormapname` (a single sequential map), `pos_colormap` /
`neg_colormap`, `splitcolor`, `maxcolor` / `mincolor`, and `color` (solid). So the
same options color surface data and volume blobs.

### Managed display (`fmridisplay`)

A surface object lives in a stateful `fmridisplay` as a **surface-native layer**, so
the same managed-display experience works for surface and volume data. The simplest
entry point is `surface(obj)` itself — it now **returns a managed `fmridisplay`**
with the matching native surfaces added and the data painted:

```matlab
o2 = surface(surf_stat);                           % managed display, native meshes
o2 = set_colormap(o2, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);  % recolors in place
o2 = removeblobs(o2);                              % restores the anatomy

% Or add to a display you already have:
o2 = fmridisplay;
o2 = surface(o2, 'hcp inflated left');  o2 = surface(o2, 'hcp inflated right');
o2 = addblobs(o2, surf_stat, 'colormap', 'hot');   % paints the fs_LR meshes directly
% equivalently, add + paint the object's own matching surfaces in one call:
o2 = surface(o2, surf_stat);
```

`addblobs` (and `surface(o2, obj)`) detect the `fmri_surface_data` and paint matching
cortical meshes **directly from the per-vertex data** (no volume resampling), using
the same central `canlab_colormap` value→color map as montages, so colors match. Each
hemisphere shows its own data (resolved by tag, falling back to vertex x-position when
`addbrain` relabels the `foursurfaces` patches). The layer participates in
`set_colormap` / `set_opacity` / `rethreshold` / `removeblobs` / the controller like a
volume layer (`rethreshold` on a surface layer is a magnitude cutoff applied per
vertex — no volume/`volInfo` needed). It is **skipped with a clear `spacemismatch`
warning** on any registered surface whose mesh is a different surface space than the
object (add the matching surfaces instead — `surface(obj)` does so automatically).

**Mixed grayordinate objects (cortex + subcortex).** For an object that also has
subcortical voxels (a `91k`-style dscalar), `surface(obj)` adds a *second* managed
layer: the subcortical grayordinates go in as a standard **volume layer**, which paints
the subcortical anatomical meshes the `foursurfaces` set draws (thalamus, brainstem,
cerebellum, …) and is confined to them (it does not bleed onto the cortical meshes /
medial wall). Both layers respond to `set_colormap` / `rethreshold` / `set_opacity` /
the controller together. `montage(o2, obj)` builds the slice montage and renders those
same subcortical grayordinates on the slices (cortical surface data has no volume and
is not shown on slices — use `surface(obj)` for the cortex).

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
