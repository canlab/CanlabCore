# `fmri_surface_data` — methods and usage reference

`fmri_surface_data` is the CANlab object for **cortical-surface and grayordinate
(HCP CIFTI) data**. It is the surface analogue of `fmri_data`: data are stored
flat as a `[grayordinates × maps]` matrix in `.dat`, with a `brain_model`
describing the cortical-surface vertices and subcortical voxels, so the familiar
CANlab method names (`mean`, `threshold`, `apply_mask`, `surface`, `write`, …)
work on surface data and interoperate with the rest of the toolbox.

Everything runs **natively in MATLAB** — no gifti toolbox, FieldTrip, Connectome
Workbench, or FreeSurfer is required at runtime.

- **Walkthrough (runnable):** [`fmri_surface_data_walkthrough.m`](fmri_surface_data_walkthrough.m)
- **Design rationale / roadmap:** [`fmri_surface_data_design_plan.md`](fmri_surface_data_design_plan.md)
- Class folder: `CanlabCore/@fmri_surface_data/`; native I/O helpers: `CanlabCore/Surface_tools/`

---

## Contents

1. [Concepts](#1-concepts)
2. [Construction](#2-construction)
3. [Properties](#3-properties)
4. [Method reference](#4-method-reference)
   - [Import / export](#import--export)
   - [Geometry & space](#geometry--space)
   - [Volume ↔ surface mapping](#volume--surface-mapping)
   - [Data operations](#data-operations)
   - [Rendering & QC](#rendering--qc)
5. [Native I/O functions](#5-native-io-functions-surface_tools)
6. [Surface spaces](#6-surface-spaces)
7. [Design notes & limitations](#7-design-notes--limitations)

---

## 1. Concepts

**Grayordinates.** The HCP "grayordinate" model represents the cortex as
*surface vertices* and the subcortex/cerebellum as *voxels*. The standard "91k"
space has 91,282 grayordinates = left + right cortex vertices (the non-medial-wall
subset of the 32,492 fs_LR-32k vertices per hemisphere) plus ~31,870 subcortical
2 mm MNI voxels. This is already a flat `[grayordinates × maps]` matrix — exactly
what `.dat` expects.

**`brain_model` replaces `volInfo`.** Cortical surface vertices have no voxel
coordinates, so `brain_model` (mirroring the CIFTI BrainModels) is the geometry
source of truth. The inherited `volInfo` slot is populated to describe **only**
the subcortical voxel sub-block (so `to_fmri_data` and volInfo-aware code work),
and is empty for surface-only objects.

**No empty-squeezing.** Unlike `fmri_data` (which drops empty voxels to save
space on sparse volumes), grayordinate data is already compact, so `.dat` is
**always** the full `[grayordinates × maps]` set. `remove_empty` / `replace_empty`
are no-ops, and masking simply zeros out-of-mask rows.

---

## 2. Construction

```matlab
obj = fmri_surface_data                       % empty object
obj = fmri_surface_data(cifti_filename)       % .dscalar/.dtseries/.dlabel.nii
obj = fmri_surface_data(gifti_filename)       % .func/.shape/.label/.surf.gii
obj = fmri_surface_data(cifti_struct)         % a canlab_read_cifti output struct
obj = fmri_surface_data(gifti_struct)         % a canlab_read_gifti output struct
obj = fmri_surface_data('dat', X, 'brain_model', bm, 'surface_space', 'fsLR_32k', ...)
obj = fmri_surface_data(image_vector_obj)     % recast a matching object
```

The constructor auto-detects CIFTI vs GIFTI by extension and reads it with the
native readers. `surface_space` and `grayordinate_type` are inferred from the
per-hemisphere vertex count.

```matlab
s = fmri_surface_data(which('transcriptomic_gradients.dscalar.nii'));
size(s.dat)            % [96854 3]
s.surface_space       % 'fsLR_32k'
```

---

## 3. Properties

| Property | Description |
|---|---|
| `dat` | `[nGrayordinates × nMaps]` single. The CIFTI cdata matrix (cortex L, cortex R, subcortical voxels). Always full (never squeezed). |
| `brain_model` | Geometry source of truth (CIFTI BrainModels). `.models{i}` has `.struct`, `.type` `'surf'`/`'vox'`, `.start`, `.count`, `.numvert`, `.vertlist` (0-based), `.voxlist` (3×N); plus `.vol` (`.dims`, `.sform`), `.grayordinate_type`, `.cluster`. |
| `geom` | Cortical mesh cache (faces/vertices), loaded for rendering. |
| `imagetype` | `'dscalar'` / `'dtseries'` / `'dlabel'` / `'func'` / `'shape'` / `'label'`. |
| `series_info` | For `.dtseries`: `.start`/`.step`/`.unit`/`.exponent`. |
| `label_table` | For `.dlabel`/`.label`: a MATLAB table (variables `key`, `name`, `rgba`). |
| `surface_space` | `'fsLR_32k'`, `'fsaverage_164k'`, … (gatekeeps `compare_space`; drives mesh/warp choice). |
| `volInfo` | Inherited; populated for the subcortical voxel sub-block only (empty for surface-only). |
| `mask` | Optional `[nGray × 1]` logical (or another same-space object) for `apply_mask`. |
| `X` / `Y` / `Y_names` / `covariates` / `covariate_names` / `images_per_session` / `metadata_table` / `image_metadata` / `additional_info` | Per-map (column) annotations, same names/roles as `fmri_data`. |
| `removed_voxels` / `removed_images` | Inherited but **vestigial** (always all-false); grayordinate data is never squeezed. |
| `image_names` / `fullpath` / `history` / `source_notes` / `dat_descrip` | Inherited provenance / file metadata. |

---

## 4. Method reference

> Type `methods(fmri_surface_data)` for the full list. Methods inherited from
> `image_vector` that operate purely on `.dat` (`get_wh_image`, `mahal`,
> `pca`, …) also apply. volInfo-dependent display/QC methods (`descriptives`,
> `montage`, `orthviews`, `flip`) are not yet surface-aware — run them on the
> volumetric part via `to_fmri_data(obj)` for now.

### Import / export

#### `write(obj, filename)`
Write to a CIFTI-2 (`.nii`) or GIFTI (`.gii`) file natively; dispatches on the
extension. Rebuilds the CIFTI maps dimension from `intent` + `image_names` /
`label_table` / `series_info`; re-emits the original XML faithfully when the
layout is unchanged, else regenerates it.

```matlab
write(s, '/tmp/out.dscalar.nii');                 % CIFTI-2
write(s, 'fname', '/tmp/out.func.gii');           % GIFTI
```

#### `to_fmri_data(obj)`
Export the subcortical/volumetric grayordinates as a standard `fmri_data` object
in MNI space (writeable to `.nii`, montageable, etc.). Cortical surface
grayordinates are dropped (use `surf2vol`).

```matlab
vol = to_fmri_data(s);     % subcortex as fmri_data
```

### Geometry & space

#### `reconstruct_image(obj)`
Returns a struct with dense per-hemisphere vertex arrays (medial wall = `NaN`)
and a `[X × Y × Z × maps]` subcortical volume.

```matlab
r = reconstruct_image(s);
r.cortex_left          % [numvert × nMaps], medial wall NaN
r.cortex_right
r.volume               % subcortical volume
r.volume_volInfo       % its SPM-style volInfo
```

#### `compare_space(obj, obj2)`
Surface analogue of `image_vector.compare_space`, preserving the integer contract:
`0` same · `1` different space (tag/layout) · `2` missing `brain_model` · `3` same
space but different in-data grayordinates.

#### `reparse_contiguous(obj, 'which_image', k)`
Label contiguous clusters of active (nonzero, non-NaN) grayordinates: cortex via
the mesh edge graph (`graph`/`conncomp`), subcortex via 26-connectivity
(`spm_clusters`). Writes integer labels to `brain_model.cluster`.

```matlab
[obj, ncl] = reparse_contiguous(threshold(s, 1.0, 'positive'));
```

### Surface resampling and volume ↔ surface mapping

#### `resample_surface(obj, target_space, 'interp', 'linear')`
Resample the **cortical** data to another surface space, natively. Target-space
keywords (also `resample_surface(obj,'list')`): `fsaverage_164k` (aliases
`fsaverage`, `164k`), `fsaverage6`/`fsaverage5`/`fsaverage4`, `fs_LR_32k` (aliases
`fsLR`, `hcp`, `32k`). fsaverage↔fs_LR uses the vendored HCP registration
("deformed") sphere so both meshes share one spherical frame and are resampled by
barycentric (or nearest) interpolation; fsaverage down-sampling is the exact nested
icosahedral subset. Interpolation weights depend only on geometry, so they are
built once and applied to all maps (a sparse matrix-multiply) — a 50-map object
costs ~the same as one map. Continuous data defaults to barycentric (`'linear'`);
binary masks and `.dlabel` images default to `'nearest'`. Subcortical voxels are
carried through unchanged.

```matlab
s32 = resample_surface(vol2surf(t), 'fsLR_32k');   % fsaverage -> HCP fs_LR-32k
s6  = resample_surface(vol2surf(t), 'fsaverage6'); % nested downsample (exact)
```

#### `vol2surf(volume_obj, 'interp', 'linear')`  *(method on `image_vector` / `fmri_data`)*
Project a volumetric image (MNI152) onto the **fsaverage-164k** cortical surface,
returning an `fmri_surface_data`. **This is the CBIG registration-fusion mapper,
natively:** it samples the vendored CBIG RF-ANTs MNI→fsaverage per-vertex
coordinates (Wu et al. 2018, *Hum Brain Mapp*) with `interpn` — a line-for-line
reimplementation of CBIG's `CBIG_RF_projectMNI2fsaverage.m` using the identical
warp, driven by SPM's `volInfo.mat` instead of FreeSurfer's `MRIread`, so no
FreeSurfer/Workbench is required. Use `'interp','nearest'` for label maps.

```matlab
ssurf = vol2surf(fmri_data('weights.nii'));        % fsaverage_164k object
ssurf = vol2surf(t, 'interp', 'nearest');          % e.g. an integer atlas
```

#### `surf2vol(obj, 'reference', ref)`
Inverse of `vol2surf`: project an **fsaverage-164k** object back to an MNI152
volume (`fmri_data`), scattering the same CBIG coordinates (`accumarray` mean).
`'reference'` (an `fmri_data`/`image_vector`) sets the target grid; default is MNI
2 mm `[91 109 91]`. The vol→surf→vol round-trip correlates ~1.0 on cortical voxels.

```matlab
backvol = surf2vol(ssurf);            % -> fmri_data in MNI
backvol = surf2vol(ssurf, 'reference', some_fmri_data);
```

### Data operations

#### `mean(obj, ['omitnan'])`
Average across maps (columns), returning a single-map object (geometry preserved).

#### `apply_mask(obj, mask, ['invert'])`
Keep a subset of grayordinates, **zeroing** the rest (no shrinking; D5b). `mask`
is a `[nGray × 1]` logical/numeric vector or another same-space `fmri_surface_data`.

```matlab
s2 = apply_mask(s, s.dat(:,1) > 0);    % zero non-positive grayordinates
```

#### `threshold(obj, thresh, ['positive'|'negative'], ['k', N])`
Zeros grayordinates inside the threshold range. `thresh` is a scalar `t` (keeps
`|value| ≥ t`) or `[lo hi]` (keeps `value ≤ lo | value ≥ hi`). With `'k', N`,
applies a **cluster-extent threshold** after the raw threshold — removing
contiguous clusters smaller than `N` grayordinates (mesh-graph connectivity for
cortex, 26-connectivity for subcortex; see `reparse_contiguous`).

```matlab
t = threshold(s, 2.0, 'positive', 'k', 20);   % >2, clusters >= 20 grayordinates
```

#### `rebuild_like(obj, newdat)`
Wrap a new `[nGray × K]` matrix into an `fmri_surface_data` carrying `obj`'s
geometry. Used internally by data-transforming methods.

### Analysis

These mirror the `fmri_data` analysis surface. `predict`, `ttest`, and `ica`
delegate to the corresponding `fmri_data`/`image_vector` methods (treating each
grayordinate as a feature), so the algorithms are identical; geometry-bearing
results are remapped back to `fmri_surface_data`.

#### `cat(obj1, obj2, ...)` / `[obj1, obj2, ...]`
Concatenate objects along the map (image) dimension — e.g. building a group
dataset from per-subject maps. All objects must share the same grayordinate
space (`compare_space == 0`). Per-map fields (`X`, `Y`, `covariates`,
`image_names`, `metadata_table`) are concatenated too.

```matlab
all_subs = cat(sub1, sub2, sub3);     % or [sub1 sub2 sub3]
```

#### `ttest(obj, [pthresh], [k])`
Grayordinate-wise one-sample t-test across maps. Returns an `fmri_surface_data`
with the t-map in `.dat` and `.p`/`.ste`/`.sig`/`.dfe` in
`.additional_info.statistic`.

```matlab
t = ttest(cat(sub1, sub2, sub3));
surface(t, 'clim', [-5 5]);
```

#### `regress(obj, [X])`
Grayordinate-wise OLS regression of the maps onto design `X` (`[nMaps × p]`;
defaults to `obj.X`; no intercept added automatically). Returns coefficients in
`.dat` (`[nGray × p]`) and `.t`/`.p`/`.se`/`.dfe` in `.additional_info.statistic`.
For contrasts/diagnostics, use `to_fmri_data` + `fmri_data.regress`/`glm_map`.

```matlab
obj.X = [ones(n,1) age(:)];
b = regress(obj);
surface(get_wh_image(b, 2));          % the age effect
```

#### `predict(obj, 'algorithm_name', ..., 'nfolds', ...)`
Cross-validated multivariate prediction (set `obj.Y` first). Returns
`[cverr, stats, optout]` as `fmri_data.predict`; `stats.weight_obj` is remapped
to an `fmri_surface_data` you can `surface`/`write`.

```matlab
obj.Y = scores(:);
[err, stats] = predict(obj, 'algorithm_name', 'cv_lassopcr', 'nfolds', 5);
surface(stats.weight_obj);
```

#### `ica(obj, [nIC])`
Spatial ICA; returns an `fmri_surface_data` whose maps are the spatial
components. **Requires the GIFT/`icatb` toolbox (`icatb_fastICA`)** — the one
method that is not fully self-contained.

### Parcellation & regions

#### `apply_parcellation(obj, parcels, ['area'])`
Average each map within the parcels of a surface atlas. `parcels` is another
`fmri_surface_data` (a `.dlabel` on the same space; its `label_table` supplies
names) or an integer key vector `[nGray × 1]`. Key 0 / NaN = background (medial
wall) and is excluded. `'area'` area-weights by per-vertex surface area.

```matlab
atl = fmri_surface_data(which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii'));
[pm, labels, tbl] = apply_parcellation(s, atl);   % pm: [nMaps × nParcels]
```

Returns `parcel_means` `[nMaps × nParcels]`, a `labels` cell, and a summary
`table` (`key`, `label`, `n_grayordinates`, `total_weight`).

#### `surface_region(obj, ['which_image', k])`
Summarize contiguous clusters (active = nonzero, non-NaN) as a struct array — the
surface analogue of `region()`. Each element has `.struct`, `.type`,
`.cluster_id`, `.grayord_rows`, `.vertex_indices`, `.XYZmm` (centroid),
`.numVox`, `.val`. (For a full CANlab `region` object on the subcortical part,
use `region(to_fmri_data(obj))`.)

```matlab
reg = surface_region(threshold(t, 3, 'positive', 'k', 20));
[reg.numVox]            % cluster sizes
```

### Rendering & QC

#### `surface(obj, ...)`
Render on cortical surfaces. Three modes:

- **Native (default):** builds a **managed, stateful `fmridisplay`** whose surface
  views are a four-surface set (L/R lateral + medial) **chosen to match the object's
  `surface_space`** (table below), and paints the data as a managed surface-native
  layer. Because the returned object is under a controller, `set_colormap` /
  `set_opacity` / `rethreshold` / `removeblobs` / `refresh` act on the surfaces (this
  is how you change the colormap after rendering). Each hemisphere shows its own data.
- **`'existingsurface', handles`:** color patch handles you already have (returns a
  struct with `.figure`/`.axes`/`.surfaces`).
- **`'mni_surface', name`:** create an `addbrain` surface (e.g. `'left'`,
  `'hcp inflated'`) and render onto it, projecting through a volume when the
  surface is not the object's native mesh (returns a struct).

**Default four-surface view per surface space.** With no explicit surface argument,
`surface(obj)` picks the view that matches the object's space, so the data renders at
native (or near-native) fidelity automatically:

| `surface_space` | verts/hemi | default four-surface view | how it paints |
|---|---|---|---|
| `fsLR_32k` | 32,492 | `foursurfaces_hcp` | **native** (direct, no resampling) |
| `fsaverage_164k` | 163,842 | `foursurfaces_freesurfer` | **native** (direct, no resampling) |
| `fsaverage6` / `fsaverage5` / `fsaverage4` | 40,962 / 10,242 / 2,562 | `foursurfaces_freesurfer` | resampled **up** to fsaverage-164k (nearest, cached) |
| `onavg_41k` / `onavg_10k` | 40,962 / 10,242 | `foursurfaces_hcp` | resampled to fs_LR-32k (nearest, cached) |

Only `fsLR_32k` and `fsaverage_164k` ship inflated display meshes (in
`canlab_canonical_brains/Canonical_brains_surfaces/`). The other spaces need no new
mesh assets: the nested fsaverage spaces are subsets of fsaverage-164k, and onavg is
fs_LR-aligned via the bundled onavg spheres, so each renders on its parent/aligned
display space through a fast, cached nearest-neighbour resample (a one-line
`render_layer_surfaces: resampling …` message prints the first time).

You can also add a surface object to an **existing** managed display:
`o2 = surface(o2, obj)` adds the space-matched native surfaces and paints it;
`o2 = surface(o2, obj, 'foursurfaces_hcp')` (or any `addbrain` surface keyword)
overrides that with a named surface set. Requesting a **different standard** mesh than
the data (e.g. fsaverage meshes for fs_LR data) auto-resamples the object onto that
mesh (the same fast, cached nearest-neighbour resample, with the printed message); an
**arbitrary** (non-standard) MNI mesh is handled by projecting the object to a volume.
The `fmridisplay:render_layer_surfaces:spacemismatch` warning fires only when a mesh is
neither a recognized standard cortical mesh nor projectable.

| Option | Meaning |
|---|---|
| `'surftype'` | `'inflated'` (default), `'midthickness'`, `'sphere'` (fs_LR). |
| `'which_image'` | map (column) to render (default 1). |
| `'clim'` / `'cmaprange'` | `[lo hi]` color limits (default symmetric from data). |
| `'colormap'` | named map (`'hot'`, `'parula'`, …) or `[n × 3]` matrix. |
| `'pos_colormap'` / `'neg_colormap'` | `[n × 3]` colormaps (default hot / cool). |

```matlab
o2 = surface(s, 'which_image', 1);                 % managed fs_LR, 4 views
o2 = set_colormap(o2, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);  % recolor live
o2 = removeblobs(o2);                              % back to anatomy

o2 = montage(fmridisplay, 'axial');                % existing managed display
o2 = surface(o2, s);                               % add + paint its native surfaces

surface(ssurf, 'mni_surface', 'left');             % on an addbrain MNI surface
han = addbrain('hcp inflated left');               % an fs_LR mesh
surface(s, 'existingsurface', han);                % color it directly
```

#### `render_on_surface(obj, handles, ...)`
Lower-level: color given patch `handles`. Patches whose vertex count matches a
hemisphere are colored **directly** (true native-space rendering — works for any
matching-space mesh, e.g. all `addbrain` fs_LR/fsaverage surfaces); other surfaces
are colored by projecting the object to a volume and reusing
`image_vector.render_on_surface`. Options: `'which_image'`, `'clim'`,
`'pos_colormap'`, `'neg_colormap'`.

#### `plot(obj, ['norender'])`
Quick QC panel: value histogram, per-map mean ± sd, coverage, and a mean-map
surface render.

---

## 5. Native I/O functions (`Surface_tools/`)

These standalone functions back the object and can be used directly. They depend
only on core MATLAB (+ the JVM for gzip) — no external toolbox.

| Function | Purpose |
|---|---|
| `canlab_read_cifti(file)` | Read CIFTI-2 (`.dscalar`/`.dtseries`/`.dlabel`/`.dconn.nii`) → struct (`.cdata`, `.diminfo`, `.intent`, `.hdr`, `.xml`). |
| `canlab_write_cifti(file, cii)` | Write CIFTI-2 (faithful re-emit or regenerate). |
| `canlab_read_gifti(file)` | Read GIFTI (`.surf`/`.func`/`.shape`/`.label.gii`) → struct (`.vertices`, `.faces`, `.cdata`, `.labels`). |
| `canlab_write_gifti(file, gii)` | Write GIFTI (ASCII / Base64 / GZipBase64). |
| `canlab_surface_vertexcolors(vals, clim, poscm, negcm, graycolor)` | Map per-vertex values → truecolor (split hot/cool; NaN/zero = `graycolor`, default gray). |
| `canlab_cbig_warp_path(name)` | Resolve the vendored CBIG RF warp paths. |

---

## 6. Surface spaces

| Space | Per-hemi vertices | Source | Bundled display mesh |
|---|---|---|---|
| `fsLR_32k` | 32,492 | native CIFTI (HCP grayordinates) | inflated (`foursurfaces_hcp`) |
| `fsaverage_164k` | 163,842 | `vol2surf` (CBIG RF mapping) | inflated (`foursurfaces_freesurfer`) |
| `fsaverage6` | 40,962 | `resample_surface` (nested subset of 164k) | — (resamples to 164k) |
| `fsaverage5` | 10,242 | `resample_surface` (nested subset of 164k) | — (resamples to 164k) |
| `fsaverage4` | 2,562 | `resample_surface` (nested subset of 164k) | — (resamples to 164k) |
| `onavg_41k` | 40,962 | `resample_surface` (equal-area, fs_LR-aligned) | — (resamples to fs_LR) |
| `onavg_10k` | 10,242 | `resample_surface` (equal-area, fs_LR-aligned) | — (resamples to fs_LR) |

Native CIFTI is `fsLR_32k`; `vol2surf` produces `fsaverage_164k`. fs_LR and fsaverage
have different mesh topologies, but `resample_surface` maps natively between **all** of
the spaces above (nested fsaverage downsampling is exact; every other direction uses a
spherical barycentric/nearest interpolation on a shared fs_LR-aligned frame — see
`resample_surface(obj, 'list')`). Only `fsLR_32k` and `fsaverage_164k` ship inflated
display meshes in `canlab_canonical_brains/Canonical_brains_surfaces/`; the other
spaces render by resampling onto their parent/aligned display space (see the default
four-surface view table under [`surface`](#4-method-reference)). Native rendering uses
the matching bundled mesh directly (`addbrain('hcp inflated')` for fs_LR,
`addbrain('inflated')` for fsaverage).

---

## 7. Design notes & limitations

- **Group-template mapping = CBIG registration fusion.** `vol2surf`/`surf2vol`
  are native reimplementations of the CBIG RF-ANTs MNI152↔fsaverage mappers
  (Wu et al. 2018), using the vendored warp under `Cifti_plotting/
  CBIG_registration_fusion_surf2vol_vol2surf/` — `vol2surf` ≡
  `CBIG_RF_projectMNI2fsaverage`. So you are already using CBIG registration
  fusion, with no FreeSurfer/CBIG-toolbox dependency. It is a fixed group
  correspondence (correct for group MNI maps; not a per-subject ribbon mapper).
  CBIG's heavier `fsaverage2Vol` ribbon-fill needs FreeSurfer + the MARS toolbox +
  external mask geometry (not bundled) and is intentionally not wired in.
- **fs_LR ↔ fsaverage.** Not yet bridged; keep a workflow in one space, or go
  through a volume.
- **Planned:** a colorbar/outline overlay on the native render, and surface-aware
  overrides of `montage`/`orthviews`/`flip`/`descriptives`; the volumetric escape
  hatch (`to_fmri_data`) covers those cases today. (Cluster-extent thresholding is
  already implemented — `threshold(obj, thr, 'k', N)`.)
- See [`fmri_surface_data_design_plan.md`](fmri_surface_data_design_plan.md) for the
  full milestone roadmap and the rationale behind these choices.
