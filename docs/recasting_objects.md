# Recasting (converting) between object types

CanlabCore has several closely related object classes. Many analyses naturally produce one type and consume another, so most conversions are short, well-defined idioms. This page lists the canonical recipes.

The relevant classes are:

```
image_vector  (abstract base; rarely used directly)
├── fmri_data           generic image data + .X, .Y, covariates
├── statistic_image     stat maps with t / p / sig
├── atlas               labeled parcellation with probability maps
└── fmri_mask_image     binary mask (legacy)

region                  list of contiguous clusters (not a subclass of image_vector)
fmridisplay             figure-handle container (not image data)
```

See [Object_methods.md](Object_methods.md) for the full per-class method listings.

## Quick reference

| From | To | Idiom |
|---|---|---|
| `fmri_data` | `statistic_image` | `t = ttest(obj)`, `t = regress(obj, ...)`, `t = signtest(obj)` |
| `fmri_data` | `region` | `r = region(ttest(obj))` (via stat map) |
| `fmri_data` | `atlas` | use `apply_atlas` to get values per region; not a direct conversion |
| `fmri_data` | `image_vector` | `iv = image_vector('dat', obj.dat, 'volInfo', obj.volInfo, ...)` |
| `image_vector` | `fmri_data` | `obj = fmri_data(iv)` |
| `statistic_image` | `region` | `r = region(t)` (uses .sig to keep only suprathreshold voxels) |
| `statistic_image` | `fmri_data` | `obj = fmri_data(t)` (drops .sig / .p; keeps `.dat`) |
| `statistic_image` | `fmri_mask_image` / mask | `m = convert2mask(t)` |
| `atlas` | `region` | `r = atlas2region(atl)` (one element per atlas region) |
| `atlas` | `fmri_data` | `obj = fmri_data(atl)`; or `parcel_data2fmri_data(atl, values)` to broadcast per-parcel values to voxels |
| `atlas` | `statistic_image` | `parcel_stats2statistic_image(atl, stat_per_parcel, ...)` |
| `region` | `fmri_data` | `obj = region2fmri_data(r)` |
| `region` | `image_vector` | `iv = region2imagevec(r)` |
| `region` | `atlas` | `atl = region2atlas(r)` |
| `region` | `struct` | `s = region2struct(r)` (plain MATLAB struct, e.g. for saving) |
| `brainpathway` | `fmri_data` | `obj = brainpathway2fmri_data(bp)` |
| any | display layer | `o2 = canlab_results_fmridisplay(obj)` (creates `fmridisplay`) |

## Common patterns

### Group t-test workflow: `fmri_data` to `statistic_image` to `region`

The most common chain in the toolbox:

```matlab
imgs = load_image_set('emotionreg');     % fmri_data, 30 contrasts
t    = ttest(imgs);                      % statistic_image
t    = threshold(t, 0.005, 'unc');       % statistic_image with .sig set
r    = region(t);                        % region: one entry per blob
table(r);                                % atlas-labeled results table
montage(r, 'regioncenters');             % per-blob mini-montage
```

`ttest` does the recast from `fmri_data` to `statistic_image`. `region(t)` does the recast from `statistic_image` to `region` and uses the threshold state to define cluster boundaries.

### From an atlas to per-region data

If you have voxel-level data and want to analyze it parcel-wise, you typically don't recast the data object. Instead you keep the `fmri_data` and apply the atlas to it:

```matlab
atl     = load_atlas('canlab2024');                  % atlas
imgs    = load_image_set('emotionreg');              % fmri_data
parcel_means = apply_atlas(imgs, atl);               % images x parcels matrix
% or
parcel_obj   = apply_parcellation(imgs, atl);        % parcel-level fmri-data-like object
```

If you want the atlas itself as `fmri_data`-like data, use the constructor:

```matlab
obj = fmri_data(atl);                                 % one image, integer labels in .dat
```

If you have **per-parcel statistics** (e.g., a vector of one t-value per region) and want to project them back to voxel space:

```matlab
voxel_obj  = parcel_data2fmri_data(atl, parcel_values);    % voxelwise broadcast
voxel_stat = parcel_stats2statistic_image(atl, t_per_parcel, p_per_parcel);
```

### From regions back to whole-brain images

Sometimes you build a `region` object (for example, by clustering or by manual ROI definition) and then need to feed it into a method that expects an `fmri_data` or an `image_vector`. The two recipes:

```matlab
obj = region2fmri_data(r);     % full fmri_data with .dat at the region voxels
iv  = region2imagevec(r);      % bare image_vector (lighter; no .X / .Y)
```

`region2atlas(r)` is useful if your `region` object came from clustering and you want to assign each cluster an integer label and call it an `atlas` going forward.

### From a statistic_image to a binary mask

Many CanlabCore functions accept either a binary mask or a statistic_image as a "mask" argument, but if you specifically want a clean binary mask (e.g., to write to disk or apply to a different dataset), use:

```matlab
m = convert2mask(t);           % binary fmri_mask_image / fmri_data
```

This is equivalent to thresholding and then keeping only the suprathreshold voxels.

### Display layer: any image to fmridisplay

`fmridisplay` is not an image-data class — it is a container of figure handles for slice montages and surfaces. The way you usually create one is:

```matlab
o2 = canlab_results_fmridisplay(t);            % creates the figure with t-map blobs
o2 = removeblobs(o2);                          % strip the blobs, keep the figure
o2 = addblobs(o2, region(t2));                 % overlay a different result
```

So strictly speaking you don't recast an image object to `fmridisplay`; you build an `fmridisplay` from your image data via `canlab_results_fmridisplay` (or by hand) and then keep using it as a layered display.

## Notes on `.dat` and `volInfo` consistency

When you recast image classes by passing one to another's constructor (e.g., `fmri_data(image_vector_obj)`), the `.dat`, `.volInfo`, `removed_voxels`, and `removed_images` fields are copied directly. If the source object had been compressed (some voxels in `removed_voxels`), the destination inherits that state. If you want a clean, padded copy, call `replace_empty(obj)` before recasting:

```matlab
obj_padded = replace_empty(obj);
new_obj    = fmri_data(obj_padded);
```

See the discussion of `removed_voxels` / `removed_images` and `replace_empty` / `remove_empty` in [`fmri_data_methods.md`](fmri_data_methods.md) and `CanlabCore/CLAUDE.md` — most recast bugs come from forgetting this step.
