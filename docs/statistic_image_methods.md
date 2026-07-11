# `statistic_image` methods, organized by area

This is a functional index of methods available on a `statistic_image`
object. `statistic_image` is a subclass of `image_vector` for storing
and manipulating statistical images (t-maps, p-maps, thresholded
results). Multiple images can be stored in a single object, thresholding
is reversible (via the `threshold` method) without changing the
underlying values, and standard visualizations (`orthviews`, `surface`,
`montage`) honor the `.sig` field so only suprathreshold voxels are
shown. `orthviews_niivue` pops the thresholded map open as an interactive
web viewer in your browser.

This class also inherits all `image_vector` methods listed in
[image_vector_methods.md](image_vector_methods.md). Only methods owned
by `statistic_image` (defined in `@statistic_image/`) are listed below;
where a method also exists in `@image_vector/`, the `statistic_image`
override wins. Type `methods(my_obj)` in MATLAB for the live list on
any instance.

## Properties

`statistic_image` inherits all `image_vector` properties (`dat`,
`volInfo`, `removed_voxels`, `removed_images`, `image_names`, `fullpath`,
`files_exist`, `history`, `dat_descrip`, `source_notes`) and adds the
following:

| Property | Description |
|---|---|
| `type` | String with image type: `generic`, `t`, `p`, `r`, `robreg` |
| `p` | Matrix of p-values for images, [voxels x images] |
| `p_type` | String with source info for p-values |
| `ste` | Matrix of standard error values for images, [voxels x images] |
| `threshold` | Latest statistical threshold applied |
| `thr_type` | Information about the statistical threshold applied |
| `sig` | Logical matrix of which voxels are significant, [voxels x images] |
| `N` | Sample size |
| `dfe` | Error degrees of freedom for the test |
| `image_labels` | Semantic labels for each image, cell array |

## Basic image math and operations

| Method | From | One-liner |
|---|---|---|
| `cat` | `@statistic_image` | Concatenate statistic_image objects, resampling to first as needed |
| `select_one_image` | `@statistic_image` | Select a single image, reducing `.p`, `.ste`, `.sig`, `.dat` accordingly |
| `reparse_contiguous` | `@statistic_image` | Re-build the contiguous-cluster index for the statistic_image (override) |
| `convert2mask` | `@statistic_image` | Convert to fmri_mask_image based on the `.sig` field |

## Display and visualization

| Method | From | One-liner |
|---|---|---|
| [`orthviews`](individual_functions/fmri_data_orthviews.md) | `@statistic_image` | Orthviews honoring `.sig` so only suprathreshold voxels show (override) |
| `orthviews_niivue` | `@statistic_image` | Web "orthviews": write a self-contained NiiVue `.html` of the (thresholded) map and open it in the browser — a one-line wrapper over `canlab_niivue` |
| [`riverplot`](individual_functions/statistic_image_riverplot.md) | `@statistic_image` | Riverplot of relationships among images in the object |

**Interactive viewers.** `orthviews_niivue(t)` is the quickest way to inspect a thresholded map in a browser — it writes a portable NiiVue page (crosshair coordinate + value + **atlas region** readout, with a single-region outline/shade) to a temp folder and opens it. The underlying [`canlab_niivue(t)`](canlab_niivue_guide.md) exposes the full option set (colormaps, color limits, output folder, embedding). In MATLAB, `canlab_orthviews(t)` opens the enhanced SPM-style three-plane viewer (multiple blob layers, region tables, and the same crosshair atlas region-name readout).

## Statistics

| Method | From | One-liner |
|---|---|---|
| [`threshold`](individual_functions/statistic_image_threshold.md) | `@statistic_image` | Threshold based on statistical p-values, FDR, extent, etc. (reversible; override) |
| [`multi_threshold`](individual_functions/statistic_image_multi_threshold.md) | `@statistic_image` | Multiple-threshold visualization for nested significance |
| `conjunction` | `@statistic_image` | Conjunction of two thresholded statistic_images (positive and negative separately) |
| `estimateBayesFactor` | `@statistic_image` | Voxel-wise Bayes Factors for t-tests, correlations, or proportions |

## Tables

| Method | From | One-liner |
|---|---|---|
| [`table`](individual_functions/statistic_image_table.md) | `@statistic_image` | Print a table of labeled regions from a thresholded statistic_image (override) |

## Misc utilities

| Method | From | One-liner |
|---|---|---|
| `check_properties` | `@statistic_image` | Fill in empty `.p`, `.ste`, `.sig` fields to the full in-mask grid (under construction) |
