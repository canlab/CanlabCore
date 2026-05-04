# `statistic_image` methods, organized by area

This is a functional index of methods available on a `statistic_image`
object. `statistic_image` is a subclass of `image_vector` for storing
and manipulating statistical images (t-maps, p-maps, thresholded
results). Multiple images can be stored in a single object, thresholding
is reversible (via the `threshold` method) without changing the
underlying values, and standard visualizations (`orthviews`, `surface`,
`montage`) honor the `.sig` field so only suprathreshold voxels are
shown.

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
| `orthviews` | `@statistic_image` | Orthviews honoring `.sig` so only suprathreshold voxels show (override) |
| `riverplot` | `@statistic_image` | Riverplot of relationships among images in the object |

## Statistics

| Method | From | One-liner |
|---|---|---|
| `threshold` | `@statistic_image` | Threshold based on statistical p-values, FDR, extent, etc. (reversible; override) |
| `multi_threshold` | `@statistic_image` | Multiple-threshold visualization for nested significance |
| `conjunction` | `@statistic_image` | Conjunction of two thresholded statistic_images (positive and negative separately) |
| `estimateBayesFactor` | `@statistic_image` | Voxel-wise Bayes Factors for t-tests, correlations, or proportions |

## Tables

| Method | From | One-liner |
|---|---|---|
| `table` | `@statistic_image` | Print a table of labeled regions from a thresholded statistic_image (override) |

## Misc utilities

| Method | From | One-liner |
|---|---|---|
| `check_properties` | `@statistic_image` | Fill in empty `.p`, `.ste`, `.sig` fields to the full in-mask grid (under construction) |
