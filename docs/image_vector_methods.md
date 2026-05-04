# `image_vector` methods, organized by area

This is a functional index of methods available on an `image_vector` object.
`image_vector` is the superclass for the CANlab neuroimaging-data classes
(`fmri_data`, `statistic_image`, `atlas`, `fmri_mask_image`); you will
rarely create an `image_vector` directly. The methods listed here are
inherited by all of those subclasses, and most appear with `From:
@image_vector` in `fmri_data_methods.md`.

`image_vector` stores brain image data in a flat (2-D) voxels x images
matrix together with the meta-data (`volInfo`) needed to round-trip back
to 3-D image space. It supports image math, masking, resampling,
visualization, statistics, multivariate prediction, region/atlas
extraction, and I/O. Type `methods(my_obj)` in MATLAB for the live list
on any instance.

## Properties

| Property | Description |
|---|---|
| `source_notes` | Free-text notes about the data source |
| `dat` | Image data, a [voxels x images] single-precision matrix |
| `dat_descrip` | String description of the dataset |
| `volInfo` | Brain mask + voxel-to-world mapping (mat, xyzlist, cluster) |
| `removed_voxels` | Logical vector of empty in-mask voxels removed (saves space; see `remove_empty`/`replace_empty`) |
| `removed_images` | Logical vector of images removed (saves space; see `remove_empty`/`replace_empty`) |
| `image_names` | List of image names loaded into the object, no paths |
| `fullpath` | List of image names with full paths; used by `write` |
| `files_exist` | Logical vector: do files in `fullpath` exist on disk? |
| `history` | Cell array history of object processing, for provenance |

## Basic image math and operations

Combining images, set-like operations, voxel-level arithmetic, and the
"removed/replaced empty" state machine that several methods rely on.

| Method | From | One-liner |
|---|---|---|
| `get_wh_image` | `@image_vector` | Pick a subset of images by index |
| `mean` | `@image_vector` | Voxel-wise mean across images |
| `prctile` | `@image_vector` | Voxel-wise percentile thresholds across images |
| `image_math` | `@image_vector` | General Boolean / arithmetic between objects |
| `union` | `@image_vector` | Union/intersection masks for two objects |
| `flip` | `@image_vector` | L-R flip of all images in the object |
| `trim_mask` | `@image_vector` | Drop empty/zero rows from `.dat` |
| `select_voxels_by_value` | `@image_vector` | Keep voxels whose value matches a predicate |
| `winnerTakeAll` | `@image_vector` | Per-voxel argmax across images |
| `replace_empty` | `@image_vector` | Re-expand `.dat` to full padded voxel space |
| `remove_empty` | `@image_vector` | Compress `.dat` by dropping empty voxels |
| `reparse_contiguous` | `@image_vector` | Recompute contiguous-voxel groupings |
| `rebuild_volinfo_from_dat` | `@image_vector` | Reconstruct `volInfo` from `.dat` (recovery) |
| `reconstruct_image` | `@image_vector` | 2-D `.dat` to 3-D / 4-D MATLAB array |
| `expand_into_atlas_subregions` | `@image_vector` | Replicate into atlas-defined subregions |
| `subdivide_by_atlas` | `@image_vector` | Split into sub-objects, one per atlas region |

## Display and visualization

Most of these create a figure or modify an existing `fmridisplay`. Many
require a graphics environment.

| Method | From | One-liner |
|---|---|---|
| `histogram` | `@image_vector` | Histogram of `.dat` values, per image |
| `montage` | `@image_vector` | Slice montage on canonical anatomy |
| `slices` | `@image_vector` | One-slice-per-image montage |
| `display_slices` | `@image_vector` | 3-pane (ax/cor/sag) compact slice view |
| `slice_movie` | `@image_vector` | Movie of slices through the volume |
| `rmssd_movie` | `@image_vector` | Movie of frame-to-frame RMSSD |
| `orthviews` | `@image_vector` | SPM-style orthviews (requires SPM graphics) |
| `surface` | `@image_vector` | Render on cortical surface |
| `render_on_surface` | `@image_vector` | Lower-level surface render with options |
| `render_on_cerebellar_flatmap` | `@image_vector` | Cerebellar SUIT flatmap rendering |
| `isosurface` | `@image_vector` | 3-D isosurface from voxel volume |
| `pattern_surf_plot_mip` | `@image_vector` | Axial maximum-intensity-projection surface |
| `wedge_plot_by_atlas` | `@image_vector` | Polar/wedge plot by atlas region |
| `plot_current_orthviews_coord` | `@image_vector` | Print MNI of current orthviews crosshair |

## Resampling and interpolation

Anything that changes the voxel grid or compares grids.

| Method | From | One-liner |
|---|---|---|
| `resample_space` | `@image_vector` | Resample to another object's space |
| `resample_space_simple_reference` | `@image_vector` | Simpler resampling against a reference |
| `resample_time` | `@image_vector` | Resample / interpolate along the image (time) axis |
| `interpolate` | `@image_vector` | Fill missing values via 3-D linear interpolation |
| `compare_space` | `@image_vector` | Diagnostic: are two objects in the same space? |
| `define_space_mapping` | `@image_vector` | Build the mapping between two spaces |

## Statistics

Voxel-wise inference, prediction, and multivariate analyses.

| Method | From | One-liner |
|---|---|---|
| `searchlight` | `@image_vector` | Spherical-searchlight prediction/classification |
| `searchlightLukas` | `@image_vector` | Variant of searchlight |
| `ica` | `@image_vector` | Spatial ICA |
| `pca` | `@image_vector` | Spatial PCA |
| `mahal` | `@image_vector` | Mahalanobis distance per image vs. set |

## Tables

Producing tabular reports from images and stat maps.

| Method | From | One-liner |
|---|---|---|
| `table` | `@image_vector` | Atlas-labeled table of regions in a stat map |
| `table_of_atlas_regions_covered` | `@image_vector` | Coverage table against an atlas |
| `print_publication_table` | `@image_vector` | Pre-formatted publication-style results table |

## Annotation with spatial similarity

Comparing images or regions to reference atlases / signatures /
meta-analytic maps for interpretation.

| Method | From | One-liner |
|---|---|---|
| `image_similarity_plot` | `@image_vector` | Cosine/correlation similarity vs. a basis set |
| `image_similarity_plot_bucknermaps` | `@image_vector` | Convenience wrapper for Buckner-network maps |
| `hansen_neurotransmitter_maps` | `@image_vector` | Hansen neurotransmitter map similarity |

## Data extraction

Pulling values out of images, by mask / atlas / parcellation / coordinate.

| Method | From | One-liner |
|---|---|---|
| `apply_mask` | `@image_vector` | Restrict object to voxels in a mask |
| `apply_parcellation` | `@image_vector` | Mean-or-pattern expression per parcel |
| `extract_roi_averages` | `@image_vector` | Average per contiguous region (or mask values) |
| `extract_gray_white_csf` | `@image_vector` | Mean + top-5 components in GM, WM, CSF |
| `get_xyzmm_coordinates` | `@image_vector` | In-mask voxel indices to MNI mm coordinates |

## Data processing

Workflows for transforming and processing data objects.

| Method | From | One-liner |
|---|---|---|
| `preprocess` | `@image_vector` | Many preprocessing options (filter / outliers / scale) |

## Quality control

Diagnosing and cleaning a dataset before analysis.

| Method | From | One-liner |
|---|---|---|
| `descriptives` | `@image_vector` | Print summary stats for the dataset |
| `qc_metrics_second_level` | `@image_vector` | QC metrics across a 2nd-level set |
| `outliers` | `@image_vector` | Detect outlier images |
| `outliers_xval` | `@image_vector` | Outlier detection with cross-validated threshold |
| `jackknife_similarity` | `@image_vector` | Leave-one-out spatial similarity |

## Misc utilities

I/O, type management, provenance, threshold helpers.

| Method | From | One-liner |
|---|---|---|
| `write` | `@image_vector` | Write `.dat` back to NIfTI / Analyze on disk |
| `read_from_file` | `@image_vector` | Re-read pixel data from `.fullpath` into `.dat` |
| `check_image_filenames` | `@image_vector` | Validate `.fullpath` entries exist on disk |
| `enforce_variable_types` | `@image_vector` | Cast `.dat` and friends to canonical types |
| `history` | `@image_vector` | Show the `.history` provenance log |
| `isempty` | `@image_vector` | True iff `.dat` is empty |
| `threshold` | `@image_vector` | Threshold values (typically used on statistic_image) |
| `unstack_by_condition` | `@image_vector` | Split into sub-objects by a condition vector |
