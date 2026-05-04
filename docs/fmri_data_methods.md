# `fmri_data` methods, organized by area

This is a functional index of methods available on an `fmri_data` object.
Many methods are inherited from `image_vector` (the superclass) and apply to
its other subclasses (`statistic_image`, `atlas`, `fmri_mask_image`) too —
the inheritance source is shown in the column. The on-disk organization is
per-class (`@fmri_data/`, `@image_vector/`); this document is the
cross-cutting functional view.

Type `methods(my_obj)` in MATLAB for the live list on any instance.

## Basic image math and operations

Combining images, set-like operations, voxel-level arithmetic, and the
"removed/replaced empty" state machine that several methods rely on.

| Method | From | One-liner |
|---|---|---|
| `cat` | `@fmri_data` | Concatenate fmri_data objects along the image (4th) dim |
| `horzcat` | both | `[a, b]` operator on objects, image-wise |
| `split` | `@fmri_data` | Split object into N subsets by image index |
| `get_wh_image` | `@image_vector` | Pick a subset of images by index |
| `mean` | `@image_vector` | Voxel-wise mean across images |
| `prctile` | `@image_vector` | Voxel-wise percentile across images |
| `plus`, `minus`, `power` | `@image_vector` | Element-wise arithmetic on objects |
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
| `reconstruct_image` | `@image_vector` | 2-D `.dat` → 3-D / 4-D MATLAB array |
| `expand_into_atlas_subregions` | `@image_vector` | Replicate into atlas-defined subregions |
| `subdivide_by_atlas` | `@image_vector` | Split into sub-objects, one per atlas region |

## Display and visualization

Most of these create a figure or modify an existing `fmridisplay`. Many
require a graphics environment (won't fully work in headless `matlab -batch`).

| Method | From | One-liner |
|---|---|---|
| `plot` | `@fmri_data` | Custom multi-panel QC plot for an fmri_data object |
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
| `riverplot` | `@fmri_data` | Riverplot visualization |
| `wedge_plot_by_atlas` | `@image_vector` | Polar/wedge plot by atlas region |
| `plot_current_orthviews_coord` | `@image_vector` | Print MNI of current orthviews crosshair |
| `saveplots` | `@fmri_data` | Save the current plots to disk |

## Resampling and interpolation

Anything that changes the voxel grid or compares grids.

| Method | From | One-liner |
|---|---|---|
| `resample_space` | `@image_vector` | Resample to another object's space |
| `resample_space_simple_reference` | `@image_vector` | Simpler resampling against a reference |
| `resample_time` | `@image_vector` | Resample / interpolate along the image (time) axis |
| `interpolate` | `@image_vector` | General interpolation in voxel space |
| `compare_space` | `@image_vector` | Diagnostic: are two objects in the same space? |
| `define_space_mapping` | `@image_vector` | Build the mapping between two spaces |
| `spm_coregister` | `@fmri_data` | Wrapper around SPM's coregistration |

## Statistics

Voxel-wise inference, prediction, and multivariate analyses.

| Method | From | One-liner |
|---|---|---|
| `ttest` | `@fmri_data` | Voxelwise one-sample t-test → statistic_image |
| `signtest` | `@fmri_data` | Voxelwise non-parametric sign test |
| `regress` | `@fmri_data` | Voxelwise multiple regression (uses `obj.X`) |
| `robfit_parcelwise` | `@fmri_data` | Robust regression at parcel level |
| `searchlight` | `@image_vector` | Spherical-searchlight prediction/classification |
| `searchlightLukas` | `@image_vector` | Variant of searchlight |
| `ica` | `@image_vector` | Spatial ICA |
| `pca` | `@image_vector` | Spatial PCA |
| `mahal` | `@image_vector` | Mahalanobis distance per image vs. set |
| `fitlme_voxelwise` | `@fmri_data` | Voxelwise mixed-effects model |
| `rsa_regression` | `@fmri_data` | Representational-similarity regression |
| `hrf_fit` | `@fmri_data` | HRF estimation per voxel |
| `test_generalizability` | `@fmri_data` | Cross-cohort generalization test |
| `canlab_connectivity_preproc` | `@fmri_data` | Connectivity-prep pipeline (despike / nuisance / bandpass) |

## Multivariate prediction

| Method | From | One-liner |
|---|---|---|
| `predict` | `@fmri_data` | Cross-validated multivariate prediction |
| `predict_test_suite` | `@fmri_data` | Battery of `predict` sanity checks |
| `evaluate_spatial_scale` | `@fmri_data` | Spatial-scale evaluation for patterns |
| `bootstrap_structure_coeff_diff` | `@fmri_data` | Bootstrap diff in structure coefficients |
| `structure_coefficient_map` | `@fmri_data` | Per-voxel structure coefficient map |
| `structure_coefficients` | `@fmri_data` | Structure coefficients |
| `get_model_encoding_map` | `@fmri_data` | Encoding map for a fitted model |
| `model_brain_pathway` | `@fmri_data` | Pathway-modeling on connectivity |
| `model_mpathi` | `@fmri_data` | Multi-pathway model |
| `dual_regression` | `@fmri_data` | Dual regression on group ICA components |

## Tables

Producing tabular reports from images and stat maps.

| Method | From | One-liner |
|---|---|---|
| `table` | `@image_vector` | Atlas-labeled table of regions in a stat map |
| `table_of_atlas_regions_covered` | `@image_vector` | Coverage table against an atlas |
| `print_publication_table` | `@image_vector` | Pre-formatted publication-style results table |
| `ttest_table_and_lateralization_test` | `@fmri_data` | Combined t-table + L/R lateralization test |

## Annotation with spatial similarity

Comparing images or regions to reference atlases / signatures /
meta-analytic maps for interpretation.

| Method | From | One-liner |
|---|---|---|
| `image_similarity_plot` | `@image_vector` | Cosine/correlation similarity vs. a basis set |
| `image_similarity_plot_bucknermaps` | `@image_vector` | Convenience wrapper for Buckner-network maps |
| `annotate_binary_results_map` | `@fmri_data` | Annotate a thresholded map with atlas/signature labels |
| `annotate_continuous_neuroimage_maps` | `@fmri_data` | Annotation for continuous (unthresholded) maps |
| `neurosynth_feature_labels` | `@fmri_data` | Look up Neurosynth feature labels |
| `neurosynth_lexical_plot` | `@fmri_data` | Lexical plot of Neurosynth labels |
| `hansen_neurotransmitter_maps` | `@image_vector` | Hansen neurotransmitter map similarity |

## Data extraction

Pulling values out of images, by mask / atlas / parcellation / coordinate.

| Method | From | One-liner |
|---|---|---|
| `apply_mask` | `@image_vector` | Restrict object to voxels in a mask |
| `apply_atlas` | `@image_vector` | Mean-or-pattern expression per atlas region |
| `apply_parcellation` | `@image_vector` | Mean-or-pattern expression per parcel |
| `extract_roi_averages` | both | Average per contiguous region (or mask values) |
| `extract_gray_white_csf` | `@image_vector` | Mean + top-5 components in GM, WM, CSF |
| `extract_measures_batch` | `@fmri_data` | Batch-extract multiple measures at once |
| `get_xyzmm_coordinates` | `@image_vector` | Voxel indices → MNI mm coordinates |
| `normalize_gm_by_wm_csf` | `@fmri_data` | Normalize GM signal by WM/CSF |
| `runRestMetrics` | `@fmri_data` | Resting-state derived metrics |

## Data processing

Workflows for transforming and processing data objects

| `denoise_timeseries_pipeline` | `@fmri_data` | End-to-end denoising for timeseries data |
| `preprocess` | `@image_vector` | Many preprocessing options (filter / outliers / scale) |
| `rescale` | `@fmri_data` | Per-image / per-voxel rescaling |

## Quality control

Diagnosing and cleaning a dataset before analysis.

| Method | From | One-liner |
|---|---|---|
| `descriptives` | `@image_vector` | Print summary stats for the dataset |
| `qc_metrics_second_level` | `@image_vector` | QC metrics across a 2nd-level set |
| `outliers` | `@image_vector` | Detect outlier images |
| `outliers_xval` | `@image_vector` | Outlier detection with cross-validated threshold |
| `windsorize` | `@fmri_data` | Voxel-wise windsorization |
| `sim_data` | `@fmri_data` | Simulate fmri_data for testing |
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
| `create` | `@fmri_data` | Internal helper to assemble fields into an object |
| `unstack_by_condition` | `@image_vector` | Split into sub-objects by a condition vector |
| `validate_object` | `@fmri_data` | Sanity-check internal invariants |
