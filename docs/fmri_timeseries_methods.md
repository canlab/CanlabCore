# `fmri_timeseries` methods, organized by area

`fmri_timeseries` is a subclass of `fmri_data` (and therefore of
`image_vector`) specialized for 4-D time series data. It inherits all
`fmri_data` and `image_vector` properties and methods (data storage in
`.dat` as `[voxels x images]`, `.volInfo`, masks, history, etc.) and adds
a few timeseries-specific fields: a TR, an embedded
`fmri_glm_design_matrix` for first-level GLM specification, and a
processing-status table that tracks slice-timing, realignment,
denoising, filtering, normalization, and smoothing state.

Because the class inherits from `fmri_data`, the vast majority of methods
relevant to a timeseries (preprocess, regress, predict, ica, mahal,
montage, plot, write, etc.) live on the parent classes; see
[`fmri_data` methods](fmri_data_methods.md). The table below lists only
methods defined on `@fmri_timeseries` itself. Type
`methods(my_obj)` in MATLAB for the live list on any instance.

## Properties

| Property | Description |
|---|---|
| `glm_design_obj` | Embedded `fmri_glm_design_matrix` for first-level model specification |
| `processing_status_table` | Table tracking slice-timing, realignment, denoising, filter cutoffs, normalization template, smoothing FWHM |
| `TR` | Repetition time in seconds |

Inherited from `fmri_data` / `image_vector`: `dat`, `volInfo`, `mask`,
`fullpath`, `image_names`, `removed_images`, `removed_voxels`, `history`,
`X`, `Y`, `covariates`, `images_per_session`, `additional_info`, etc.

## Methods specific to `@fmri_timeseries`

| Method | From | One-liner |
|---|---|---|
| `fmri_timeseries` | `@fmri_timeseries` | Constructor: load images with a mask, attach TR and a default GLM design |
| `filloutliers` | `@fmri_timeseries` | Per-voxel moving-median outlier replacement with spline interpolation (60-s window) |

For everything else (resampling, preprocessing, regression, prediction,
extraction, display, I/O), see the inherited methods documented in
[`fmri_data` methods](fmri_data_methods.md).
