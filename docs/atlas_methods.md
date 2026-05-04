# `atlas` methods, organized by area

The `atlas` class is a CANlab data class for brain atlases and parcellations.
It is a subclass of `image_vector` (via `fmri_data`), so an `atlas` object
behaves like an image-vector object with extra fields and methods specifically
designed for region-labeled / probabilistic parcellations. The integer
parcel-index volume is stored in `.dat` (one integer per voxel),
voxels-by-regions probability maps live in `.probability_maps`, and per-region
text labels (with optional hierarchical levels and long-form descriptions) live
in `.labels`, `.labels_2`...`.labels_5`, and `.label_descriptions`.

In addition to the methods listed below, `atlas` inherits all `fmri_data` and
`image_vector` methods (see `fmri_data_methods.md` and
`image_vector_methods.md`). That means standard operations such as `montage`,
`surface`, `apply_mask`, `write`, `descriptives`, `flip`,
`image_similarity_plot`, `image_math`, `resample_space`, `extract_roi_averages`,
etc. are all available on atlas objects too. Only methods that override or are
unique to `@atlas/` are documented here in detail. Use `methods(my_atlas)` in
MATLAB for the live list; `load_atlas` returns named atlases distributed with
CANlab.

## Properties

`atlas` inherits all `image_vector` and `fmri_data` properties (e.g., `.dat`,
`.volInfo`, `.history`, `.image_names`, `.fullpath`, etc.). The properties
defined on the `atlas` class itself are:

| Property | Description |
|---|---|
| `atlas_name` | Short description or name of the atlas |
| `probability_maps` | Voxels x regions matrix of probability values for each region (sparse) |
| `labels` | Cell array of text strings, one per region (primary labels) |
| `label_descriptions` | Regions x 1 cell array of long-form descriptions |
| `labels_2` | Optional secondary label cell array (e.g., coarser parcellation) |
| `labels_3` | Optional tertiary label cell array |
| `labels_4` | Optional quaternary label cell array |
| `labels_5` | Optional quinary label cell array |
| `references` | String matrix of associated publications |
| `space_description` | Description of atlas space/template (e.g., `MNI152NLin2009cAsym`); used by `render_on_surface` for automatic surface projection |
| `property_descriptions` | Internal cell array of human-readable descriptions for each property (legacy / introspection) |
| `additional_info` | Free-form struct for attaching arbitrary metadata (legacy / extension slot) |

Note: `dat` (inherited from `image_vector`) holds an integer vector with one
parcel index per voxel for atlas objects; this is enforced by
`check_properties`.

## Basic image math and operations

Combining atlases, splitting / merging regions, and reshaping the
parcellation.

| Method | From | One-liner |
|---|---|---|
| `horzcat` | `@atlas` | `[a, b]` operator on atlas objects (sequential `merge_atlases`) |
| `merge_atlases` | `@atlas` | Add regions from one atlas to another, with replace/no-replace options |
| `remove_atlas_region` | `@atlas` | Remove region(s) by name or integer index |
| `reorder_atlas_regions` | `@atlas` | Reorder atlas regions, optionally grouping by label patterns |
| `select_atlas_subset` | `@atlas` | Select a subset of regions by name or integer code |
| `split_atlas_by_hemisphere` | `@atlas` | Divide bilateral regions into separate L and R regions |
| `split_atlas_into_contiguous_regions` | `@atlas` | Split each labeled region into separate contiguous blobs |
| `downsample_parcellation` | `@atlas` | Remap to a coarser nested parcellation using `labels_2`...`labels_5` |
| `probability_maps_to_region_index` | `@atlas` | Rebuild integer index `.dat` from `.probability_maps` (winner-take-all) |
| `assign_vals` | `@atlas` | Assign numeric values to regions and produce an `fmri_data` object |
| `parcel_data2fmri_data` | `@atlas` | Expand parcel-wise values into a voxelwise `fmri_data` object |
| `parcel_stats2statistic_image` | `@atlas` | Expand parcel-wise t / p / dfe into a voxelwise `statistic_image` |

## Display and visualization

| Method | From | One-liner |
|---|---|---|
| `montage` | `@atlas` | Display atlas regions on a standard slice montage (forwards to `region/montage`) |
| `isosurface` | `@atlas` | Render each atlas region as a colored 3-D isosurface |

## Statistics

| Method | From | One-liner |
|---|---|---|
| `atlas_similarity` | `@atlas` | Annotate atlas regions with labels + Dice / coverage stats from a reference atlas |
| `match_atlas_labels` | `@atlas` | Per-region Dice coefficients and best-match labels against another atlas |

## Tables

| Method | From | One-liner |
|---|---|---|
| `label_table` | `@atlas` | Build a MATLAB `table` from the atlas's label / description fields |

## Data extraction

| Method | From | One-liner |
|---|---|---|
| `extract_data` | `@atlas` | Atlas-wise parcel means and local pattern responses from an `fmri_data` object |
| `atlas_get_probability_maps` | `@atlas` | Return probability maps (or indicator maps) as an `fmri_data` object |
| `get_region_volumes` | `@atlas` | Per-region volume in mm^3 and raw voxel counts |
| `num_regions` | `@atlas` | Total regions, regions with data, and any missing region indices |
| `find_closest_region` | `@atlas` | Find the closest atlas region to a given [x y z] mm coordinate |
| `get_regions_at_crosshairs` | `@atlas` | Atlas labels at the current SPM orthviews crosshair |
| `select_regions_near_crosshairs` | `@atlas` | Select regions within X mm of the SPM orthviews crosshair |
| `atlas2region` | `@atlas` | Convert an atlas object to a `region` object |

## Data processing

Workflows for transforming and labeling atlas objects.

| Method | From | One-liner |
|---|---|---|
| `threshold` | `@atlas` | Threshold atlas regions based on values in `.probability_maps` |
| `atlas_add_L_R_to_labels` | `@atlas` | Standardize lateralization suffixes (`_L`/`_R`) on atlas labels |

## Misc utilities

| Method | From | One-liner |
|---|---|---|
| `create` | `@atlas` | Internal helper to populate atlas fields from name/value pairs (used by constructor) |
| `check_properties` | `@atlas` | Validate / cast atlas fields and (optionally) compress the index |

## Inherited methods

`atlas` also inherits all `fmri_data` and `image_vector` methods (see
`fmri_data_methods.md` and `image_vector_methods.md`). Commonly used
inherited methods on atlas objects include `apply_mask`, `apply_atlas`,
`apply_parcellation`, `extract_roi_averages`, `resample_space`, `write`,
`descriptives`, `flip`, `image_similarity_plot`, `image_math`, `surface`,
`render_on_surface`, and `orthviews`.
