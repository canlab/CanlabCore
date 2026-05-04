# `region` methods, organized by area

This is a functional index of methods available on a `region` object.
A `region` is a class for groups of voxels defined as anatomical or
functional regions (contiguous voxels above a threshold in an analysis,
parcels from an atlas, etc.). The class replaces the older `clusters`
struct from the scnlab toolbox; the class definition adds error
checking and structured methods (visualization, tabling, conversion to
other image types).

`region` is a standalone class (not a subclass of `image_vector`), so
it does NOT inherit `image_vector` methods. Instead, several region
methods convert to and from `image_vector` / `fmri_data` / `atlas` so
that those classes' tooling can be applied. Type `methods(my_obj)` in
MATLAB for the live list on any instance.

## Properties

| Property | Description |
|---|---|
| `title` | Title of the region object |
| `shorttitle` | Short title (used as label) |
| `descrip1` | Description string 1 |
| `descrip2` | Description string 2 (lists key methods) |
| `XYZ` | [3 x voxels] voxel coordinates within the region |
| `XYZmm` | [3 x voxels] world (mm) coordinates within the region |
| `val` | Per-voxel values (typically the input mask values) |
| `val_descrip` | Description of values stored in `val` |
| `Z` | Per-voxel max-stat / Z-scores (legacy field, still widely used) |
| `Z_descrip` | Description of values stored in `Z` (legacy) |
| `threshold` | (legacy) Threshold associated with the region, kept for compatibility |
| `voxSize` | Voxel size in mm (3-vector from the affine) |
| `M` | Affine matrix (4x4) inherited from the source mask |
| `dim` | Image dimensions (3-vector) |
| `numVox` | Number of voxels in the region |
| `numpeaks` | Number of peaks within the region |
| `center` | Voxel-coordinate center of mass |
| `mm_center` | mm-coordinate center of mass |
| `timeseries` | Per-region timeseries (if extracted) |
| `contrastdata` | Per-region contrast data (if extracted) |
| `dat` | Per-region averaged data (images x 1 vector when extracted) |
| `all_data` | Full per-voxel data extracted from a source image set |
| `source_images` | Filenames of images data were extracted from |
| `custom_info1` | Custom info slot 1 (often the source mask filename) |
| `custom_info1_descrip` | Description for `custom_info1` |
| `custom_info2` | Custom info slot 2 |
| `custom_info2_descrip` | Description for `custom_info2` |

## Basic image math and operations

| Method | From | One-liner |
|---|---|---|
| `merge` | `@region` | Merge two or more regions into one, combining fields appropriately |
| `posneg_separate` | `@region` | Split regions into positive- and negative-valued sub-regions |
| `subdivide_by_atlas` | `@region` | Subdivide each blob by anatomical atlas parcels |
| `subdivide_by_local_max` | `@region` | Subdivide regions by local peak Z-score / maxima |
| `reparse_continguous` | `@region` | Re-define regions based on contiguous blobs (note: spelling is `continguous` in source) |
| `select_coordinates_near_regions` | `@region` | Filter MNI mm coordinates by minimum distance to the region object |

## Display and visualization

| Method | From | One-liner |
|---|---|---|
| `montage` | `@region` | Slice montage of region object on canonical anatomy |
| `orthviews` | `@region` | SPM orthviews via `cluster_orthviews` with hot/cool colormap |
| `surface` | `@region` | Render region blobs on cutaway / canonical surfaces |
| `isosurface` | `@region` | One 3-D isosurface per region, with optional L/R color matching |
| `labelled_surface` | `@region` | Transparent isosurfaces with centroid labels and text annotations |
| `match_colors_left_right` | `@region` | Assign matched colors to symmetric L/R regions |

## Tables

| Method | From | One-liner |
|---|---|---|
| `table` | `@region` | Print and return a labeled table of all regions (positive and negative) |
| `table_simple` | `@region` | Simple one-row-per-region results table (useful after atlas subdivision) |
| `table_of_atlas_regions_covered` | `@region` | Tabulate atlas parcels covered by the regions |
| `ttest_table_by_condition` | `@region` | Per-region means + stats tables for an `fmri_data` object across conditions |
| `autolabel_regions_using_atlas` | `@region` | Auto-label regions using an atlas object |

## Data extraction

| Method | From | One-liner |
|---|---|---|
| `extract_data` | `@region` | Extract data and apply local patterns from an image_vector for the region's voxels |
| `check_extracted_data` | `@region` | Verify region-average data via re-extraction from source images |

## Misc utilities

Conversion to and from other CANlab object classes, and basic predicates.

| Method | From | One-liner |
|---|---|---|
| `region2atlas` | `@region` | Convert a region object to an atlas object |
| `region2fmri_data` | `@region` | Convert a region object to an fmri_data object |
| `region2imagevec` | `@region` | Convert a region object to an image_vector object |
| `region2struct` | `@region` | Convert a region object to a plain MATLAB struct array (legacy compatibility) |
| `isempty` | `@region` | True iff region array is empty or its first XYZ is empty |
