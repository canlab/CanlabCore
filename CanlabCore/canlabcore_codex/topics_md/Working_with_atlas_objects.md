# Working_with_atlas_objects

This topic summarizes the atlas class and common operations for selecting, merging, splitting, labeling, and converting atlas regions to other CANlab object types.

## Object methods

### atlas (constructor)
Creates an atlas object from an integer-labeled image or probabilistic maps.
Major options: 'mask' (mask image), 'noverbose', 'sample2mask' (resample images to mask space vs. mask to image space).

### select_atlas_subset
Selects regions by label or integer index, optionally collapsing them into a single mask.
Major options: 'flatten', 'labels_2' (or any label field), 'exact', 'regexp', 'deterministic'/'mostprob', 'conditionally_ind'.

### atlas2region
Converts atlas objects to region objects using label indices or probability maps.
Major options: 'use_probabilities', 'nocontiguous'/'unique_mask_values', 'noverbose'.

### merge_atlases
Adds regions from one atlas into another with control over voxel replacement.
Major options: region selection args for select_atlas_subset, 'noreplace', 'always_replace', 'noverbose'.

### split_atlas_by_hemisphere
Splits bilateral regions into left/right labels using x-coordinates.
Major options: none (behavior is based on midline x = 0 convention).

### split_atlas_into_contiguous_regions
Splits multi-blob regions into separate contiguous regions.
Major options: none (uses contiguity in voxel space).

### atlas_similarity
Annotates atlas regions with labels and overlap metrics from a reference atlas.
Major options: none; outputs include dice, coverage, and label tables.

### get_region_volumes
Returns region volumes and voxel counts for each atlas label.
Major options: none.

### atlas_get_probability_maps
Returns probability maps (or indicator maps) as an fmri_data object.
Major options: none.

### parcel_stats2statistic_image / parcel_data2fmri_data
Map parcel-level statistics or values back to a voxelwise image in atlas space.
Major options: parcel stats and optional significance inputs (see function headers).

### label_table
Returns a table of atlas region labels and metadata.
Major options: none.

### isosurface / montage
Visualization methods for atlas surfaces and slice montages.
Major options: display and coloring options (see function headers).

## Stand-alone functions

### parcel_images / parcel_clusters (Parcellation_tools)
Utilities for building and inspecting parcellations from image data.
Major options: clustering parameters and output controls (see function headers).

### parcel_complete_sets
Organizes parcel sets and metadata for downstream analyses.
Major options: see function header.

### parcel_cl_nmds / parcel_cl_nmds_plots
Runs NMDS on parcel-level data and generates diagnostic plots.
Major options: NMDS parameters and plotting controls (see function headers).
