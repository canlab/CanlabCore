# Clustering

This topic focuses on identifying, thresholding, and summarizing contiguous voxel clusters in statistical maps, plus utilities for cluster-level statistics, reporting, and visualization.

## Object methods

### statistic_image.threshold
Thresholds statistical images and supports cluster extent correction (GRF) and voxel-wise corrections.
Major options: threshold type ('unc', 'fdr', 'fwe', 'bfr', 'extent'/'cluster_extent', raw range types), 'k' (extent), 'mask', 'noverbose', and residual images/df for FWE.

### statistic_image.multi_threshold
Applies multiple thresholds and tracks cluster persistence across levels.
Major options: threshold vectors and cluster size vectors; supports cluster-based pruning (see function header).

### image_vector.threshold
Thresholds image_vector objects with optional cluster extent constraints.
Major options: threshold type and 'k' (extent) (see function header).

### statistic_image.orthviews
Converts images to clusters for orthviews display using cluster_orthviews.
Major options: passes inputs through to cluster_orthviews (overlay, colors, add, etc.).

## Stand-alone functions

### image2clusters
Menu-driven conversion from a thresholded image to clusters, with optional cross-image masking.
Major options: interactive threshold/extent selection and optional cross-contrast masking (pos/neg/both).

### mask2clusters
Converts a mask image or logical volume into a cluster structure.
Major options: mask input and optional affine (see function header).

### clusters2mask
Converts cluster structures into a mask image or matrix.
Major options: output format controls and mask size (see function header).

### cluster_table
Prints a cluster summary table with optional Talairach/Carmack labels and file output.
Major options: subcluster flag, label flag, 'tal_info', 'talairach', 'writefile', 'noverbose', extra fields to print.

### cluster_local_maxima
Finds local maxima separated by a minimum distance in mm.
Major options: dthresh (distance in mm), verbose, optional output of subcluster assignment.

### merge_clusters / merge_nearby_clusters
Combines clusters into a single structure or merges nearby clusters by distance.
Major options: distance thresholds and merge behaviors (see function headers).

### cluster_intersection / cluster_set_intersection
Computes intersections between cluster sets or between clusters and masks.
Major options: intersection criteria and output selection (see function headers).

### cluster_sigregions / cluster_ttest
Performs cluster-wise t-tests and reports significant regions with optional behavioral regressors.
Major options: pthr, optional behavioral regressor, nonparametric permutation settings (via npm_ttest).

### cluster_names
Assigns or edits cluster labels, with optional overlays and surface plots.
Major options: addflag, overlay, rename, and display toggles.

### cluster_tool
Interactive GUI for cluster editing, subdivision, and reporting.
Major options: GUI-driven; see cluster_tool_* helpers for specific behaviors.
