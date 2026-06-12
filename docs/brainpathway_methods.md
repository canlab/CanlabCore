# `brainpathway` methods, organized by area

`brainpathway` is a handle-class container for brain connectivity and
pathway analyses. It stores neuroimaging data at multiple spatial scales
(voxel, node, region, network, partition) and a parcellation
(`region_atlas`, an `atlas` object) that ties them together. Listeners
make it reactive: assigning new `voxel_dat`, `region_dat`, `node_weights`,
or `connectivity_properties` automatically recomputes downstream region
averages, node responses, and connectivity matrices. Connectivity is
estimated by a user-pluggable function handle (default `@corr`) and is
stored in `obj.connectivity.regions` and `obj.connectivity.nodes`.

`brainpathway_multisubject` is a subclass of `brainpathway` for
group-level analyses; it adds a subject dimension to connectivity and
related fields and a `subject_metadata` table. See the dedicated section
below.

Type `methods(my_obj)` in MATLAB for the live list on any instance.

## Properties

| Property | Description |
|---|---|
| `region_atlas` | An `atlas` object defining k regions used for all extractions and connectivity |
| `voxel_dat` | `[voxels x images]` data (single-precision); cell-of-subjects in `brainpathway_multisubject` |
| `node_dat` | `[images x nodes]` aggregated node responses |
| `region_dat` | `[images x regions]` region-average data |
| `network_dat` | `[images x networks]` network-level data |
| `partition_dat` | `[images x partitions]` partition-level data |
| `node_weights` | `{1 x n}` cell of pattern weights (voxels x nodes) per region |
| `node_labels` | `{1 x n}` cell of node names |
| `node_clusters` | `int32` vector of cluster assignments for each node |
| `node_cluster_labels` | Cell of cluster names |
| `region_indx_for_nodes` | `int32` vector mapping each node to its region index |
| `connectivity` | Struct with `.regions` and `.nodes`, each holding `r`, `p`, and within/between summaries |
| `graphstruct` | Struct with `within_network_degree`, `between_network_degree` |
| `graph_properties` | Struct with `regions` and `nodes` tables of graph-theoretic metrics |
| `connectivity_properties` | Struct with `c_fun_han` (e.g. `@corr`) and `c_fun_arguments`; setting it triggers re-estimation |
| `connections_apriori` | Logical `k x k x n` matrices of a priori connections per network |
| `additional_info` | Free-form struct for user metadata |
| `listeners` | Property listeners that drive automatic recalculation |
| `verbose` | Verbose-output flag |
| `data_quality` | Struct with `tSNR`, `tSTD`, `median_corr` per region |

## Basic operations

| Method | From | One-liner |
|---|---|---|
| `copy` | `@brainpathway` | Deep copy of a handle-class brainpathway object |
| `select_atlas_subset` | `@brainpathway` | Restrict object to a subset of atlas regions, propagating to data |
| `attach_voxel_data` | `@brainpathway` | Resample an `fmri_data` object into the atlas space and attach as `voxel_dat` |
| `nan2zero` | `@brainpathway` | Replace NaN-containing region columns with zeros |
| `find_node_indices` | `@brainpathway` | Resolve names/integers to a logical index of nodes |
| `reorder_regions` | `@brainpathway` | Reorder atlas regions by index, cluster, or label-group order |
| `reorder_regions_by_node_cluster` | `@brainpathway` | Convenience reorder using `node_clusters` |

## Display and visualization

| Method | From | One-liner |
|---|---|---|
| `plot_connectivity` | `@brainpathway` | Plot region- and node-level connectivity matrices with optional partitions |

## Statistics and connectivity

| Method | From | One-liner |
|---|---|---|
| `seed_connectivity` | `@brainpathway` | Voxel-/region-/node-wise correlation maps with one or more seed regions |
| `threshold_connectivity` | `@brainpathway` | Threshold connectivity matrices (FDR / uncorrected / Bonferroni) |
| `degree_calc` | `@brainpathway` | Within-, between-, and total-cluster degree from `node_clusters` |
| `cluster_regions` | `@brainpathway` | Ward clustering of regions on `region_dat`; sets `node_clusters` |
| `cluster_region_subset_by_connectivity` | `@brainpathway` | Sub-cluster a target group based on cross-group connectivity profiles |
| `cluster_voxels` | `@brainpathway` | Correlate every voxel with every region average and cluster voxels |

## Internal listeners and helpers

These are static methods invoked by listeners. You normally do not call
them directly; assigning to the relevant property triggers them.

| Method | From | One-liner |
|---|---|---|
| `intialize_nodes` | `@brainpathway` | Create one node per region with uniform voxel weights |
| `update_region_data` | `@brainpathway` | Recompute region averages from `voxel_dat` |
| `update_node_data` | `@brainpathway` | Recompute node responses from `voxel_dat` and `node_weights` |
| `update_region_connectivity` | `@brainpathway` | Recompute `connectivity.regions` from `region_dat` |
| `update_node_connectivity` | `@brainpathway` | Recompute `connectivity.nodes` from `node_dat` |

## Misc utilities

| Method | From | One-liner |
|---|---|---|
| `brainpathway2fmri_data` | `@brainpathway` | Map region/voxel-level values from a brainpathway field into an `fmri_data` object for visualization |

## brainpathway_multisubject (group-level extension)

`brainpathway_multisubject` inherits all of the above. The key difference
is that connectivity matrices and per-subject region data carry an extra
subject dimension (e.g., `connectivity.regions.r` is `k x k x n`). The
class adds a `subject_metadata` table and `HDIs` (hub-disruption index)
struct, and disables the auto-update listeners after construction.

### Additional properties

| Property | Description |
|---|---|
| `subject_metadata` | Table of subject-level metadata (ID, motion, age, condition, etc.) |
| `HDIs` | Struct with `regions` and `nodes` tables of hub-disruption indices |

### Methods specific to `brainpathway_multisubject`

| Method | From | One-liner |
|---|---|---|
| `add_subject` | `@brainpathway_multisubject` | Append a single-subject `brainpathway` to the group object |
| `get_wh_subjects` | `@brainpathway_multisubject` | Subset the group object by a logical subject mask |
| `flatten_conn_matrices` | `@brainpathway_multisubject` | 3-D `k x k x n` matrices to subjects-by-edges 2-D matrix |
| `bct_toolbox_undirected_graph_metrics` | `@brainpathway_multisubject` | Per-subject BCT graph metrics into `graph_properties.regions` |
| `compute_HDIs` | `@brainpathway_multisubject` | Hub-disruption indices vs. a reference group (Achard 2012) |
| `ISC` | `@brainpathway_multisubject` | Brain-wide inter-subject correlation across parcel means |
| `ttest` | `@brainpathway_multisubject` | One- or two-sample edge-wise t-test across subjects with FDR |
| `qcfc` | `@brainpathway_multisubject` | QC-FC correlations between motion (or other QC) and edge strength |
| `plot_qcfc` | `@brainpathway_multisubject` | Plot QC-FC correlations and (optionally) partial-out covariates |
| `plot_connectivity` | `@brainpathway_multisubject` | Group-mean region/node connectivity plus within/between summaries |
| `update_region_connectivity` | `@brainpathway_multisubject` | Listener stub overriding the superclass behavior |
