# Annotating_images_with_spatial_similarity_including_tools_for_calculating_similiarity_with_resting_state_networks_neurochemical_binding_maps_neurosynth_maps_meta_analysis_maps_and_others

This topic covers tools for annotating neuroimages by spatial similarity to reference maps, including resting-state networks, neurochemical binding maps, Neurosynth topic/term maps, and other meta-analytic or signature mapsets.

## Object methods

### image_vector.image_similarity_plot
Core method for computing similarity between input images and a reference mapset (e.g., Buckner/Yeo resting-state networks, signature libraries, or custom mapsets).
Major options: 'mapset' (named mapsets or custom fmri_data), 'networknames', 'cosine_similarity', 'binary_overlap', 'plotstyle' (polar/wedge), 'average', 'compareGroups', 'nofigure', 'noplot', 'colors', 'dofixrange', 'Error_STD'.

### image_vector.image_similarity_plot_bucknermaps
Convenience wrapper that compares images to Buckner Lab 7-network resting-state maps with polar plotting.
Major options: 'average', plot controls (see function header).

### image_vector.hansen_neurotransmitter_maps
Computes similarity between images and Hansen 2022 PET neurotransmitter binding maps.
Major options: 'cosine_similarity' (or correlation), 'noplot'/'nofigure', 'colors', 'dofixrange', 'average', 'compareGroups'.

### fmri_data.neurosynth_feature_labels
Computes similarity between images and Neurosynth topic/term maps and returns top matching labels.
Major options: 'topics_fi', 'topics_ri', 'images_are_replicates', 'noverbose', 'display_output'.

### fmri_data.neurosynth_lexical_plot
Plots correlations between images and Neurosynth lexical features using neurosynth_feature_labels.
Major options: inherits neurosynth_feature_labels options.

### fmri_data.annotate_binary_results_map
Annotates thresholded (binary) results maps with gradients, resting-state networks, neurochemical maps, and Neurosynth matches.
Major options: see function header; relies on image_similarity_plot and neurosynth_feature_labels.

### fmri_data.annotate_continuous_neuroimage_maps
Annotates continuous-valued maps with resting-state networks, neurochemical maps, and Neurosynth matches.
Major options: see function header; relies on image_similarity_plot and neurosynth_feature_labels.

### image_vector.outliers_xval
Computes similarity of each image to the mean of others using repeated k-fold CV (MAD, cosine, and correlation similarity).
Major options: see function header for cross-validation settings.

## Stand-alone functions

### canlab_pattern_similarity (Statistics_tools/canlab_pattern_similarity.m)
Low-level similarity engine (dot product, cosine, correlation, or binary overlap) used by apply_mask and image_similarity_plot.
Major options: 'dot_product' (default), 'cosine_similarity', 'correlation', 'binary_overlap', 'posterior_overlap', 'ignore_missing', 'treat_zero_as_data', 'exclude_zero_mask_values'.

### plot_correlation_matrix (Visualization_functions/plot_correlation_matrix.m)
Plots similarity or correlation matrices with optional clustering and labeling (useful for mapset similarity summaries).
Major options: 'input_is_r', 'names', 'reorder_by_clustering', 'doplot'.

### distosim (Statistics_tools/Support_functions/distosim.m)
Converts distance/dissimilarity matrices into similarity matrices and validates matrix form.
Major options: none beyond input matrix.

### riverplot (fmri_data/riverplot.m)
Visualizes similarity matrices between two sets of maps with ribbon plots; can compute similarity via image_similarity_plot.
Major options: 'similarity_matrix' (precomputed), 'r'/'dice', thresholding modes, and plotting controls.
