# Thresholding_images

This topic covers statistical and raw-value thresholding for neuroimaging maps, including FDR/FWE correction, cluster-extent thresholds, and utilities for preparing residuals for smoothness estimation.

## Object methods

### statistic_image

#### threshold
Statistical thresholding for statistic_image objects with reversible masking and multiple correction types.
Major options: thresh_type ('fdr', 'fwe', 'bfr', 'unc', 'extent'/'cluster_extent', 'raw-between', 'raw-outside'), 'k' (cluster extent), 'mask' (ROI/brain mask), 'df' (for FWE), 'noverbose'.

#### multi_threshold
Applies multiple thresholds for visualization and returns positive/negative cluster sets at each level.
Major options: 'thresh' (p-value vector), 'sizethresh' (cluster sizes), 'poscolors', 'negcolors', 'nodisplay', 'o2' (reuse fmridisplay), 'writestats', 'noplot', 'wh_montages'.

### image_vector

#### threshold
Raw-value thresholding for image_vector/fmri_data objects (irreversible), with optional extent pruning.
Major options: thresh_type ('raw-between'/'raw-outside'), 'k' (cluster extent), 'trim_mask', 'noverbose'.

## Stand-alone functions

### FDR (Image_thresholding/FDR.m)
Computes FDR thresholds for a vector of p-values using independence/positive dependence and a nonparametric variant.
Major options: none beyond p-values and q.

### threshold_imgs (Image_thresholding/threshold_imgs.m)
Thresholds image files by height and optional cluster extent, saving filtered images to disk.
Major options: k (extent threshold), direction ('pos', 'neg', 'both').

### clusterSizeMask (Image_thresholding/clusterSizeMask.m)
Creates a binary mask of clusters that meet an extent threshold given a height-thresholded map.
Major options: sizeThresh and height_mask inputs.

### cl_ext_make_resid (Image_thresholding/cl_ext_make_resid.m)
Creates residual images and a mask for smoothness estimation in cluster-extent correction workflows.
Major options: 'mask' (explicit mask image), 'outputdir'.

### cl_ext_spm_grf (Image_thresholding/cl_ext_spm_grf.m)
Estimates cluster-extent thresholds using SPM GRF correction and residual images.
Major options: 'doplot', 'twotail'.

### cl_ext_3dClustSim (Image_thresholding/cl_ext_3dClustSim.m)
Estimates cluster-extent thresholds using AFNI 3dClustSim and residual images.
Major options: 'iter' (simulation iterations), 'twotail', 'fwhm' (manual smoothness).

### cl_ext_spm_spm (Image_thresholding/cl_ext_spm_spm.m)
Modified SPM GLM estimation used to retain residual images for cluster-extent workflows.
Major options: uses SPM struct inputs and masking fields (see function header).
