# Extracting_brain_data_from_regions_of_interest

This topic covers extracting ROI averages and local pattern expression from masks or atlases, including gray/white/CSF summaries and atlas-based parcel extraction.

## Object methods

### fmri_data.extract_roi_averages
Extracts ROI averages or local pattern expression from an fmri_data object using a mask or atlas.
Major options: 'unique_mask_values' vs. 'contiguous_regions', 'pattern_expression', 'nonorm', 'cosine_similarity', 'correlation'.

### image_vector.extract_roi_averages
Extracts ROI averages from an image_vector object; similar to fmri_data version but without pattern-expression mode.
Major options: 'unique_mask_values' vs. 'contiguous_regions', 'noverbose' (see function header).

### image_vector.apply_parcellation
Computes parcel means and optional local pattern expression for atlas/parcellation objects.
Major options: 'pattern_expression', 'correlation', 'norm_mask', 'ignore_missing', 'cosine_similarity'.

### image_vector.apply_mask
Applies a mask or weight map to an image_vector object; can return pattern expression or correlations.
Major options: 'pattern_expression', 'correlation', 'norm_mask', 'ignore_missing', 'invert', 'cosine_similarity'.

### image_vector.extract_gray_white_csf
Extracts global GM/WM/CSF means and top components using canonical tissue masks.
Major options: 'eval' (custom summary function), 'masks' (custom tissue masks).

### atlas.extract_data
Atlas method that wraps apply_parcellation to extract parcel means and local pattern expression from fmri_data.
Major options: same as apply_parcellation, including 'pattern_expression', 'correlation', 'norm_mask', 'ignore_missing', 'cosine_similarity'.

## Stand-alone functions

### extract_image_data (Data_extraction/extract_image_data.m)
Standalone ROI extraction from image files with optional averaging across contiguous regions or unique mask values.
Major options: 'contiguous_regions' (default) or 'unique_mask_values'.

### tor_extract_rois (Data_extraction/tor_extract_rois.m)
Extracts timeseries from clusters in SPM results, optionally fitting a design matrix.
Major options: SPM/VOL/xX inputs for reusing loaded results, and verbosity control via the 4th argument.
