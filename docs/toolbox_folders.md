# CanlabCore folder map

This page lists the directories under `CanlabCore/` and what each one contains. The toolbox is a flat collection of folders rather than packaged MATLAB modules; running `canlab_toolbox_setup` adds them all to your path.

## Class folders

These are MATLAB class folders (`@<classname>/`). Each contains the class definition and the methods you can call on instances of that class. See [Object_methods.md](Object_methods.md) for the cross-cutting documentation index.

| Folder | Class | Purpose |
|---|---|---|
| `@image_vector/` | `image_vector` | Abstract base class for objects that hold image data as `[voxels x images]`. Most other image classes inherit from it. |
| `@fmri_data/` | `fmri_data` | Workhorse class for fMRI / PET / contrast images. Adds `.X`, `.Y`, `.covariates`, prediction, regression, ICA, etc. |
| `@statistic_image/` | `statistic_image` | Stat maps with t / p / effect-size values, with thresholding state. |
| `@atlas/` | `atlas` | Brain atlases / parcellations with probability maps and labels. |
| `@region/` | `region` | Lists of contiguous clusters / ROIs as a unit of analysis. |
| `@fmridisplay/` | `fmridisplay` | Container holding handles for slice montages and surfaces, used to keep figure layers around for blob editing. |
| `@brainpathway/` | `brainpathway` | Connectivity / pathway-modeling object for one subject. |
| `@brainpathway_multisubject/` | `brainpathway_multisubject` | Group-level extension of `brainpathway` for many subjects. |
| `@fmri_timeseries/` | `fmri_timeseries` | Specialized container for raw timeseries data. |
| `@canlab_dataset/` | `canlab_dataset` | Generic subject x variable behavioral / clinical data container. |
| `@fmri_glm_design_matrix/` | `fmri_glm_design_matrix` | Design matrix container for first-level GLMs. |
| `@predictive_model/` | `predictive_model` | Holds artifacts of a fitted multivariate prediction model. |
| `@fmri_mask_image/` | `fmri_mask_image` | Binary mask container (legacy; many newer methods accept `fmri_data` or `image_vector` directly). |

## Function library folders

| Folder | Purpose |
|---|---|
| `canlab_canonical_brains/` | Canonical underlay images, surface meshes, and standard masks (gray-matter, white-matter, ventricles, brainmask). Used by `montage`, `surface`, `extract_gray_white_csf`, etc. |
| `Cifti_plotting/` | Tools to read CIFTI dlabel files and render data on surfaces. Includes `render_cifti_on_brain`, `cifti_struct_2_region_obj`, `make_surface_figure`. |
| `Cluster_contig_region_tools/` | Tools for finding, merging, and labeling contiguous clusters in image volumes. The bridge between voxelwise stat maps and `region`-object outputs. |
| `Data_extraction/` | Loaders and ROI extraction. Includes `load_atlas`, `load_image_set`, `canlab_load_ROI`, `extract_image_data`, `extract_from_rois`. |
| `Data_processing_tools/` | Timeseries filtering, FFT, distance matrices, mean imputation, and similar low-level helpers (e.g., `hpfilter`, `canlab_compute_similarity_matrix`, `canlab_extract_ventricle_wm_timeseries`). |
| `diagnostics/` | Quality-control and image-diagnostic functions: `canlab_qc_metrics1`, `scnlab_norm_check`, intensity histograms, image effect-size maps. |
| `External/` | Vendored third-party toolboxes (matlab_bgl, spider, lasso, BCT, umap, export_fig, boundedline, etc.). Generally treat as read-only. |
| `Filename_tools/` | Filename and BIDS helpers: `canlab_list_files`, `canlab_list_subjects`, `extract_bids_info`, `expand_4d_filenames`, `gunzip_image_names_if_gz`. |
| `fmridisplay_helper_functions/` | Internal helpers used by `@fmridisplay` to lay out, color, and update slice / surface views. |
| `GLM_Batch_tools/` | SPM-driven first- and second-level GLM batch system. Driven by a `DSGN` struct; see `canlab_glm_subject_levels`, `canlab_glm_group_levels`, `canlab_glm_README.txt`. |
| `hewma_utility/` | Hierarchical Exponentially Weighted Moving Average change-point analysis utilities for fMRI timeseries. |
| `HRF_Est_Toolbox4/` | HRF estimation (canonical + temporal / dispersion derivatives, FIR / smooth FIR, IL-model, NL-gamma, splines). |
| `Image_computation_tools/` | Voxelwise computation utilities used by methods (e.g., `image_eval_function`, `image_eval_function_multisubj`, `mean_image`, `scn_num_volumes`, `scn_write_plane`, `tor_global`). |
| `Image_space_tools/` | Reslicing and voxel-space transformations (e.g., `mat2voxel`, `mm2voxel`, `voxel2mm`, `scn_map_image`). |
| `Image_thresholding/` | Statistical thresholding utilities (`FDR`, cluster-extent thresholds, GRF, `clusterSizeMask`). |
| `Index_image_manip_tools/` | The `iimg_*` family — internal index-image / volInfo manipulation that backs many class methods (`iimg_read_img`, `iimg_get_data`, `iimg_reconstruct_vols`, `iimg_reslice`, etc.). |
| `Misc_utilities/` | General-purpose helpers and the canonical `documentation_template.m`. Includes things that don't fit elsewhere — `framewise_displacement`, `canlab_print_legend_text`, etc. |
| `mlpcr/` | Multi-level Principal Component Regression, used as a `predict` algorithm option. |
| `Model_building_tools/` | First-level GLM design construction: `onsets2fmridesign`, `getPredictors`, `create_random_er_design`, `fmri_spline_basis`, basis-set tools. |
| `OptimizeDesign11/` | Genetic-algorithm fMRI experimental-design optimization. |
| `Parcellation_tools/` | Parcellation utilities: cluster solutions, splitting / merging parcels, parcel-wise summary. |
| `peak_coordinates/` | Peak coordinate detection and extraction (`tal_pos_neg`, peak utilities). |
| `Reporting/` | Reporting helpers for analysis output (e.g., `summarize_regression`). |
| `ROI_drawing_tools/` | Interactive ROI drawing utilities. |
| `Sample_datasets/` | Small test fixtures shipped with the toolbox: Wager 2008 EmotionReg contrasts, Pinel localizer timeseries, Jepma single-trial dataset, Woo 2015 BMRK3 pain, Atlas 2012 REMI behavioral. See [sample_datasets.md](sample_datasets.md). |
| `Statistics_tools/` | Statistical functions: `correlation`, `glmfit_*`, classifiers, bootstrap, FDR, mediation helpers, igls, predictor selection, etc. |
| `Unit_tests/` | The `matlab.unittest`-based test suite. Run with `canlab_run_all_tests`. See [Unit_tests/README.md](../CanlabCore/Unit_tests/README.md). |
| `Visualization_functions/` | The plotting library: `barplot_columns`, `addbrain`, `canlab_results_fmridisplay`, `canlab_canonical_brain_surface_cutaways`, color maps, histogram tools, etc. |
| `web_repository_tools/` | Helpers to clone, update, and sync the CANlab GitHub repositories (e.g., `canlab_clone_github_repository`). |

## Top-level files

| File | Purpose |
|---|---|
| `canlab_toolbox_setup.m` | Entry-point script users run from the directory above their CANlab clones; finds and adds CanlabCore + sibling repos to the MATLAB path. |
| `canlab_clone_github_repository.m` | Helper for `canlab_toolbox_setup`; clones missing CANlab repos. |
| `fmri_data_methods.md` | Cross-cutting functional index of `@fmri_data` + `@image_vector` methods (lives at the toolbox root, not in `docs/`). |
| `poorly_documented_functions.txt` | Ranked survey of files whose help blocks differ from the documentation template; used to scope docstring rewrites. |
