# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

CanlabCore is a MATLAB toolbox for MRI/fMRI/PET analysis from the Cognitive and Affective Neuroscience Lab (PI: Tor Wager). The core abstraction is a small set of object types that wrap neuroimaging data and provide a consistent, high-level method surface (`plot`, `predict`, `ica`, `threshold`, `apply_atlas`, `montage`, `surface`, ...). There is **no build system, no test runner, and no linter** — code is loaded onto the MATLAB path and exercised interactively or via user scripts.

## Setup and "running" the toolbox

- From a directory **above** the cloned repos, run `canlab_toolbox_setup` in MATLAB. It searches for sibling CANlab repos (CanlabCore, Neuroimaging_Pattern_Masks, CANlab_help_examples, MediationToolbox, RobustToolbox, etc.), adds each to the path with subfolders, and offers to `git clone` any that are missing.
- Required dependencies: MATLAB + Statistics Toolbox + Signal Processing Toolbox + **SPM12** (https://www.fil.ion.ucl.ac.uk/spm/). SPM is used heavily for image I/O (`spm_vol`, `spm_read_vols`, `spm_orthviews`).
- The companion repo `Neuroimaging_Pattern_Masks` (already added as a working directory here) provides the atlases, signatures, and meta-analysis maps that `load_atlas` and `load_image_set` resolve by keyword. Many methods silently depend on its files being on the path.

## Tests

There is no test harness. What exists:
- `CanlabCore/Unit_tests/` — three standalone scripts (`check_roi_extraction.m`, `jackknife_similarity_unit_test.m`, `resampling_pattern_expression_unit_test1.m`). Run by `cd`'ing in MATLAB and calling the function name.
- `@fmri_data/predict_test_suite.m` and `@fmri_data/validate_object.m` — broader sanity checks invoked on a constructed object (e.g. `validate_object(my_fmri_data)`).
- For ad-hoc verification, the canonical pattern is to load a sample dataset (`load_image_set('emotionreg')` or files under `CanlabCore/Sample_datasets/`) and run the method end-to-end.

## Object architecture (the part you must understand to be productive)

Almost everything is built around a single design idea: **brain images are stored flat as a `[voxels × images]` matrix in the object's `.dat` field, with `volInfo` carrying the inverse mapping back to 3-D space.** This lets generic statistical/ML code operate on `.dat` without knowing about neuroimaging, while `reconstruct_image`, `orthviews`, `montage`, `surface`, etc. recover the spatial view on demand.

### Class hierarchy

`image_vector` is the abstract superclass. The classes you will actually instantiate are its subclasses:

- **`fmri_data`** — the workhorse. Holds 4-D fMRI/PET/contrast data plus `.X` (predictors), `.Y` (outcomes), `.covariates`, `.images_per_session`, etc. Most analysis methods (`predict`, `ica`, `regress`, `searchlight`, `mahal`, `preprocess`) live here.
- **`statistic_image`** — t/p/effect-size maps. Knows about thresholds; `threshold(...)` re-thresholds without losing the underlying values.
- **`atlas`** — labeled parcellations. Has `.probability_maps`, `.labels`, `.label_descriptions`, and methods like `select_atlas_subset`, `merge_atlases`, `downsample_parcellation`, `atlas2region`.
- **`fmri_mask_image`** — binary masks (mostly legacy; many newer methods accept a plain `fmri_data` or `image_vector` as a mask).

Other top-level classes (not subclasses of `image_vector`):
- **`region`** — list of contiguous clusters. Produced by `region(statistic_image)` and consumed by `montage`, `table`, `surface`, `extract_data`. The bridge between voxelwise maps and ROI summaries.
- **`fmridisplay`** — a registered handle bag for a montage/surface figure. Workflow is `o2 = canlab_results_fmridisplay(...)`, then `addblobs(o2, region(t))`, `removeblobs(o2)`, `addblobs(o2, ..., 'nolegend')`. The point is that the figure persists; you swap blob layers in/out without re-rendering anatomy.
- **`brainpathway` / `brainpathway_multisubject`** — connectivity / pathway-modeling objects.
- **`canlab_dataset`** — generic subject × variable behavioral/clinical data container with its own `glm`, `mediation`, `scatterplot`, etc.
- **`fmri_glm_design_matrix`**, **`fmri_timeseries`**, **`predictive_model`** — specialized containers for design matrices, raw timeseries, and ML model artifacts.

### MATLAB `@class` directories

Each class lives in `CanlabCore/@<classname>/`. Files in that directory are methods of that class, dispatched via the first argument. The constructor is `@classname/classname.m`. **Adding a method = dropping a `.m` file into the `@class/` folder** with `function out = methodname(obj, ...)`. There is no methods block to edit; `methods(obj)` discovers them from disk. Because `fmri_data`, `statistic_image`, and `atlas` all subclass `image_vector`, methods defined in `@image_vector/` are inherited by all of them — a method only needs to be redefined in a subclass directory if its behavior differs.

### Provenance and "removed" bookkeeping

Two invariants that recur across nearly every method:

1. **`history`** — a cell array of strings appended to by methods. New methods that mutate the object should push a one-line description.
2. **`removed_voxels` / `removed_images`** — when voxels or images are dropped (e.g. `remove_empty`, `apply_mask`), the object shrinks `.dat` and records which rows/columns were removed. `replace_empty(obj)` re-expands `.dat` to the original shape, padded with zeros, so downstream code that needs full-space indexing can rely on it. Many bugs in this codebase historically came from forgetting to call `replace_empty` or `remove_empty` at the right point — when in doubt, call `replace_empty` before reasoning about voxel positions and `remove_empty` before doing math across `.dat` rows.

## Canonical workflows (use these as templates)

```matlab
% Group analysis end-to-end
imgs = load_image_set('emotionreg');         % fmri_data with sample images
plot(imgs); descriptives(imgs);              % QC
t = ttest(imgs);                             % statistic_image
t = threshold(t, .005, 'unc', 'k', 10);      % cluster-extent threshold
r = region(t);                               % region object, one per blob
table(t);                                    % printed/atlas-labeled results
o2 = canlab_results_fmridisplay(t, 'full');  % registered montage+surface figure
montage(r, 'regioncenters', 'colormap');     % per-blob mini-montage
```

```matlab
% ROI extraction against an atlas
atl = load_atlas('canlab2024');              % keyword-resolved atlas
parcel_means = apply_parcellation(imgs, atl); % images x parcels
```

```matlab
% Cross-validated prediction
[cv, stats, optout] = predict(imgs, 'algorithm_name','cv_lassopcr', 'nfolds',5);
```

## Layout

- `CanlabCore/@*/` — the object classes described above.
- `CanlabCore/Statistics_tools/`, `Visualization_functions/`, `Data_processing_tools/`, `Image_thresholding/`, `Model_building_tools/`, `Reporting/` — function libraries called by the class methods. Edit here when a method's logic is not class-specific.
- `CanlabCore/Data_extraction/` — `load_atlas.m`, `load_image_set.m` (keyword resolvers), `extract_*` helpers.
- `CanlabCore/GLM_Batch_tools/` — `canlab_glm_subject_levels` / `canlab_glm_group_levels`, an SPM12-driven first/second-level batch system. Driven by a `DSGN` struct; see `canlab_glm_dsgninfo.txt` and `canlab_glm_README.txt`.
- `CanlabCore/HRF_Est_Toolbox2/` and `HRF_Est_Toolbox4/` — Lindquist-lab HRF estimation (Logit, sFIR, spline, canonical).
- `CanlabCore/OptimizeDesign11/` — genetic-algorithm fMRI design optimization.
- `CanlabCore/Cifti_plotting/`, `Parcellation_tools/`, `Cluster_contig_region_tools/`, `ROI_drawing_tools/`, `Image_space_tools/`, `Image_computation_tools/` — domain-specific helpers.
- `CanlabCore/External/` — vendored third-party toolboxes (`matlab_bgl`, `spider`, `lasso`, `boundedline`, `export_fig`, `BCT`, `umap`, ...). Treat as read-only; don't refactor.
- `CanlabCore/Sample_datasets/` — small datasets used by examples and walkthroughs.
- `CanlabCore/Unit_tests/` — sparse standalone test scripts (see Tests above).
- `nipype/`, `docs/`, `docs_sphinx_old/` — Python wrappers and old docs; rarely touched.

## Conventions worth knowing

- **First argument is always the object** (`function out = method(obj, varargin)`); methods are typically called as `method(obj, ...)` rather than `obj.method(...)`, though both work.
- **`varargin` keyword pairs** are the universal option style. Existing methods use a hand-rolled `for i=1:length(varargin), switch varargin{i}, case 'foo', foo = varargin{i+1};` loop. **New functions should use `inputParser` instead** — see the next section.
- **Many methods accept either a filename, an `fmri_data`, or another `image_vector` subclass** as their "image-like" argument and dispatch via `isa(...)`. Preserve that polymorphism when editing.
- **Spatial alignment is not implicit.** Methods that combine two image objects generally either error or call `resample_space(a, b)` first; if you write a new combiner, do the same — don't assume two objects share `volInfo`.
- **`.asv` files are MATLAB autosave artifacts** and are gitignored; ignore them. A few committed `*_old.m` files are intentional legacy fallbacks (e.g. `region2imagevec_old.m`, `predictive_model_old.m`) — don't delete them without checking callers.

## When writing new functions

These rules apply to new code. Existing code does not need to be retrofitted.

1. **Name stand-alone functions `canlab_<function_name>`.** This namespaces the function so it does not collide with future external toolboxes the user may add to their path. Class methods (files inside `@class/`) are exempt — they're already namespaced by the class.

2. **Match the documentation format in `CanlabCore/Misc_utilities/documentation_template.m`.** That template defines the section ordering (Usage, Inputs, Outputs, Examples, References, etc.) and comment style readthedocs expects. Open it before writing the help block; copy the structure rather than improvising.

3. **Use `inputParser` for variable input arguments**, following the `INPUT PARSER TEMPLATE` section of `documentation_template.m`. Retain the explanatory comments inside that block — they're a teaching scaffold for future readers, not noise. Implement the `'plot'`, `'verbose'`, `'doplot'`, and `'doverbose'` parameters whenever the function has plotting or chatter that the caller might want to suppress.

4. **Include a runnable example in the help block** that loads or creates a test dataset (e.g. `load_image_set('emotionreg')`, `sim_data`, or a synthetic array) and demonstrates the function with a few of the most common options. The example should be copy-pasteable: someone with CanlabCore on their path should be able to highlight it and run it.

## Documentation pointers

- Function-by-function reference: https://canlabcore.readthedocs.org/en/latest/
- Walkthroughs and batch-script examples (the best way to learn the API): https://github.com/canlab/CANlab_help_examples
- Lab landing page: https://canlab.github.io
