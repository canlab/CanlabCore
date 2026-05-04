# Object methods in CanlabCore

CanlabCore is built around a small set of MATLAB classes that wrap neuroimaging data. Each class has **properties** (the data fields it stores — image values, masks, design matrices, region labels, etc.) and **methods** (the things you can do with that data — plot it, threshold it, run a t-test, render it on a brain surface, save it to disk). This page is the entry point to per-class documentation.

The design philosophy is interactive analysis with simple commands. A typical group analysis is a handful of one-liners:

```matlab
imgs = load_image_set('emotionreg');     % fmri_data with 30 images
plot(imgs);                              % QC summary
t    = ttest(imgs);                      % voxelwise t-test -> statistic_image
t    = threshold(t, 0.005, 'unc');       % re-threshold
r    = region(t);                        % connect blobs -> region object
table(r);                                % atlas-labeled results table
montage(r);                              % brain montage
```

Most user-facing image classes inherit from a common abstract base, `image_vector`, which stores image data as a flat `[voxels x images]` matrix in a `.dat` field with the inverse mapping back to 3-D space in `.volInfo`. This is what lets generic statistical / ML code operate on `.dat` while the class methods handle the spatial reconstruction transparently.

## Tutorials and walkthroughs

The fastest way to learn the toolbox is by example.

- **[Walkthroughs](https://canlab.github.io/walkthroughs/)** — step-by-step analysis tutorials with code.
- **[Tutorials](https://canlab.github.io/tutorials/)** — longer-form tutorials.
- **[CANlab_help_examples](https://github.com/canlab/CANlab_help_examples)** repository — runnable MATLAB scripts (`example_help_files/`) and HTML output with figures (`published_html/`), plus the second-level batch script system.
- **[canlab.github.io](https://canlab.github.io/)** — top-level entry point with Setup, Repositories, and Interactive fMRI sections.

## Class hierarchy

```
image_vector  (abstract base; rarely used directly)
├── fmri_data           generic image data + .X, .Y, covariates
├── statistic_image     stat maps with t / p / sig
├── atlas               labeled parcellation with probability maps
└── fmri_mask_image     binary mask (legacy)

region                  list of contiguous clusters
fmridisplay             figure-handle container for layered brain montages
brainpathway            connectivity / pathway model (one subject)
brainpathway_multisubject   group-level extension of brainpathway
fmri_timeseries         specialized container for raw timeseries
canlab_dataset          subject x variable behavioral / clinical data
fmri_glm_design_matrix  first-level GLM design matrix
predictive_model        artifacts of a fitted multivariate prediction model
```

## Object classes

Listed in roughly the order most users encounter them. Click a class name for the full per-class page (intro, properties, methods grouped by category).

| Class | Description |
|---|---|
| **[`fmri_data`](fmri_data_methods.md)** | The workhorse. Holds fMRI / PET / contrast images plus optional predictor matrix `.X`, outcome vector `.Y`, covariates, and metadata. Most analysis methods (`predict`, `regress`, `ica`, `searchlight`, `ttest`, `signtest`) live here. |
| **[`image_vector`](image_vector_methods.md)** | Abstract superclass. You rarely create one directly, but most of the methods you call on an `fmri_data`, `statistic_image`, or `atlas` are inherited from here (`apply_mask`, `resample_space`, `montage`, `surface`, `extract_roi_averages`, etc.). |
| **[`statistic_image`](statistic_image_methods.md)** | Stat maps (t / p / effect-size) with thresholding state. Produced by `ttest`, `regress`, etc. The `threshold` method re-thresholds without losing the underlying values. |
| **[`atlas`](atlas_methods.md)** | Brain atlases / parcellations. Has `.probability_maps`, `.labels`, `.label_descriptions`, `.references`. Methods include `select_atlas_subset`, `merge_atlases`, `downsample_parcellation`, `atlas2region`, `apply_atlas`. Use `load_atlas` to load by keyword. |
| **[`region`](region_methods.md)** | List of contiguous clusters / ROIs as a unit of analysis. Produced by `region(t)` from a thresholded `statistic_image`. Consumed by `montage`, `table`, `surface`, `extract_data`. |
| **[`fmridisplay`](fmridisplay_methods.md)** | Container holding figure handles for slice montages and surfaces. Built by `canlab_results_fmridisplay`; lets you swap blob layers in / out without re-rendering the anatomy underneath. |
| **[`brainpathway`](brainpathway_methods.md)** | Connectivity / pathway-modeling object for one subject. The `brainpathway_multisubject` extension is documented on the same page. |
| **[`fmri_timeseries`](fmri_timeseries_methods.md)** | Specialized container for raw timeseries data. |
| **[`canlab_dataset`](canlab_dataset_methods.md)** | Generic subject x variable behavioral / clinical data container with its own `glm`, `mediation`, `scatterplot`, `get_var`, `add_vars` methods. |
| **[`fmri_glm_design_matrix`](fmri_glm_design_matrix_methods.md)** | Holds GLM design matrices (X) for first-level fMRI analyses. Methods like `build`, `add`, `replace_basis_set`. |
| **[`predictive_model`](predictive_model_methods.md)** | Holds artifacts of a fitted multivariate prediction model — cross-validated predictions, weight maps, performance summaries. |
| **[`fmri_mask_image`](fmri_mask_image_methods.md)** | Binary mask container. Mostly legacy; many newer methods accept an `fmri_data` or plain `image_vector` as a mask directly. |

## Cross-cutting topics

- **[Recasting (converting) between object types](recasting_objects.md)** — `region2fmri_data`, `atlas2region`, `region(t)`, etc., and when to call `replace_empty` before converting.
- **[fmri_data_methods.md](fmri_data_methods.md)** is also the cross-cutting *functional* index of `@fmri_data` + `@image_vector` methods, organized by area (basic math / display / resampling / statistics / multivariate prediction / tables / annotation / data extraction / data processing / quality control / misc utilities). The same category structure is used in every per-class page.
- **[Atlases, regions, and patterns](atlases_regions_and_patterns.md)** — registry of available atlases, named regions, and multivariate signature patterns, with paper citations.
- **[Sample datasets](sample_datasets.md)** — the small datasets that ship with CanlabCore plus the `load_image_set` keyword registry, with paper citations.
- **[Toolbox folder map](toolbox_folders.md)** — what lives in each subfolder of `CanlabCore/`.

## Conventions

- **State management.** Many methods rely on the "removed voxels / removed images" bookkeeping. If you are manipulating `.dat` directly, call `replace_empty(obj)` to expand to the full padded voxel space before reasoning about voxel positions, and `remove_empty(obj)` before doing math across `.dat` rows. See the relevant section of [`fmri_data_methods.md`](fmri_data_methods.md).
- **Provenance.** The `.history` cell array tracks transformations applied to an object. Many methods append to it automatically.
- **Polymorphism.** Most "image-like" arguments accept either a filename, an `fmri_data`, or any `image_vector` subclass and dispatch with `isa(...)`. Spaces are reconciled with `resample_space` when objects don't already match.
