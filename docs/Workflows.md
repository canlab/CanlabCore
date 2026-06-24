# CANlab workflows

A **workflow** is an end-to-end recipe that chains several CanlabCore methods together to accomplish a common analysis goal. Where the per-class method pages ([`fmri_data`](fmri_data_methods.md), [`atlas`](atlas_methods.md), [`region`](region_methods.md), …) document *individual* methods, a workflow shows how those pieces fit together — and which method to reach for when several can do a similar job.

Each workflow comes in two complementary parts:

- a **roadmap** — a conceptual overview that names the available approaches, explains how they relate, and helps you choose the right one; and
- a **walkthrough** — a runnable, didactic guide with copy-pasteable sample code on built-in datasets, including the figures it produces.

Start with the roadmap to orient yourself, then follow the walkthrough to run it.

## Available workflows

| Workflow | What it does | Roadmap (overview) | Walkthrough (code) |
|---|---|---|---|
| **ROI / atlas data extraction** | Pull region-of-interest, pattern, parcel, tissue-compartment, and sphere/coordinate summaries out of brain images, then visualize them (bar plots, line plots, multi-subject slope plots). | [ROI extraction roadmap](workflows/ROI_extraction_methods_roadmap.md) | [Extract & visualize ROI data — how-to](workflows/extract_roi_data_howto.md) |
| **First-level fMRI time-series modeling** (`glm_map`) | Build a single-subject event-related model from onsets (FSL tables, SPM-style cells, or an `SPM.mat`), choose a basis set, add nuisance covariates and contrasts, screen the design (VIF/cVIF, efficiency, high-pass filter), simulate data, fit (AR errors), threshold, and visualize. | [First-level roadmap](workflows/glm_map_first_level_roadmap.md) | [First-level how-to](workflows/glm_map_first_level_howto.md) |
| **Second-level fMRI group analysis** (`glm_map`) | Group regression on contrast images: build the object via `fmri_data.regress` or the `glm_map` estimator API, handle outliers and WM/CSF nuisance signals (`normalize_gm_by_wm_csf`), fit OLS or robust, screen the design, threshold, and visualize. | [Second-level roadmap](workflows/glm_map_second_level_roadmap.md) | [Second-level how-to](workflows/glm_map_second_level_howto.md) |

*More workflows will be added here over time.*

The two `glm_map` workflows share the [`glm_map` object page](glm_map_methods.md) (full method reference) and the [`fmri_glm_design_matrix`](fmri_glm_design_matrix_methods.md) design-building object.

## Visualizing results

Most workflows end at a thresholded `statistic_image` or a `region` object. The same map can be rendered several ways — pick by output medium:

**Static figures (MATLAB):**

- **`montage(t)`** — slice montage on a canonical anatomical underlay; the workhorse for figures and reports.
- **`surface(t)`** — render activation on 3-D cortical surfaces. Style presets include `'foursurfaces_hcp'` (lateral + medial views of both hemispheres with brainstem, on HCP pial surfaces), `'inflated'`, and cutaways; `isosurface` gives 3-D blobs.
- **[`canlab_results_fmridisplay(t)`](individual_functions/canlab_results_fmridisplay.md)** — pre-built montage + surface scaffolds (`'compact2'`, `'full'`, …) returning a registered `fmridisplay` whose blob layers you can swap without re-rendering the anatomy.
- **`table(t)` / `region(t)`** — atlas-labeled results table and the `region` objects behind it.

**Interactive viewers (point-and-click):**

- **`canlab_orthviews(t)`** — SPM-style three-plane viewer in MATLAB. Click/drag to navigate; the bottom strip names the **atlas region under the crosshair** (attach any atlas with `canlab_orthviews('AddAtlasLabel', atl)`).
- **[`canlab_niivue(t)`](canlab_niivue_guide.md)** — a portable, self-contained `.html` web viewer (NiiVue) with colormap/threshold/opacity controls, a crosshair coordinate + value readout, and an atlas region readout that outlines the region under the crosshair. Email it or embed it in an HTML report. **`orthviews_niivue(t)`** is a one-liner shortcut that writes the page to a temp folder and opens it in your browser.

See [Visualizing images and results](Object_methods.md#visualizing-images-and-results) for the full set of visualization entry points.
