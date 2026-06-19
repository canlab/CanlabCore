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
