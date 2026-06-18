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

*More workflows will be added here over time.*
