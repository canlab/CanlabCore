# Changelog

All notable changes to CanlabCore are documented in this file.

## [Unreleased]

### Visualization: stateful `fmridisplay` (display overhaul, Phase 1–2)

- **`fmridisplay` is now a handle class** (`classdef fmridisplay < handle`). The
  value-style call contract is fully preserved — `o2 = addblobs(o2, ...)` and every
  existing call site keep working (a full audit found no copy-then-reuse aliasing).
  The object instance is now a single source of truth that figures and a controller
  widget can reference.
- **Blob layers retain their source data and render options**, so they can be
  re-rendered in place. New `@fmridisplay` methods: `rethreshold`, `set_colormap`,
  `set_opacity`, and `refresh` (live re-threshold / colormap / opacity on montages).
- **`controller(obj)`** opens a `uifigure` control panel bound to the instance
  (per-layer opacity / colormap / threshold / visibility), the MATLAB-side analog of
  the `canlab_niivue` web control panel.
- **Fix:** a bare `'colormap'` flag (e.g. `montage(r, o2, 'colormap')`) no longer
  errors in `render_blobs`; `addblobs` strips it before forwarding.
- New tests: `Unit_tests/image_vector/canlab_test_fmridisplay_handle.m`.

## [v2.1.0] - 2026-05-18

This release is dominated by infrastructure work: a real automated test harness,
continuous integration, and a structured documentation set. There are no breaking
API changes for typical users; a handful of long-deprecated and legacy files were
removed.

### Highlights

- **First automated test suite and CI.** `matlab.unittest`-based harness with ~70
  tests, run on every push/PR via GitHub Actions and nightly against the
  walkthrough scripts.
- **Structured documentation set under `docs/`.** Per-class method indexes, 35
  per-function "code map" diagrams, 59 per-function help pages with runnable
  Quick examples and sample output PNGs, and curated atlas / sample-dataset
  pages with DOI citations.
- **136 docstrings reformatted** to the sphinx / read-the-docs format used by
  canlabcore.readthedocs.org.
- **Forward-compatible SPM handling.** 26 SPM version-check sites in 21 files no
  longer enumerate accepted versions; SPM25, SPM26 and later now work without
  code changes.
- **Substantial cleanup.** `HRF_Est_Toolbox2/`, `docs_sphinx_old/`, 35
  `_old`/`_backup`/tmp files, and 18 `.asv` autosave files removed.

### Added — Testing infrastructure

- `CanlabCore/Unit_tests/` reorganized into per-class subfolders (`fmri_data/`,
  `image_vector/`, `statistic_image/`, `atlas/`, `region/`, `walkthroughs/`,
  `workflows/`).
- New entry point `canlab_run_all_tests.m` — works interactively and in CI,
  supports tag filtering, JUnit XML output, and a `'Walkthroughs'` parameter
  (`'skip'` default, `'only'`, `'include'`).
- Shared fixtures in `Unit_tests/helpers/`: `canlab_get_sample_fmri_data.m`,
  `canlab_get_sample_thresholded_t.m`, and a `skip_on_environment_error`
  helper for headless / missing-data CI runners.
- ~70 test cases covering:
  - Core object lifecycle: `fmri_data` construction, `load_image_set`,
    `replace_empty` / `remove_empty` round-trip, `apply_mask`, `ttest` →
    `threshold` → `statistic_image.sig`, `region(t)`, atlas constructor.
  - Functional coverage: `cat`/`split`/`get_wh_image`/`mean`, `regress`,
    display (`montage`, `slices`, `orthviews`, `surface`), `resample_space`,
    similarity, ROI extraction, QC, table, `flip`, `enforce_variable_types`.
  - 30 help-example tests mirroring the per-function docs
    (`canlab_test_help_examples.m`) so the documentation cannot silently
    bit-rot.
  - 10 walkthrough wrapper tests, one per `publish_canlab_help_set1` script,
    run nightly.
- Legacy ad-hoc test scripts moved to `Unit_tests/old_to_integrate/` pending
  rewrite.
- `Unit_tests/README.md` written for the new conventions.

### Added — Continuous integration

- `.github/workflows/test.yml` — runs the fast suite on every push/PR
  (Ubuntu, R2024b) with sibling checkouts of `canlab/Neuroimaging_Pattern_Masks`
  and SPM25 (pinned to tag 25.01.02).
- `.github/workflows/tests-walkthroughs.yml` — nightly cron + `workflow_dispatch`,
  runs the walkthrough tier and uploads JUnit results as an artifact.
- README gains tests status badges.

### Added — Documentation

#### New under `docs/`

- **Top-level reference:** `Object_methods.md` (entry point with class hierarchy
  diagram and curated index), `recasting_objects.md` (idioms for converting
  between class types), `toolbox_folders.md` (one-liner per CanlabCore subfolder).
- **Per-class method indexes (12 pages):** `fmri_data_methods.md`,
  `image_vector_methods.md`, `statistic_image_methods.md`, `region_methods.md`,
  `atlas_methods.md`, `fmridisplay_methods.md`, `brainpathway_methods.md`,
  `fmri_timeseries_methods.md`, `canlab_dataset_methods.md`,
  `fmri_glm_design_matrix_methods.md`, `predictive_model_methods.md`,
  `fmri_mask_image_methods.md`.
- **Per-function code maps (35):** editable PowerPoint sources in
  `docs/code_maps_pptx/` and cropped PNG renders in `docs/code_maps_png/`.
- **Per-function help pages (59):** `docs/individual_functions/<name>.md`, each
  with the reformatted help block, embedded code map, Inputs/Outputs tables, a
  runnable Quick example, a sample output PNG, and See-also cross-links.
- **Atlases / patterns / datasets:**
  - `atlases_regions_and_patterns.md` — combined whole-brain, cortical,
    subcortical, thalamus, brainstem, cerebellum, 7T, networks, signatures;
    ~84 DOI hyperlinks.
  - `sample_datasets.md` — built-in `Sample_datasets/` files and the full
    `load_image_set` keyword registry; ~36 DOI hyperlinks.
- **Agent / build tooling:** `docs/_codemap_tools/` (Python DSL for rendering
  code maps), `docs/help_doc_agent_instructions*` (conventions used to produce
  the docs).

#### Docstring formatting

- 102 class methods (`@fmri_data`, `@image_vector`, `@atlas`, `@region`,
  `@statistic_image`, `@fmridisplay`) and 17 high-traffic stand-alones
  (`load_image_set`, `load_atlas`, `canlab_results_fmridisplay`, the
  `canlab_glm_*` family, `canlab_get_underlay_image`, `canlab_list_files`, ...)
  reformatted to the sphinx layout (`:Usage:`, `:Inputs:`, `:Optional Inputs:`,
  `:Outputs:`, `:Examples:`, `:See also:`). Existing content preserved; missing
  sections added; undocumented options surfaced.
- Standard `'plot'` / `'doplot'` / `'verbose'` / `'doverbose'` flags documented
  wherever the function actually parses them.

#### Project-level

- `CLAUDE.md` added — architecture and conventions guide (also useful for human
  contributors).
- `README.md` expanded with resources, an "Object methods reference" pointer,
  status badges, and editorial clarifications.

### Changed — SPM compatibility

- The brittle pattern
  `switch spm('Ver'); case {SPM5, SPM8, SPM12, SPM25} ... otherwise error('Unknown version!')`
  is replaced with `case {SPM2, SPM99} legacy; otherwise modern` across 26
  sites in 21 files (`@fmri_data/fmri_data`, `@atlas/atlas`,
  `@fmri_timeseries/fmri_timeseries`, `@image_vector/read_from_file`, the
  `Data_extraction/extract_*` family, `Index_image_manip_tools/iimg_*`,
  `Image_computation_tools/*`, `Visualization_functions/montage_clusters*`,
  `Statistics_tools/classify_naive_bayes*`,
  `Filename_tools/expand_4d_filenames`, `Filename_tools/scan_get_files`,
  `diagnostics/scnlab_*`, `Cluster_contig_region_tools/clusters2mask`). SPM26-dev
  and future releases now use the modern code path automatically.

### Fixed

- `@image_vector/flip.m` — multi-image input previously iterated only the first
  image's Z-planes and then errored from `rebuild_volinfo_from_dat` on the
  wrong-size vector. Now errors with id `image_vector:flip:multiImage` pointing
  at `mean(obj)` / `get_wh_image(obj, k)`, and calls `replace_empty` at the top
  so it works on either compressed or padded input.
- `fmridisplay_helper_functions/render_blobs.m` — reshape bug that broke
  `canlab_help_2b_basic_image_visualization` and
  `canlab_help_5_regression_walkthrough` at the `region`/`montage` call.
- `canlab_force_directed_graph.m` — small robustness fix.

### Removed

- `CanlabCore/HRF_Est_Toolbox2/` removed in entirety (superseded by
  `HRF_Est_Toolbox4/`).
- `docs_sphinx_old/`, `docs/.gitignore`, `docs/.idea/` removed.
- `@image_vector/apply_atlas.m` removed (was always-warning deprecated; use
  `apply_parcellation`).
- 4 dead operator overloads removed: `@image_vector/{plus, minus, power, horzcat}.m`.
- `@fmri_data/bootstrap_structure_coeff_diff.m` removed (dead).
- 35 `_old` / `_backup` / `tmp` files removed (each verified to have an active
  sibling): `resample_space_old`, `predictive_model_old`, `region2imagevec_old`,
  `region2imagevec2tmp`, `cluster_table_old`, `canlab_glm_subject_levels_old`,
  `canlab_glm_subject_levels_run1subject_old`, `weighted_reg_old2`,
  `weighted_reg_oldglmfit_old`, `igls_old`, `igls_old2`,
  `getVertexColors_old_backup{,2}`, `construct_model_tmpwork`, `tmp_onsets2dx`,
  and HRF_Est_Toolbox `Old_stuff/More_recent_old_stuff/Example_old`.
- 18 git-tracked `.asv` autosave artifacts removed (`External/spider/`,
  `External/umap/`, `Image_computation_tools/`).

### Added — Small additions

- `@fmridisplay/removeblobs.m`, `@fmridisplay/removepoints.m` — explicit
  blob/point removal methods.
- `@atlas/get_regions_at_crosshairs.m` — query the atlas at the current
  orthviews crosshairs.
- `@fmri_data/saveplots.m` — convenience save method.
- `documentation_template.m` updated to reflect the codified docstring +
  `inputParser` conventions.

### Known follow-ups

Flagged during this work, intentionally not fixed here:

- `@brainpathway/reorder_regions_by_node_cluster.m` appears to duplicate
  `reorder_regions.m`.
- `@fmri_glm_design_matrix/import_onsets.m` references an undefined
  `DesginmatrixTable` and is broken.
- `@predictive_model/`: `crossval` is a literal empty function; `train`, `test`,
  `report`, `plot`, `montage`, `select_features`, `bootstrap`, `error_analysis`,
  `permutation_test`, `confusionchart`, `rocplot` error when called. Documented
  as forthcoming stubs.
- `@brainpathway_multisubject/update_region_connectivity` is a listener stub
  marked "DOES NOT WORK".
- `@statistic_image/check_properties.m` marked "under construction, do not use
  yet" by its own header.
- 15 documentation TODOs remain (3 in atlases — Carmack 2004 midbrain, Brooks
  RVM, inferior olive; 12 in datasets — mostly single-trial dataset citations
  and Ke 2024 topic-label DOI).
- 4 moderate Dependabot vulnerabilities on vendored JS in old docs.

## [v.2.0.0] - 2026-04-29

Previous stable release. See `git log v1.0.0..v.2.0.0` for the full history
between the 2014 import and this tag.

## [v1.0.0] - 2014-08-12

Initial import of the SCN Core Support.
