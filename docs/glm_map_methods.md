# `glm_map` methods, organized by area

`glm_map` is a scikit-learn-style object for **mass-univariate GLM / multiple
regression** on brain images. A single object bundles three things:

1. the **design specification** — either an event/first-level design (onsets,
   durations, parametric modulators, basis set, wrapped in an
   [`fmri_glm_design_matrix`](fmri_glm_design_matrix_methods.md) in `.design`),
   or a pre-built design matrix `.X` supplied directly;
2. the **fitted result maps** — `betas`, `t`, `contrast_estimates`,
   `contrast_t` as [`statistic_image`](statistic_image_methods.md) objects, and
   `df`, `sigma`, `residuals` as [`fmri_data`](fmri_data_methods.md) objects;
3. the **design diagnostics** — variance inflation factors (VIF), contrast VIFs
   (cVIF), leverage, Cook's distance, condition number, design efficiency, and
   (for time-series designs) a recommended high-pass filter, all in the nested
   `.diagnostics` struct.

It is the canonical output type of [`fmri_data.regress`](fmri_data_methods.md),
which now **returns a `glm_map`** (not a plain struct). Type `methods(my_obj)`
in MATLAB for the live list on any instance.

## Two kinds of model: first-level (time series) and second-level (group)

The same object class covers both ends of the standard fMRI analysis pipeline.

- **First-level (within-run time series).** `glm_map` can *store and build*
  single-subject time-series models. Give it event **onsets** (and optional
  durations, parametric modulators, and a basis set) and it convolves them with
  a hemodynamic response to build the design matrix `X` by HRF convolution
  (`build_design`). Mark `g.is_timeseries = true` so that `fit` can use
  autoregressive (`'AR'`) error models, which are appropriate for serially
  correlated BOLD time series. You can assemble the design from onsets
  (`import_onsets`, several formats) or import a complete SPM first-level model
  (`import_SPM`).

- **Second-level (group analysis).** `glm_map` also works with a
  **static / pre-specified design matrix** — the usual situation in a group
  analysis, where each "observation" is one subject's contrast image and `X`
  is a small design (e.g. an intercept plus a covariate). Here you supply `X`
  directly (`glm_map('X', X, 'level', 2)`); there are no onsets to convolve.

A `glm_map` knows which mode it is in via `.level` (1 = first-level,
2 = second-level) and `.is_timeseries`.

## Ways to create one

You never build a `glm_map` field by field. There are a few entry points, and
they all return the **same** object type.

1. **From `fmri_data.regress` (quick path).** Put a design matrix in `dat.X`
   and call `regress`; the result is a `glm_map`:

   ```matlab
   dat   = load_image_set('emotionreg');          % 30 contrast images
   dat.X = [zscore((1:30)') ones(30,1)];          % covariate + intercept
   g     = regress(dat, .001, 'unc');             % g is a glm_map
   montage(g.t);                                  % thresholded t maps
   ```

   For backward compatibility the historical struct field names still work as
   aliases: `g.b` → `g.betas`, `g.con_t` → `g.contrast_t`,
   `g.contrast_images` → `g.contrast_estimates`, `g.resid` → `g.residuals`,
   `g.variable_names` → `g.regressor_names`, `g.C` → `g.contrasts`.

2. **As an estimator (the scikit-learn-style API).** Specify the design first,
   screen it, add contrasts, then `fit`:

   ```matlab
   g = glm_map('X', X, 'level', 2, 'regressor_names', {'cov','intercept'});
   g = add_contrasts(g, [1 0], {'cov_effect'});
   g = run_diagnostics(g);                        % VIF / cVIF / efficiency report
   g = fit(g, dat);                               % runs fmri_data.regress
   table(g, 'contrast'); montage(g, 'contrast_t');
   ```

3. **From event onsets (first-level).** Wrap an `fmri_glm_design_matrix`, or
   import onsets directly into the object:

   ```matlab
   d = fmri_glm_design_matrix(TR, 'nscan', nscan, 'units', 'secs', ...
           'onsets', onsets, 'condition_names', {'A','B'});
   g = glm_map(d);  g.is_timeseries = true;  g = build_design(g);

   % or, in one call (FSL event table, SPM-style cells, or a file):
   g = import_onsets(glm_map, 'events.csv', 'TR', 2, 'nscan', 200);
   ```

4. **From an SPM first-level model.** `import_SPM` copies an SPM12/SPM25
   `SPM.mat` design into the object and flags the event regressors as of
   interest:

   ```matlab
   g = import_SPM(glm_map, '/path/to/SPM.mat');
   ```

5. **Re-cast a regress-style struct.** `glm_map(out_struct)` maps the fields of
   a results structure onto the object's properties.

## Design

- **Value semantics.** Every mutating method returns a NEW object — write
  `g = run_diagnostics(g)`, not `run_diagnostics(g)`.
- **Property names mirror `fmri_data.regress`.** Related outputs are grouped
  into nested structs — `.input_parameters` (fit options),
  `.input_image_metadata` (provenance of the fitted images), and `.diagnostics`
  — so a `glm_map` looks just like the `out` struct `regress` used to return.
- **Regressor roles.** `.wh_interest`, `.wh_nuisance`, and `.wh_intercept` are
  logical indicators over the columns of `X`. For event designs they come from
  the design partition (event regressors are of interest, covariates such as
  motion parameters are nuisance); for direct designs, non-intercept columns
  are of interest by default and you mark covariates of no interest with
  `.nuisance_columns`. Diagnostics use these to report VIFs **with and without**
  nuisance covariates, and contrast/efficiency tools operate on the regressors
  of interest only.
- **The design matrix is exposed at the top level.** `g.X` reads through to the
  built design (`g.design.xX.X` in event mode); `disp(g)` and `summary(g)` print
  the number of observations and regressors, broken down into of-interest,
  nuisance, and intercept.
- **Intercept and contrast handling at `fit`.** `fit` ensures the model has an
  intercept before estimating: a single overall intercept is enforced as the
  **last** column (an existing one is moved, never duplicated; one is added if
  absent), while multi-run designs that already carry per-run baseline columns
  are left as they are. The model is fit first (all betas estimated) and
  contrasts are applied afterward to form the contrast maps. A contrast vector
  that is one (or more) elements **short** is padded with 0 weights for the
  intercept, so you can write contrasts over just the regressors of interest.
  When you supply `regressor_names` for a design that includes an explicit
  intercept column, **include a name for the intercept** so names stay aligned
  with the columns (a missing one defaults to `'Intercept'`).
- **A note on naming.** `diagnostics` is a *property* (the results struct); the
  method that computes it is the verb **`run_diagnostics`** (call it as
  `g = run_diagnostics(g, ...)`; results land in `g.diagnostics`).

## Methods, by area

### Construct / import a design

| Method | What it does |
|---|---|
| `glm_map(...)` | Constructor: wrap an `fmri_glm_design_matrix`, take `'X', X` directly, re-cast a regress-style struct, or set any property. |
| `import_onsets(g, source, ...)` | Build the wrapped design from event onsets — an **FSL-style** event table (`.csv`/`.xlsx`/table with onset/duration/name or integer event-code columns), **SPM-style** cell arrays (onsets/durations/parametric modulators, one cell per condition), or a file. Bootstraps a design from `'TR'`/`'nscan'` and builds it when enough info is present. |
| `import_SPM(g, SPM)` | Import an SPM12/SPM25 first-level model (`SPM` struct, `SPM.mat` path, or directory); optionally load betas. Derives the of-interest / nuisance partition from SPM's event semantics. |
| `build_design(g)` | Convolve onsets/durations with the basis set to build `X` (event/first-level mode). |
| `replace_basis_set(g, cond, xBF)` | Swap the basis set for one condition (e.g. canonical HRF → spline or HRF+derivatives), rebuild, clear stale results, and re-fit if data are supplied. |

### Contrasts

| Method | What it does |
|---|---|
| `add_contrasts(g, C, names)` | Append one or more linear contrasts (one row per contrast, over the regressors). |
| `create_orthogonal_contrast_set(g)` | Assign an orthogonal, zero-sum contrast set spanning the regressors of interest (placeholder names). |

### Diagnostics

| Method | What it does |
|---|---|
| `run_diagnostics(g, ...)` | Compute and (by default) print a narrative report: VIFs/cVIFs for the **full design** and for the **regressors of interest only**, leverage, Cook's distance, scaled condition number, design **efficiency** (`calcEfficiency`), and — for time-series designs — cumulative power by frequency with a recommended high-pass-filter cutoff. Stored in `g.diagnostics`. |
| `plot_design(g)` | Plot the design: one panel per event type showing the actual basis-convolved regressor(s) with event/duration boxes, plus a heat map of the full design matrix. When diagnostics are available, also draws VIF/cVIF plots with severity reference lines. |

### Fit and report results

| Method | What it does |
|---|---|
| `fit(g, data, ...)` | Run the regression on an `fmri_data` object (the compute engine is `fmri_data.regress`); populate `betas`, `t`, contrasts, `df`, `sigma`, and diagnostics. Options: `'robust'`, `'AR', order` (time series), `'residuals'`, `'pthresh'`, `'thresh_type'`. |
| `threshold(g, p, type, ...)` | Re-threshold the stored `t` / `contrast_t` maps without refitting (delegates to `statistic_image.threshold`). |
| `table(g, which_map, ...)` | Atlas-labeled results table for a chosen map (`'t'`, `'betas'`, `'contrast'`, `'contrast_t'`). |
| `montage(g, which_map, ...)` | Brain montage of a chosen result map. |
| `summary(g)` | Narrative summary: analysis name, model (level + input variables), design diagnostics, and (once fitted) how many maps, the threshold, and the number of significant voxels per regressor and contrast. |

### Housekeeping

| Method | What it does |
|---|---|
| `disp(g)` | One-screen object listing (also runs when you type the variable name). |
| `validate_object(g)` | Ensure every nested struct (`input_parameters`, `input_image_metadata`, `diagnostics`, `fit_parameters`) exposes its full set of fields. |
| `check_properties(g)` | Coerce field types and check contrast bookkeeping. |

## Key properties

| Property | Holds |
|---|---|
| `design` | the wrapped `fmri_glm_design_matrix` (event/first-level mode) |
| `X` *(dependent)* | the `[observations × regressors]` design matrix |
| `level` / `is_timeseries` | 1 = first-level, 2 = second-level; AR errors allowed when `is_timeseries` |
| `regressor_names`, `contrasts`, `contrast_names` | column names and the `[regressors × contrasts]` contrast matrix |
| `wh_interest` / `wh_nuisance` / `wh_intercept` | logical role indicators over the columns of `X` |
| `betas`, `t`, `contrast_estimates`, `contrast_t` | result maps (`statistic_image`) |
| `df`, `sigma`, `residuals` | per-voxel error df, residual SD, residuals (`fmri_data`) |
| `diagnostics` | nested struct of VIF/cVIF, leverage, Cook's D, condition number, efficiency, HP filter |
| `input_parameters`, `input_image_metadata`, `fit_parameters` | options used and input-image provenance |
| `analysis_name`, `notes`, `history`, `warnings` | provenance / metadata |

## Workflows

Two end-to-end walkthroughs build a design, screen it, fit it, and report
results with `glm_map`:

- **[First-level fMRI time-series modeling](workflows/glm_map_first_level_roadmap.md)** — an event-related design with several event types, nuisance covariates, and multiple runs; importing onsets, choosing a basis set, simulating data, contrasts, diagnostics, fitting (with AR models), and visualizing results. ([how-to](workflows/glm_map_first_level_howto.md))
- **[Second-level fMRI group analysis](workflows/glm_map_second_level_roadmap.md)** — group regression on contrast images: creating the object two ways, outlier handling, gray/white/CSF covariates and `normalize_by_wm_csf`, OLS vs robust fitting, diagnostics, and visualizing results. ([how-to](workflows/glm_map_second_level_howto.md))

## See also

- [`fmri_data`](fmri_data_methods.md) — the data object and `regress` (the compute engine)
- [`fmri_glm_design_matrix`](fmri_glm_design_matrix_methods.md) — the wrapped first-level design object
- [`statistic_image`](statistic_image_methods.md) — the result-map type (`threshold`, `montage`, `table`)
