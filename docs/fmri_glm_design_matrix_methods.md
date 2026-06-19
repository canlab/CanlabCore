# `fmri_glm_design_matrix` methods, organized by area

`fmri_glm_design_matrix` is a container for fMRI first-level design
matrices. Its fields mirror SPM's `SPM.xY`, `SPM.nscan`, `SPM.xBF`,
`SPM.Sess`, and `SPM.xX`, so it interoperates with SPM-style workflows
while adding CANlab-specific affordances such as per-condition basis
sets, parametric modulators added by name, and convenient `build` /
`plot` methods. Once `Sess` (onsets, names, parametric modulators) and
`xBF` (basis set) are filled in, `build` constructs the design-matrix
columns; `plot` shows the basis sets, predictors, and the resulting `X`.

Construct with a TR (`fmri_glm_design_matrix(2)`) and optionally pass
`fieldname, value` pairs (e.g. `'nscan'`, `'units'`, `'onsets'`,
`'condition_names'`, `'pm'`, `'pmnames'`) to populate the model. Type
`methods(my_obj)` in MATLAB for the live list on any instance.

## Properties

| Property | Description |
|---|---|
| `TR` | Repetition time in seconds |
| `build_method` | Build strategy, e.g. `'Separate sessions'` |
| `history` | Provenance log of operations applied |
| `xY` | SPM-style data structure (with `RT` field set from TR) |
| `nscan` | `[1 x s]` number of scans per session |
| `xBF` | Basis-function struct (`name`, `length`, `order`, `bf`, `dt`, `T`, `T0`, `UNITS`, `Volterra`) |
| `Sess` | `[1 x s]` session-array with `U` (onsets), `C` (user covariates), `row`, `col`, `Fc` |
| `xX` | Design-matrix struct with `X`, `iH`, `iC`, `iB`, `iG`, `name` |

## Basic operations

| Method | From | One-liner |
|---|---|---|
| `fmri_glm_design_matrix` | `@fmri_glm_design_matrix` | Constructor; fills `xBF` with default spline basis from TR |
| `add` | `@fmri_glm_design_matrix` | Add fields by name (onsets, condition_names, pm, etc.) |
| `Add_Event_Info` | `@fmri_glm_design_matrix` | Append event onsets / durations / names from an Excel sheet |
| `import_onsets` | `@fmri_glm_design_matrix` | Import onsets from a CSV/Excel design file |
| `replace_basis_set` | `@fmri_glm_design_matrix` | Replace the basis set for a single condition |
| `rotate_to_pca` | `@fmri_glm_design_matrix` | Rotate design columns within conditions to PC projection |

## Display and visualization

| Method | From | One-liner |
|---|---|---|
| `plot` | `@fmri_glm_design_matrix` | Plot basis sets, predictors, and the assembled design matrix |
| `saveplots` | `@fmri_glm_design_matrix` | Save the standard set of design-matrix plots to disk |

## Statistics

| Method | From | One-liner |
|---|---|---|
| `robustfit` | `@fmri_glm_design_matrix` | Robust fit of the design to an `fmri_data` object (with smoothing of weights) |
| `single_trial_estimates` | `@fmri_glm_design_matrix` | Write per-trial beta images using voxel-specific HRFs from a fitted model |

## Misc utilities

| Method | From | One-liner |
|---|---|---|
| `build` | `@fmri_glm_design_matrix` | Assemble `xX.X` from sessions, onsets, and basis sets |
| `build_single_trial` | `@fmri_glm_design_matrix` | Build a single-trial design matrix using one custom HRF per condition |
| `get_session_X` | `@fmri_glm_design_matrix` | Return predictors / delta / covariates / block / names for one session |
| `get_condition_assignments` | `@fmri_glm_design_matrix` | Indicator matrix mapping `X` columns to conditions and parametric modulators |
