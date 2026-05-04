# `canlab_dataset` methods, organized by area

`canlab_dataset` is a generic subject-by-variable container for
behavioral data and fMRI/study meta-data. It is not a subclass of
`image_vector`; it organizes data at five potential levels of analysis:
Experiment, Subject, Event (e.g. trial-level), Sub-Event, and
Continuous. Each level stores a `names`, `type`, `units`, `descrip`,
`data`, and `textdata` field, with cell arrays of one cell per subject
for the trial- and continuous-level data. The most-used levels are
`Subj_Level` (one row per subject) and `Event_Level` (one cell per
subject, each with a trials-by-variables matrix), which together support
multi-level mixed-effects analyses.

Constructing with the `'fmri'` keyword pre-populates the standard
demographic and fMRI event-level variable names (SessionNumber,
EventOnsetTime, etc.) for fMRI design-matrix generation. Type
`methods(canlab_dataset)` in MATLAB for the live list.

## Properties

| Property | Description |
|---|---|
| `Description` | Struct with `Experiment_Name`, `Missing_Values`, and per-level descriptive text |
| `Subj_Level` | Struct with `id`, `names`, `type`, `units`, `descrip`, `data`, `textdata` for subject-level variables |
| `Event_Level` | Struct (same fields, no `id`) for event/trial-level variables; `data` is a cell-per-subject |
| `Sub_Event_Level` | Same shape as `Event_Level`, for variable-size sub-event observations per trial |
| `Continuous` | Same shape as `Event_Level`, for time-series samples per subject |
| `wh_keep` | Struct with logical inclusion masks for filtering subjects/events |

## Basic operations

| Method | From | One-liner |
|---|---|---|
| `canlab_dataset` | `@canlab_dataset` | Constructor; `'fmri'` flag adds standard demographic + event variable names |
| `add_vars` | `@canlab_dataset` | Add or replace one or more variables at a given level |
| `concatenate` | `@canlab_dataset` | Flatten Subject- and Event-level data across all subjects |
| `replace_values` | `@canlab_dataset` | Replace values matching a predicate with a fill value |
| `select_trials_and_subjects` | `@canlab_dataset` | Filter events by string-pattern match on a text variable; optionally extract another |
| `spm2canlab_dataset` | `@canlab_dataset` | Pull Event_Level onsets/durations from a subject's `SPM.mat` |

## Display and visualization

| Method | From | One-liner |
|---|---|---|
| `bars` | `@canlab_dataset` | Bar plot (wraps `barplot_columns`) for selected variables |
| `histogram` | `@canlab_dataset` | Histogram of one variable, subject- or event-level |
| `plot_var` | `@canlab_dataset` | Mean and SE of a variable across events per subject |
| `scatterplot` | `@canlab_dataset` | Scatter of two variables, single- or multi-line per subject |
| `scattermatrix` | `@canlab_dataset` | Scatterplot matrix of pairwise event-level variables |

## Statistics

| Method | From | One-liner |
|---|---|---|
| `glm` | `@canlab_dataset` | Predict Y from X with a single-level GLM |
| `glm_multilevel` | `@canlab_dataset` | Multilevel (subject-as-random-effect) GLM |
| `mediation` | `@canlab_dataset` | Single- or multilevel mediation (wraps `mediation` toolbox) |
| `ttest2` | `@canlab_dataset` | Two-sample t-test on a subject-level variable across two `wh_keep` groups |
| `correlation` | `@canlab_dataset` | Pairwise correlations between two variables (raw and within-subject-centered) |
| `reliability` | `@canlab_dataset` | Odd/even Spearman-Brown reliability of a within-subject variable |

## Tables and descriptives

| Method | From | One-liner |
|---|---|---|
| `print_summary` | `@canlab_dataset` | Print summaries for every (or specified) variable |
| `list_variables` | `@canlab_dataset` | List variable names and descriptions |
| `get_descriptives` | `@canlab_dataset` | Return descriptive stats per variable, per level |

## Data extraction

| Method | From | One-liner |
|---|---|---|
| `get_var` | `@canlab_dataset` | Extract one or more variables in matrix and per-subject cell formats; supports conditional filters |

## Misc utilities (I/O)

| Method | From | One-liner |
|---|---|---|
| `read_from_excel` | `@canlab_dataset` | Read from an Excel design + per-subject files into the dataset |
| `write_text` | `@canlab_dataset` | Flatten and write Subject- and Event-level data to text files |
