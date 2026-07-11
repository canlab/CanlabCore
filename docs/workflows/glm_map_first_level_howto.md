# First-level fMRI time-series modeling with `glm_map` — how-to

A practical, copy-pasteable guide to building, screening, fitting, and reporting
a **single-subject (first-level)** event-related fMRI model with the
[`glm_map`](../glm_map_methods.md) object: importing onsets several ways,
choosing a basis set, simulating data, adding contrasts, diagnostics, fitting
(with AR error models), and visualizing results. This is the **code
walkthrough**; for the conceptual overview see the
[**first-level roadmap**](glm_map_first_level_roadmap.md).

The example is a synthetic **3-run** design with **4 event types**
(`Cue`, `Anticipation`, `Pain`, `Rating`) and **6 motion-like nuisance
regressors per run** — no scanner required.

## Quick reference

| Goal | Call |
|---|---|
| Import an FSL-style event table | `import_onsets(g, 'events.csv', 'TR', tr, 'nscan', n)` |
| Import SPM-style cell arrays | `import_onsets(g, onsets, durations, pmods, 'names', {...})` |
| Import a full SPM model | `import_SPM(g, 'SPM.mat')` |
| Build X from onsets | `build_design(g)` |
| Change the basis set | `replace_basis_set(g, cond, xBF)` |
| Add a contrast | `add_contrasts(g, C, names)` |
| Screen the design | `run_diagnostics(g)` |
| Fit (AR for time series) | `fit(g, data, 'AR', 4)` |
| Threshold / visualize | `threshold(g, ...)`, `montage(g)` |

---

## Setup

```matlab
addpath(genpath('/path/to/CanlabCore/CanlabCore'))
assert(~isempty(which('spm_vol')), 'SPM is not on the MATLAB path.')
TR = 1; nscan = [220 220 220];                 % 3 runs, 220 TRs each
evnames = {'Cue','Anticipation','Pain','Rating'};
```

---

## Section A — Get event information in (three routes)

All three routes flag the entered events as **of interest**; other regressors
(motion / `multiple_regressors`) become **nuisance**.

**1. FSL-style event table** — a `.csv`/`.xlsx`/table with `onset`,
`duration`, and a condition column (a name, or an integer event-type code):

```matlab
T = table([10;25;40]', ...);   % onset / duration / name columns
g = import_onsets(glm_map, T, 'TR', 1, 'nscan', 220, 'units', 'secs');
% (or a filename: import_onsets(glm_map, 'events.csv', 'TR', 1, 'nscan', 220))
```

**2. SPM-style cell arrays** — onsets, durations, and parametric modulators,
one cell per condition:

```matlab
onsets = {[10 40 70]', [25 55 85]'};     % one cell per condition
durs   = {2, 8};                         % per-condition durations
pmods  = {[1 2 3]', []};                 % optional parametric modulator on cond 1
g = import_onsets(glm_map, onsets, durs, pmods, 'TR', 1, 'nscan', 220, ...
                  'names', {'Cue','Pain'}, 'pm_names', {'intensity',''});
```

**3. A complete SPM first-level model** — copies the built design and uses
SPM's event semantics to set the interest/nuisance partition:

```matlab
g = import_SPM(glm_map, '/path/to/SPM.mat');   % struct, SPM.mat, or its folder
```

## Section B — Build the multi-run design

For the walkthrough we build the design from onsets directly, then add motion
covariates per run. Onsets/durations/condition names are laid out
**session-major** (run 1's conditions, then run 2's, …).

```matlab
rng(3);
mkons  = @() sort(randsample(8:6:200, 14))';            % ~14 events / condition
onsets = arrayfun(@(k) mkons(), 1:numel(nscan)*numel(evnames), 'unif', 0);
durs   = {2, 2, 8, 2};                                  % Pain is an 8-s epoch

d = fmri_glm_design_matrix(TR, 'nscan', nscan, 'units', 'secs', ...
        'onsets', onsets, 'condition_names', evnames);
d = add(d, 'durations', repmat(durs, 1, numel(nscan))); % share durations across runs

% 6 motion-like nuisance covariates per run (become nuisance regressors)
for s = 1:numel(nscan)
    d.Sess(s).C.C    = zscore(cumsum(randn(nscan(s), 6)));
    d.Sess(s).C.name = arrayfun(@(i) sprintf('mot%d', i), 1:6, 'unif', 0);
end

g = glm_map(d);
g.is_timeseries = true;          % enable AR error models + HP-filter diagnostics
g = build_design(g);            % onsets ⊛ canonical HRF -> X
summary(g)                       % 660 obs x 33 reg (12 interest, 18 nuisance, 3 intercept)
plot_design(g)
```

![Example first-level design](glm_map_1stlevel_design.png)

The left panels show each event type's canonical-HRF-convolved regressor with
event/duration ticks (the `Pain` epochs are wider); the right heat map shows the
full, block-structured 3-run design.

## Section C — Choosing a basis set

The default is the canonical HRF (one column per event type). To capture
response-shape variability, swap in a multi-basis set for a condition — e.g. a
spline basis for `Pain` — which the design rebuilds, giving several columns
(and several lines in that event's `plot_design` panel):

```matlab
[xBF_hires, ~] = fmri_spline_basis(TR, 'length', 20, 'nbasis', 4, 'order', 3);
g_spline = replace_basis_set(g, 3, xBF_hires);   % condition 3 = Pain -> 4 BFs
g_spline = build_design(g_spline);
```

## Section D — Contrasts and diagnostics

Add a contrast over the events, then screen the design. For a time-series model
`run_diagnostics` reports VIFs/cVIFs (with and without nuisance), efficiency,
and a recommended high-pass filter.

```matlab
% Pain vs Anticipation, over the run-1 event columns (0 elsewhere)
nm = g.regressor_names;
C  = zeros(1, g.num_regressors);
C(find(contains(nm,'Pain'),1)) = 1;  C(find(contains(nm,'Anticipation'),1)) = -1;
g  = add_contrasts(g, C, {'Pain_vs_Antic'});

g = run_diagnostics(g);        % prints the report; stored in g.diagnostics
plot_design(g)                 % adds a VIF figure
g.diagnostics.hpfilter         % recommended high-pass filter (time-series only)
```

![Interest-vs-nuisance VIFs and contrast cVIF](glm_map_1stlevel_vif.png)

The VIF figure shows the events of interest in the **full design** (left, with
motion) and **without nuisance** (middle) — here motion adds only a little
inflation — and the contrast cVIF (right), all with severity reference lines at
VIF = 1, 2, 4, 8.

## Section E — Simulate data and fit

To exercise fitting without real scans, build a synthetic BOLD data set
(`X · betas + noise`) in a template brain space; on real data you would load the
run's time series as an `fmri_data` object instead.
[`fmri_data.sim_data`](../fmri_data_methods.md) builds it in one line from the
convolved design `g.X` and a `[regressors × voxels]` coefficient map. Fit with
an AR(4) error model — recommended for serially correlated BOLD (AR(4) captures
the autocorrelation better than AR(1)).

```matlab
template = load_image_set('emotionreg', 'noverbose');          % borrow a brain space
template = apply_mask(template, fmri_data(which('gray_matter_mask.nii')));
nvox = size(template.dat, 1);

rng(7);
betas_true = zeros(g.num_regressors, nvox);                    % [regressors x voxels]
active = false(nvox,1);  active(randsample(nvox, round(.08*nvox))) = true;
betas_true(contains(g.regressor_names,'Pain'), active) = 4;     % a "Pain" response

sim = sim_data(template, 'design', g.X, 'betas', betas_true, 'noise_sigma', 3);
sim.images_per_session = g.design.nscan;                        % so run intercepts line up

g = fit(g, sim, 'AR', 4);       % autoregressive error model (use 'robust' to down-weight)
```

> **AR(p) dependencies.** AR fitting requires the **Econometrics Toolbox**
> (`fgls`) and the **Signal Processing Toolbox** (`aryule`). `fit` checks for
> both before doing any work and errors with a clear message if either is
> missing; fall back to OLS (`fit(g, sim)`) or `'robust'` in that case.
> Estimation loops over voxels, so whole-brain AR fits are slow.

## Section F — Threshold and visualize

```matlab
g = threshold(g, .001, 'unc', 'k', 5);
t_pain = get_wh_image(g.t, find(contains(g.regressor_names,'Pain'), 1));  % run-1 Pain
canlab_results_fmridisplay(t_pain, 'compact2');
% table(g, 't');   % atlas-labeled table

% Interactive inspection:
canlab_orthviews(t_pain);   % MATLAB 3-plane viewer; names the atlas region at the crosshair
canlab_niivue(t_pain);      % portable web viewer (NiiVue) — email/embed in an HTML report
```

![Thresholded result montage](glm_map_1stlevel_montage.png)

---

## Notes

- **Multi-run intercepts and contrasts.** With several runs the design carries
  one intercept column per run. When fitting a *contrast* on a multi-run design,
  define the contrast over the design's own columns and inspect
  `g.regressor_names` to place the weights; the diagnostics (cVIF, efficiency)
  operate on the design directly.
- **Real data.** Replace the simulated `sim` with `fmri_data('run1*.nii')` (or
  your preprocessed time series); set `sim.images_per_session` to the per-run
  scan counts so run intercepts line up.

## See also

- [`glm_map` methods](../glm_map_methods.md) · [first-level roadmap](glm_map_first_level_roadmap.md)
- [Second-level (group) how-to](glm_map_second_level_howto.md)
- [`fmri_glm_design_matrix`](../fmri_glm_design_matrix_methods.md)
