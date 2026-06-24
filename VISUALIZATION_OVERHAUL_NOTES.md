# Visualization Overhaul — Design Notes

Status: **Phase 1–2 + mutators + widget IMPLEMENTED** (handle class shipped); render-pipeline unification + engine consolidation still deferred (see §11)
Author: drafted with Claude Code, 2026-06-11; implemented 2026-06-24
Scope: `@fmridisplay` class and the broader CANlab visualization surface (montage / surface / orthviews / isosurface paths).

This is a resumable design document, in the same spirit as `predictive_model_DEVELOPMENT_NOTES.md`. It records what the visualization code looks like today, the core architectural decision, the target design, and a phased, non-disruptive migration path.

---

## 11. Implementation log (2026-06-24)

The load-bearing architectural change is now shipped. Done in `@fmridisplay/`:

- **`fmridisplay` is now a handle class** (`classdef fmridisplay < handle`). A full
  call-site audit across CanlabCore and CANlab_help_examples found **no** value→handle
  aliasing hazard (`o3 = o2; ... mutate(o3); ... use(o2)`); every site uses the
  self-reassigning `o2 = method(o2, ...)` idiom, which is safe. The two var-to-var
  copies that exist (`render_cifti_on_brain.m:246`, `@image_vector/montage.m:270`) are
  benign (store-at-end / copy-then-clear).
- **Blob layers retain their source + render options** (`addblobs.m`): each
  `activation_maps{}` entry now also stores `source_region`, `source_object` (the
  original, e.g. a `statistic_image`, enabling re-threshold up *and* down),
  `render_args` (the exact option set passed to `render_blobs`), `wh_montage`,
  `wh_surface`, and `applied_threshold`.
- **New in-place mutator methods** (re-render from retained source):
  `refresh`, `rethreshold`, `set_colormap`, `set_opacity`.
- **Surfaces are now unified views** (§4.2 / §7.4). Surface drawing is factored into a
  shared `render_layer_surfaces` method that derives colors from the layer's stored
  `render_args`, and `addblobs`, `surface`, and `refresh` all route through it. Concretely:
  (a) `surface(o2)` registers on the **same** handle and **pulls in** any blob layers that
  already exist, so a surface added *after* `addblobs` still shows the blobs;
  (b) `addblobs`/`removeblobs` act on montages **and** surfaces together (remove restores
  the gray anatomy via `addbrain('eraseblobs')`);
  (c) `refresh` — and therefore `rethreshold`/`set_colormap`/`set_opacity` — re-render
  surfaces too;
  (d) multi-surface keywords `'foursurfaces'` / `'foursurfaces_hcp'` /
  `'foursurfaces_freesurfer'` add the four canonical views (L/R lateral + L/R medial) as
  separate registered views on the same object instead of being silently ignored.
- **Controller widget**: `controller(obj)` opens a `uifigure` bound to the instance
  (back-pointer in figure appdata under `'fmridisplay_obj'`), with per-layer opacity
  slider, colormap dropdown, p-threshold field, and visibility toggle — the MATLAB-side
  analog of the NiiVue control panel (§6). The panel is rebuilt in place (not duplicated)
  on a repeat call, **auto-refreshes** when `addblobs`/`removeblobs` change the layer set
  (`update_controller`, via the Transient `controller_handle` property), initializes each
  control from the layer's current state, and is **type-aware**: the threshold field
  re-thresholds a `statistic_image` layer by p-value and an `fmri_data`/mask layer by raw
  value. `HandleVisibility='on'`, so `close all` closes it.
- **Robust to closed windows**: `prune_dead_views` drops montage/surface views whose
  figures the user closed (with a short note); called at the top of `addblobs`,
  `removeblobs`, and `refresh`, so closing a window no longer causes
  "Invalid or deleted object" errors. Surface erase/render also filter to live handles.
- **Layer-aware `rethreshold`** (and `montage` source retention): `image_vector.montage`
  now passes the original object via a new `addblobs` `'source_object'` override, so a
  `montage(t)` layer keeps its `statistic_image` and can be re-thresholded by p-value.
  `rethreshold` dispatches by source class: `statistic_image`→p-value, `image_vector`→raw
  value, `region`→magnitude cutoff; honors `'layers'` for per-layer thresholding.
- **Bug fix surfaced by the green-tests bar**: a *bare* `'colormap'` flag (e.g.
  `montage(r, o2, 'colormap')`) was a pre-existing crash (`render_blobs` demands an
  n×3 matrix). `addblobs` now strips a bare `'colormap'` before forwarding. This
  un-broke `canlab_help_5_regression_walkthrough`.
- **Tests**: new `Unit_tests/image_vector/canlab_test_fmridisplay_handle.m` (9 tests:
  handle identity, source retention, refresh, monotone re-threshold, set_colormap,
  set_opacity, controller binding, canonical-call-site back-compat). Verified green:
  full unit tier **157 pass / 0 fail / 2 skip**; **walkthroughs 10/10**.

**Still deferred** (the doc's own "separate cleanups", §7.3; high-risk refactors of
1000+ line files that would jeopardize the green-walkthroughs bar — left for a focused
follow-up now that the layer architecture exists to support them): unifying
`render_blobs` + `render_on_surface` into one color/threshold pipeline; routing **orthview**
refresh through the same layer mechanism (montages and surfaces are now unified; orthviews
are not yet); composited multi-layer surfaces (surfaces currently show the last-rendered
layer, since `render_on_surface` overwrites vertex colors); symmetric pull-in for montages
added after blobs; the controller's visibility toggle hides montage blobs only (not surface
coloring); collapsing the duplicate isosurface/orthviews engines; deleting deprecated
`surface_cutaway`.

---

## 1. Motivation

The visualization layer is a "Frankenstein": several independent engines that each re-implement montage/surface/orthviews and color/threshold logic, most of which produce **static** figures and return raw handles. `fmridisplay` was created to manage multi-view displays (add/remove blobs from a persisting figure), but it captures too little state to support the interaction we actually want — live re-thresholding, colormap changes, and opacity changes after a display already exists.

`canlab_niivue` (the NiiVue web viewer) is the contrasting model: each overlay/underlay is a persistent layer object with mutable `cal_min`, `colormap`, `opacity`, and a property change triggers a re-render. We want that property-centric, stateful behavior on the MATLAB side, across montages **and** surfaces, controllable from a widget — without collapsing to a single global display instance.

---

## 2. What exists today (and why it limits us)

### 2.1 `fmridisplay` is a value class that stores rasterized graphics, not data

`@fmridisplay/fmridisplay.m` is `classdef fmridisplay` — a **value** class. Fields:

- `overlay` — the underlay **filename** (string); `SPACE` — sampling grid from `spm_vol`.
- `activation_maps{}` — one struct per blob layer: `mapdata` (the **already-thresholded, resampled mask**), `cmaprange`, `mincolor`/`maxcolor`/`color`, `blobhandles`.
- `montage{}`, `surface{}`, `orthviews{}` — registered views (axis handles + orientation/direction + patch handles).
- `history`.

Decisive limitation: in `addblobs.m` → `render_blobs.m`, each slice's blob is drawn as a `surf` object whose color and alpha are **baked into CData/AlphaData at render time** (`render_blobs.m:749-781`). The continuous source values and the *threshold that produced the region* are discarded — `activation_maps` keeps only the post-threshold mask. So "re-threshold" or "change colormap" is structurally impossible without deleting and re-adding blobs from a source the caller must still be holding. The color fields that *are* stored (`cmaprange`, `mincolor`, ...) are write-only bookkeeping for the legend; nothing reads them back to re-render.

### 2.2 Two entry points, and most methods bypass both

A user obtains an `fmridisplay` either from the constructor or from `canlab_results_fmridisplay`. Of ~18 class visualization methods, only three return one: `@fmri_data/montage`, `@region/montage`, `@atlas/montage`. Everything else returns **raw handles** and never touches `fmridisplay`:

- `@image_vector/surface.m` → `[all_surf_handles, pcl, ncl]`
- `@image_vector/orthviews.m` / `@statistic_image/orthviews.m` → `region`
- `@image_vector/isosurface.m`, `render_on_surface.m`, `display_slices.m` → handles / nothing

So surfaces and orthviews are entirely outside the managed-display world.

### 2.3 Heavy duplication underneath

- **Montage**: `fmridisplay` path vs legacy `montage_clusters` (732 lines, `spm_orthviews`) vs `cluster_orthviews_montage`.
- **Orthviews**: `cluster_orthviews` (SPM backend) vs `canlab_orthviews` (1890-line SPM-free reimplementation, same API).
- **Surface**: `render_on_surface` (1138 lines), `@region/surface`, `surface_cutaway` (DEPRECATED), `addbrain` (1391 lines), `canlab_canonical_brain_surface_cutaways`, `mask2surface`.
- **Isosurface**: 4 separate implementations (`@image_vector`, `@region`, `@atlas`, `mask2surface`).
- Color/threshold logic is re-implemented in each path.

### 2.4 Why NiiVue feels better

Inverted data model: each layer is a live object with mutable `cal_min`/`colormap`/`opacity`; a property change calls `updateGLVolume()`. State lives in persistent **layer objects**, not in the pixels. That is the pattern to bring to MATLAB.

---

## 3. Core design decision

**Make the managed display a `handle` class that stores source data and display parameters as live layers, and re-renders from them — instead of a value class that stores baked graphics.** Everything else follows from this.

Why `handle` (not value):
- A handle instance can live in the workspace, be referenced by its figure(s), and be referenced by a controller widget — all pointing at **one** object, no copies. A value class cannot, which is the root reason in-place re-thresholding does not exist today.
- In-place mutation (`obj.threshold = 3` → auto-refresh) becomes possible.

---

## 4. Target architecture

### 4.1 Layer abstraction — `canlab_display_layer` (handle class)

One per overlay or underlay. Holds:

- **source**: the `statistic_image` / `fmri_data` / `atlas` / `region` at full resolution, **unthresholded**.
- **display params**: `threshold` (value + type: unc / fdr / k-extent), `colormap` / `colormap_negative`, `cmaprange`, `opacity`, `visible`, `outline`.
- **handle registry**: graphics handles drawn into each view.

`set.threshold`, `set.colormap`, etc. trigger `refresh()`, which re-derives displayed values from source and updates (or deletes+redraws) only that layer's handles in every attached view. Underlays are layers too, so underlay opacity/colormap become adjustable.

### 4.2 View abstraction

Each montage, surface, or orthview is a registered view (axes + geometry + which layers are drawn there). Layers render into all attached views; a view added later pulls in all existing layers. Generalizes today's `obj.montage{}` / `obj.surface{}`, but routes montage/surface/orthviews through **one** color+threshold pipeline instead of three.

### 4.3 Display object — `fmridisplay` v2 (handle class)

Owns `layers{}` and `views{}`, plus methods preserving today's surface:
`addblobs`, `removeblobs`, `montage`, `surface`, `orthviews`, `rethreshold`, `set_colormap`, `set_opacity`, `refresh`.

This realizes "register more properties of both overlay and underlay for multiple views": the layer holds the properties, views are many, one `refresh` propagates a change everywhere.

---

## 5. Instance vs figure — the ownership model

**Recommendation: the object instance is the single source of truth; each figure stores a back-pointer to it in appdata.** This is both candidate options used correctly, not a choice between them.

- Because it is a **handle** class, the instance can be referenced by the workspace variable, its figure(s), and its widget — one object, no copies.
- `setappdata(fig, 'canlab_display', obj)` lets a figure close, a click, or a widget callback route back to the owning instance. The figure holds a **reference**, not a duplicate — no sync problem.
- **Multiplicity is free**: every `fmri_data`/`region` displayed constructs its own handle instance with its own figures and widget. No singleton anywhere — multiple independent displays of different datasets coexist in the workspace by construction.

Do **not** store the actual data in the figure; that is what forces the brittle "figure is the source of truth for some things" split today.

---

## 6. The controller widget

A `uifigure`-based controller (`canlab_display_controller`) **bound to one display instance**. Lists its layers with per-layer: threshold slider, colormap dropdown, opacity slider, visibility/outline toggles, and layer add/remove/reorder. Every control sets a property on the layer and calls `refresh()`. Because it is bound to an instance, opening two displays gives two widgets. This is the MATLAB-side analog of the NiiVue control panel. (Use `uifigure`/`uigridlayout`; the `matlab-build-app` skill covers the patterns.)

---

## 7. Migration path (non-disruptive, phased)

Hard constraint: do not break existing `o2 = canlab_results_fmridisplay(...); addblobs(o2, region(t))` call sites.

1. **Introduce the handle class alongside the value class.** New class (e.g. `fmridisplay2` internally, or `fmridisplay` reimplemented behind a compatibility façade). `addblobs`/`montage`/`removeblobs` keep their exact current signatures and still work value-style (`o2 = addblobs(o2, ...)` returns the same handle). This value-vs-handle return semantics is the riskiest compatibility point — needs a deliberate shim plus a test script under `Unit_tests/`.
2. **Make layers store source data.** `addblobs` retains the source `region`/`statistic_image` on the layer (it already receives it), so re-threshold has something to work from, without changing what callers pass.
3. **Centralize one render pipeline** for montage + surface + orthviews; route `render_blobs` and `render_on_surface` color/threshold logic through it. Delete deprecated `surface_cutaway`; collapse duplicate isosurface/orthviews engines as separate cleanups.
4. **Funnel bypass methods through the object.** `@fmri_data/surface`, `@image_vector/orthviews`, etc. gain an option (then a default) to return a display instance with a registered view, while still returning handles where back-compat needs it.
5. **Add `rethreshold` / `set_colormap` / `set_opacity` / `refresh`**, then the widget.

Each phase is independently shippable and leaves existing scripts working.

---

## 8. Risks

- **Value→handle semantics** is the real hazard: handle objects mutate in place, so old code relying on copy-on-assignment (`o3 = o2; addblobs(o3, ...)` expecting `o2` untouched) changes behavior. Needs an explicit audit + the back-compat shim; most likely thing to bite.
- **Memory**: storing full unthresholded source per layer costs RAM (many layers × whole-brain `fmri_data`). Mitigate by storing a downsampled-to-display-space copy plus a handle to the original, or loading lazily.
- **Re-threshold from source** only works toward *lower* stringency if the unthresholded data is retained — be explicit that the layer stores pre-threshold values.

---

## 9. Suggested first step

Prototype **Phase 1**: a handle-class skeleton with one layer (storing source + display params) and one montage view, plus a `refresh()` that re-renders that layer. This de-risks the value→handle decision before committing the full migration. Pair it with a `Unit_tests/` script that exercises the back-compat `o2 = addblobs(o2, ...)` call pattern against both the old and new classes.

---

## 10. Reference: key files

| Concern | File |
|---|---|
| Display class (value) | `CanlabCore/@fmridisplay/fmridisplay.m` |
| Add/remove blobs | `CanlabCore/@fmridisplay/addblobs.m`, `removeblobs.m` |
| Montage view | `CanlabCore/@fmridisplay/montage.m` |
| Surface view | `CanlabCore/@fmridisplay/surface.m` |
| Blob rasterizer (color/alpha baked here) | `CanlabCore/fmridisplay_helper_functions/render_blobs.m` |
| Surface colorizer | `CanlabCore/@image_vector/render_on_surface.m` |
| Wrapper / 2nd entry point | `CanlabCore/Visualization_functions/canlab_results_fmridisplay.m` |
| Surface loader | `CanlabCore/Visualization_functions/addbrain.m` |
| Orthviews engines | `cluster_orthviews.m`, `canlab_orthviews.m` |
| NiiVue model to emulate | `CanlabCore/Visualization_functions/canlab_niivue/` |
