# Visualization Walkthrough — Status: COMPLETE

**Workspace:** `/Users/f003vz1/Documents/GitHub/CanlabCore`
**Branch:** `fmridisplay-controller-overhaul` (remote `origin` = github.com/canlab/CanlabCore)
**MATLAB:** R2026a via MATLAB MCP server (desktop visible to user; figures appear in MATLAB UI).

The walkthrough and its supporting color-pipeline work are finished. Nothing is
pending. History below is kept for context.

---

## What is DONE

### 1. Central colormap unification (committed & pushed)
- Commit `700f5a54` — `render_blobs.m` + `canlab_colormap.m`: all montage modes
  (split / single / solid / continuous / indexed) build a `central_cm` and color
  through `central_map_slice`.
- Commit `4d3ad933` — `canlab_colormap.m` (NaN/Inf-robust `.map()`, `indexmap`
  case), `@fmridisplay/controller.m` (indexed-atlas legend), `render_layer_surfaces.m`
  (all modes incl. indexed through the truecolor path, `'interp','nearest'`).
- Verified: `max|dRGB| = 0` across montage modes; surface indexmap vertex colors
  equal `cmap(region_index)` exactly. **Do not redo.**

### 2. Walkthrough documentation (committed)
Six pages under `docs/visualization_walkthrough/` (`index.md`, `01`…`06`), figures in
`figures/`, generators in `_gen/`. First committed in `91e9d110`.

### 3. The 8 approved revisions — ALL DONE
1. Surface figures rendered without colorbar legends (`'nobars'` in `wt_save`).
2. `canlab_orthviews` shown in §1 (new `figures/01_canlab_orthviews.png`), plus a
   pointer to `canlab_niivue` linking `docs/canlab_niivue_guide.md`.
3. §1 "How the display methods fit together" — montages/surfaces/OO API; methods
   return an `fmridisplay` object you keep manipulating.
4. §3.4 `addbrain` now points to `help addbrain` / no-arg keyword list and notes
   the `Neuroimaging_Pattern_Masks` / `canlab_load_ROI` dependency.
5. Full catalog of surfaces & composites in new §3.6; montage-sets + surface
   layouts catalog in §5.4.
6. §2.6 custom montage via `fmridisplay.montage()` (orientation / slice_range /
   spacing), new `figures/02_custom_montage.png`.
7. Master script `_gen/visualization_walkthrough.m` — one section per page
   (1.1…6.5), regenerates every figure; referenced near the top of §1 and in
   `index.md`. Ran end-to-end in R2026a with no errors.
8. §6 atlas intro — `help atlas` / Object_methods link, atlas-type table, and the
   `Neuroimaging_Pattern_Masks` repository dependency spelled out.

`_gen/gen_revise.m` (the temporary regen helper) was run, verified, and deleted.
The per-page `gen_0*.m` generators are kept for single-page regeneration.

---

## Gotchas learned
- Reload edited classdefs: `close all force; clear all; clear classes` (plain
  `clear classes` fails while instances exist).
- Capture the RIGHT figure for multi-figure methods:
  `ancestor(o.montage{1}.axis_handles(1),'figure')`, not `gcf`.
- R2026a rejects `func(...){:}` — assign the cell to a temp first
  (`cols = scn_standard_colors(n); cmap = cat(1, cols{:});`).
- `'foursurfaces'` is NOT a montage type. For atlas-on-surface:
  `sh = addbrain('foursurfaces')` then
  `render_on_surface(ctx, sh, 'colormap', cmap, 'indexmap', 'interp', 'nearest', 'nolegend')`
  where `cmap = cat(1, cols{:})`. `'indexmap'` is a flag; the color matrix goes
  through `'colormap'`.
- Dropped from figures (render badly / error), described in prose instead:
  `compact` montage (tiny), `multirow` (region() errors on arrays), `plot()` QC
  (empty axes).
- `docs/tmp_test_fmridisplay.m` is untracked scratch — intentionally NOT committed.
