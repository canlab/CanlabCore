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
  - **Later additions:** controls now **update in place** (not a full rebuild) when a layer's
    state changes from the command line, so `rethreshold`/`set_colormap`/`set_opacity` keep an
    open panel in sync (they call `update_controller`); each GUI action **echoes the
    equivalent code line** to the command window (e.g. `han = set_colormap(han, 'color',
    [1 0 0]);`), using the bound variable name from `inputname` (preserved across rebuilds via
    figure appdata; falls back to `obj`). The colormap dropdown gained more options (split,
    warm, cool, winter) plus a **`solid colour…`** entry that opens the built-in `uisetcolor`
    palette. Cosmetic: the figure is **light green**, the title shows the **bound variable
    name** (to tell multiple controllers apart), and the opacity slider's tick marks are
    removed. (True perceptual maps — inferno/viridis/etc. — are deferred; they need colormap
    matrices threaded through `render_blobs`/`render_on_surface`, riskier than the max/min
    ramps used now.)
- **Controller redesign + more colormaps.** Single-column layout per layer with a colormap
  **title stripe**, an opacity slider, a type-aware **threshold slider**, a Colors dropdown with
  a live **preview swatch**, and a visibility toggle; a footer groups **Re-render / Remove
  legend / Close**. The threshold slider is a **log-scale** p-value slider (ticks at
  .001/.005/.01/.05/.1, so .05–.1 sit close and .001–.005 spread out) for `statistic_image`
  layers, or a linear raw `|x|` slider anchored at 0 and the 99.9th percentile of `|data|` for
  `fmri_data`/region layers. Colormap menu adds `split (mango)` and `seafire` split presets
  alongside warm/cool/winter and the `solid colour…` picker. Window colour is set by the
  `FIG_COLOR` constant at the top of `controller.m` (currently `[1 .5 0]`).
  - **Later tweaks:** **mango** is now the **default** blob colormap (the `addblobs` default
    split colours changed to mango). Footer is **Re-render + Toggle legend** (legend on/off; the
    Close button was dropped). Each layer panel has a **Remove layer** button (new
    `remove_layer(obj, k)` method — removes one layer, vs `removeblobs` which removes all).
    **Clicking the colour swatch** opens the colour picker (a quick way to set/re-pick a solid
    colour; note re-selecting the already-selected `solid colour…` dropdown item can't re-fire
    in MATLAB, so the swatch is the re-pick path). The p-value slider now extends **below .001
    down to ~0** (1e-6 floor, labelled `~0`).
  - `rethreshold` already accepts a cluster-extent `'k'` argument (passed through to
    `threshold`), verified working (e.g. `rethreshold(o2, .01, 'unc', 'k', 50)`).
- **`squeeze_figure(obj)`** removes the top/bottom white space from montage figures (a montage
  often fills only ~40% of its figure's height): it shrinks the figure and remaps the slice
  axes to fill it, scaled so each slice keeps its pixel size, aspect, and horizontal placement
  — only the empty margins go. Opt-in (doesn't change montage defaults); skips figures that
  also contain surfaces.
- **Surface colorbar legends are parented to the surface figure** (`render_on_surface`), not
  to `gcf` — previously a surface legend could land on the montage window. New
  `remove_legend(obj)` method deletes the surface colorbar legends (keeps the blobs).
- **Surface colormaps now match the montage spec.** `render_on_surface` understands
  `color`/`colormap`/`splitcolor`/`pos_colormap`/`neg_colormap` but NOT `maxcolor`/`mincolor`
  (`render_on_surface.m:222`), so warm/cool/winter (max/min ramps) and solid colours were
  silently dropped on surfaces and fell back to the default split hot/cool (e.g. blue
  negatives under "warm", and cool/winter not updating). `render_layer_surfaces` now
  translates the layer's colour spec into `pos_colormap`/`neg_colormap` ramps
  (`surface_colormaps_from_args`) and strips the tokens render_on_surface can't read. Solid
  colours map to the same colour on both signs (no stray default-cool negatives). *Residual:*
  the split hot/cool map still looks **slightly** different between montage and surface
  (different colormap construction in `render_blobs` vs `render_on_surface`); fully matching
  is the deferred render-pipeline unification.
- **Multi-surface keyword uses a dedicated figure.** `surface(o2, 'foursurfaces*')` opened the
  2×2 layout in the *current* figure, so without a preceding `figure;` it drew over the
  montage and appeared to render nothing. It now reuses an empty current figure or opens a new
  one.
- **Public vs internal method surface**: the helper methods `activate_figures`,
  `prune_dead_views`, `refresh`, `render_layer_surfaces`, `update_controller` are declared
  `Hidden` (signature-only block in `fmridisplay.m`), so they no longer clutter
  `methods(obj)`/tab-completion. `Hidden` (not `Access=private`) keeps them callable, which
  matters: `activate_figures` is used by help-example scripts and `refresh` by tests. All
  user-facing methods stay public.
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
added after blobs; collapsing the duplicate isosurface/orthviews engines; deleting deprecated
`surface_cutaway`.

- **Per-layer surface visibility** (deferred; evaluated). The controller's *Visible* toggle
  hides a layer's **montage** blobs (sets `Visible` on `blobhandles`) but not its surface
  coloring, because surfaces don't store per-layer graphics — `render_on_surface` bakes one
  layer's colors into the patch vertices. A clean fix needs a per-layer `visible` flag that is
  **respected by every surface-draw path** (`addblobs`, `refresh`, `surface` pull-in): on a
  visibility change, reset the surfaces to gray and re-render only the visible layers in order
  (top visible wins). It's moderate, not a one-liner, and shares machinery with the deferred
  multi-layer surface compositing, so do it together with that.

### Central value→colour mapping, two renderers (T. Wager) — IN PROGRESS

**Step 1 done (2026-06-25): the `canlab_colormap` component.**
`CanlabCore/fmridisplay_helper_functions/canlab_colormap.m` is the central value→RGB map.
A value class with static factories `single` / `split` / `solid` / `indexed` and
`from_render_args(args, cmaprange)` (bridges the existing layer render_args), plus:
`map(values)` → N×3 RGB (matches render_blobs: split interpolates minpos→maxpos for
positives and maxneg→minneg for negatives; single ramps mincolor→maxcolor; solid is
constant; indexed rounds-and-indexes; uncoloured values are NaN), `legend_samples(n)` →
one colorbar's worth (split spans zero in a single bar — the fix for the single-range
legend), and `lut(n)` → n×3 table. Tested in
`Unit_tests/visualization/canlab_test_canlab_colormap.m` (14 cases).

**Step 2 done (2026-06-25): surfaces render in true-colour RGB via the central map.**
`render_on_surface` gained an opt-in `'truecolor', <canlab_colormap>` option: after it samples
per-vertex values, it colours each vertex with `map()` (N×3 RGB), compositing over the anatomy
gray (`anatomy_to_rgb`, which preserves curvature shading) for uncoloured vertices, and sets
`FaceVertexCData` to RGB with `CDataMapping` direct. Default behaviour is unchanged when the
option is absent (the `@image_vector/surface` bypass path still uses indexed colour).
`render_layer_surfaces` now builds `canlab_colormap.from_render_args(args, cmaprange)` and passes
it, so the managed display's surfaces match the montage colours (e.g. warm = red→yellow with no
blue negatives, verified) and `removeblobs`/erase still restore the gray. The colorbar legend is
still built by `render_on_surface` from pos/neg colormaps for now.

**Step 3 done (2026-06-25): multi-layer surface compositing.** Now that surfaces are true-colour
RGB, layers composite: `render_layer_surfaces` paints a layer onto the *current* surface colours
(no erase), so a new layer's coloured vertices win and its uncoloured vertices keep what's
underneath (lower layers / gray). render_on_surface saves the anatomy gray **once** (save-once
`UserData`) so erase/`removeblobs` still restore gray. A full recompute (reset to gray, repaint
every layer bottom→top) is `composite_surfaces`, used by `refresh` (hence rethreshold /
set_colormap / set_opacity) and `remove_layer`. Verified: a broad green mask under narrower split
stats shows green where the stats don't cover and split on top; `remove_layer` recomposites the
remainder; `removeblobs` restores gray.
**Step 4 done (2026-06-25): per-layer surface visibility + opacity blending + sign-aware
defaults.**
- Surfaces now **blend** layers by opacity rather than top-wins: `render_on_surface` takes a
  `'truecolor_alpha'` and blends a layer's coloured vertices with what's underneath
  (`a*colour + (1-a)*below`), so a semi-transparent layer lets the layer below show through —
  matching montages. `render_layer_surfaces` passes the layer's `transvalue` as the alpha
  (no more whole-patch `FaceAlpha`).
- **Per-layer visibility:** a layer with `activation_maps{k}.visible == false` is skipped in
  `render_layer_surfaces`/`composite_surfaces`; the controller's Visible toggle sets the flag,
  hides the montage blobs, and recomposites the surfaces.
- **Sign-aware default colormap** (`addblobs`, only when no explicit colour is given): a binary
  mask → solid colour, positive-only → warm, negative-only → cool, mixed +/- → mango split.
  Propagates to both montage and surface via the central map.

**Step 5 done (2026-06-25): addbrain keyword pass-through for `@fmridisplay/surface` (§7.4) +
help/region listing.**
- `@fmridisplay/surface` no longer hard-codes a per-direction `if/elseif` switch. It now parses
  three reserved keyword pairs (`'direction'`, `'orientation'`, `'axes'`), accepts a **bare
  direction token** (`surface(o2, 'thalamus')` == `surface(o2, 'direction', 'thalamus')`), and
  passes any remaining args straight through to `addbrain(dir, extra{:})`. So **any eligible
  addbrain surface/region keyword works in the managed display automatically** — new addbrain
  surfaces need no change here.
- **addbrain owns the view.** addbrain already sets the default (lateral) camera view + lighting
  per keyword, so surface() no longer duplicates that. It only (1) special-cases `'bigbrain
  left/right'` (they map to a *different* addbrain surface than `addbrain('bigbrain left')` would)
  and (2) mirrors the azimuth 180° for `'medial'` (`medial = lateral az + 180`, read back from the
  axes — works for every surface regardless of its lateral convention). Verified the managed
  multi-surface layouts (hires/hcp/freesurfer L+R × lateral/medial) keep **identical** view angles;
  only a couple of standalone lateral azimuths (`inflated right`, `surface right`) move to
  addbrain's surface-specific convention, which was the more correct one.
- **Composites centralized into addbrain (T. Wager's suggestion).** The `brainstem left/right` and
  `caudate left/right` composite handle-sets that surface() used to assemble inline now live in
  `addbrain` as cases (next to `limbic`/`bg`/`*_group`), with their oblique 3D views. surface()
  just passes the keyword through. addbrain's help region list documents them.
- **Help/region listing:** `help fmridisplay/surface` now documents the common directions, the
  bare-token form, the `'orientation'`/`'axes'` options, and points to `help addbrain` for the
  full, current keyword list.
- Tests: added `test_surface_addbrain_passthrough` (bare token + composite multi-patch) and
  `test_surface_medial_flips_azimuth` to the handle suite (now 37/37). Display 4/4 unchanged.

**Step 6 done (2026-06-25): the legend lives in the controller; figure legends default OFF.**
- **In-controller legend.** Each layer panel's colour stripe is now the legend: it is taller and
  carries numeric end labels underneath (`leglo_*`/`legmid_*`/`leghi_*`), read from the layer's
  `cmaprange` and rounded to 1 significant figure. Split (+/-) maps label both extremes with `0`
  in the centre; single ramps label just the two ends. `update_controls_in_place` refreshes the
  labels when `cmaprange` changes (rethreshold/recolor). The opacity/threshold/colors/visible rows
  shifted down a row to make room (panel 195->220 px, fig increment 205->230).
- **Figure legends OFF by default.** `@image_vector/montage` no longer draws the montage-figure
  colorbar by default (pass `'legend'` to opt back in, e.g. for export); the previously no-op
  `'nolegend'` is now effectively the default. `render_layer_surfaces`/`composite_surfaces` take a
  `show_legend` flag (default false -> pass `'nolegend'` to `render_on_surface`), so surfaces get no
  colorbar by default either. (`fmridisplay` constructor + `@image_vector/montage` now accept
  `'legend'`/`'nolegend'` without warning.)
- **Toggle legend fixed + repurposed.** The footer button toggles colorbars on the montage/surface
  FIGURES (for export). ON targets the montage figure explicitly (the old bug: `legend()` drew into
  the controller uifigure via `gcf`, so it never reappeared) and re-renders surfaces with
  `show_legend=true`; OFF removes them. State tracked in the controller's appdata; montage legend
  axes tagged `fmridisp_fig_legend` for clean removal. Verified on->off->on round-trip.
- Tests: +`test_controller_shows_legend_labels`, +`test_montage_figure_legend_off_by_default`;
  updated the surface-legend tests to force legends on (`build_montage_surface_blobs(tc,true)` /
  `composite_surfaces(o,[],true)`) and `test_surface_pulls_in_existing_blobs` to check vertex
  colouring instead of a legend. NOTE: `render_layer_surfaces` is declared in the classdef
  `methods (Hidden)` block, so its new 4-arg signature requires a fresh class load (`clear classes`)
  — verified in a clean MATLAB process; the live MCP session couldn't `clear classes` (lingering
  onCleanup objects), so the surface-legend tests there error with a stale 3-arg dispatch only.

**Bugfix (2026-06-25): mixed isosurface+isocaps surfaces (`coronal_slabs_4`).**
`@fmridisplay/surface` did `if strcmp(h(1).FaceColor, 'interp')`, which errored ("Dot indexing is
not supported for variables of this type") when addbrain returns OLD-STYLE numeric handles
(doubles) rather than graphics objects — as `coronal_slabs_4` does (8 patch handles). It also would
have blanket-grayed the whole set, destroying the isocaps that use `'interp'` to show the anatomy
cross-section. Fixed with a per-object loop using get/set (works on numeric handles too): keep any
`'interp'` colouring, set the rest to flat gray, make all opaque. Test:
`test_surface_mixed_patch_handles`.

**QA pass (2026-06-25): all unit suites green.** `canlab_test_fmridisplay_handle` 41/41,
`canlab_test_display` 4/4, `canlab_test_canlab_colormap` 14/14 (interactive). Verified in a fresh
MATLAB process for the surface-legend / 4-arg-method tests (the live MCP session can't `clear
classes` — lingering onCleanup objects — so classdef changes need a fresh class load there). Three
apparent failures in headless `-batch` are environment artifacts, not regressions: closed-figure
`ishandle` timing (`test_removeblobs_survives_closed_surface_figure`), flaky offscreen GL
(`test_surface_runs_on_thresholded_t`, passes when run alone), and the colormap suite folder not
being on the batch path — all pass in an interactive session.

**Step 7a done (2026-06-25): montage slices routed through the central `canlab_colormap`.**
`render_blobs` no longer re-derives split/single colours inline. It builds one `canlab_colormap`
(`split`/`single`) from the parsed colours + cmaprange before the slice loop, and colours each slice
via a shared `central_map_slice` helper — the SAME mapping surfaces already use, so montage and
surface colours now come from one implementation. Verified pixel-identical to the old renderer on
split/single/solid sample montages (displayed-colour maxdiff: single 0, split 4e-12 float-eps,
solid 0; the only raw-CData differences were out-of-range garbage the old code stored at invisible
alpha=0 voxels, which the central map correctly clamps). Retired the now-dead duplicated colour
logic: `splitcolor_Z_to_slicecdat` and the `cdat2`/`cdatminneg`/`cdatmaxneg` matrices.

**Next:** route the legend (`@fmridisplay/legend` + the controller stripe) through the same central
map (one colorbar spanning zero for split; retires the last copy of the colour math), then it's
fully unified.

The agreed long-term architecture for the colour pipeline. Today the value→colour logic is
duplicated and divergent: `render_blobs` (slices) and `render_on_surface` (surfaces) each
re-implement how a layer's values map to colours for split / single-ramp / solid / index
colormaps, which is exactly why a "warm" or split map looks different (or wrong) on surfaces
vs montages, and why colour fields drift out of sync (e.g. the `legend` `mincolor` bug).

Target design:
1. **One central colour-mapping module.** Given a layer's colormap *type* (single ramp,
   split +/- , solid, indexed) and its parameters (colours, cmaprange/threshold), produce a
   single canonical **value→colour map** (e.g. a function or a resolved lookup table +
   the value range), plus the matching legend spec. This is the single source of truth;
   `set_colormap` updates it and everything else reads it (montage, surface, legend,
   controller), so nothing goes stale.
2. **Pass that map to the renderers, which differ in mechanism.** Slices keep their per-slice
   `surf` objects; surfaces should render in **true-colour RGB** (N×3 `FaceVertexCData`,
   `CDataMapping` direct) computed from the central map — this also unlocks the deferred
   **multi-layer surface compositing** (composite RGB across layers, top wins per vertex) and
   **per-layer surface visibility**, neither of which is cleanly possible with the current
   indexed-colour, single-axes-colormap surface path.
3. **Outcome:** identical colours across montage/surface/legend, no field-sync bugs, and a
   natural home for new colormaps (mango, niivue inferno/viridis, …) defined once.

**Concrete bug this fixes — single-range colorbar legend.** For single-ramp colormaps
(`warm`/`cool`/`winter` = max/min colour), the legend should be **one** colorbar spanning the
**full data range (through zero)**, not a split +/- pair. Today it can't be done cleanly:
`render_blobs` computes a positive-percentile `cmaprange` (e.g. `[2.82 4.86]`) that doesn't
span zero, and `render_on_surface` draws two colorbars (pos+neg) whenever the data has both
signs — so the montage legend range is wrong and the surface legend looks split. Both follow
from the duplicated range/colorbar logic; the central map (one canonical value→colour range
per colormap type, with single-ramp → one full-range bar) is the right place to fix it.

This subsumes and is the right way to do the earlier "unify `render_blobs` + `render_on_surface`",
"composited multi-layer surfaces", and "true-colour RGB surfaces" items.

### Deferred: make the bypass methods return a managed object (§7.4)

`image_vector.surface(t, …)` and `image_vector.orthviews(t, …)` still return raw graphics
handles (unlike `montage(t)`, which returns a managed `fmridisplay`). Plan:

1. **Opt-in flag first, then default** (the §7.4 recipe; ~24 call sites + the walkthroughs
   use `surface_handles = surface(t, …)` as graphics handles, so the default must not
   change yet). In managed mode, mirror `image_vector.montage`:
   `o2 = fmridisplay; o2 = surface(o2, kw…); o2 = addblobs(o2, region(t), 'source_object', t, …)`.
   The `'source_object'` retention makes `rethreshold`/controller work on the result.
2. **Keyword parity via addbrain pass-through (T. Wager's idea).** Rather than
   re-implementing surface types inside `@fmridisplay/surface`, pass the direction keyword
   plus the remaining `varargin` straight through to `addbrain` — which `@fmridisplay/surface`
   already does for the hcp surfaces (`addbrain('hcp inflated right', varargin{:})`). That
   exposes addbrain's full surface catalog to the managed path with almost no new fmridisplay
   code: collapse the large hard-coded `if/elseif dir` switch into a default
   `h = addbrain(dir, varargin{:})` (keeping only the cases that need a specific `view()`/
   layout), and let addbrain validate / error on unknown keywords (it has 69 cases + an
   `otherwise: error('Unknown method.')`).
   - **Cutaways/slabs are already covered by addbrain** (corrects an earlier caveat):
     `addbrain.m:844` dispatches `'cutaway' / 'left_cutaway' / 'right_cutaway' /
     'right_cutaway_x8' / 'coronal_slabs' / 'coronal_slabs_4' / 'coronal_slabs_5' /
     '*_insula_slab' / 'accumbens_slab'` to `canlab_canonical_brain_surface_cutaways(meth,
     varargin{:})`. So the addbrain pass-through covers cutaways and coronal slabs in ONE
     mechanism — no separate cutaway step needed. (Only nuance: the fully *parametric*
     `surface(r, 'cutaway', 'ycut_mm', -30)` form in `@region/surface` uses the older
     `surface_cutaway` engine, a slightly different beast than addbrain's
     `canlab_canonical_brain_surface_cutaways`; the named keywords like `coronal_slabs_4`
     are the supported managed path.)
3. **The real `@fmridisplay/surface` work** (why this isn't a quick piecemeal change): the
   method currently sets `dir` only from the explicit `'direction'` keyword and never treats
   a *bare* token as a direction (except the `foursurfaces*` special-case). To accept bare
   addbrain keywords like `surface(o2, 'coronal_slabs_4')`, the option parsing must learn to
   pull a direction token out of `varargin`, strip the fmridisplay control opts
   (`direction`/`orientation`/`axes`), forward the rest to `addbrain`, and apply sensible
   post-render `view()`/lighting (some surfaces set their own). Wants per-surface-type
   testing, so do it here, not by hand-adding cases.
4. **Default surface = `left_cutaway`** (requested): change the default `dir` in
   `@fmridisplay/surface` (currently `'hires right'`) and the no-keyword default of
   `image_vector.surface` (currently `@region/surface`'s default). This is back-compat-
   affecting (changes what `surface(o2)` / `surface(t)` shows with no keyword, incl. in
   walkthroughs/scripts), so make it a deliberate decision in this phase rather than a
   standalone tweak.

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
