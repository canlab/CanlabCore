# `canlab_niivue` — interactive web viewer for brain images

[← back to Object methods index](Object_methods.md)

A lightweight, dependency-free **web orthviews** for CANlab results, built on
[NiiVue](https://github.com/niivue/niivue). It renders a NIfTI anatomical underlay plus a
statistic/mask overlay in any browser, with multiplanar / axial / coronal / sagittal / 3-D-render
layouts, transparent sub-threshold voxels (blobs only over anatomy), a colormap dropdown, a
dynamic display-threshold slider, overlay/underlay opacity sliders, and a live MNI-coordinate +
value readout that updates as you move the crosshair. There is **no plugin, no server-side
compute, and no build step** — it runs entirely in the browser and works offline.

It is the web counterpart to `orthviews` / `canlab_orthviews`: the same point-and-click slice
view, but as a portable `.html` file you can email, host, or drop into an HTML report.

---

## Live demo — emotionreg one-sample t-test

<iframe src="niivue_demo/emotionreg_ttest.html" width="100%" height="480"
        style="border:1px solid #d4d8dd; border-radius:6px;" loading="lazy"></iframe>

*One-sample t-test on the CANlab `emotionreg` sample dataset (p &lt; .005 uncorrected, k ≥ 10),
shown as warm/cool blobs over a 2&nbsp;mm MNI underlay. Use the **layout buttons** to switch views,
the **Colormap** dropdown to recolor, the **Threshold** slider to raise the display threshold, the
**Overlay/Underlay** sliders to fade layers, and click anywhere to move the crosshair (its MNI
coordinate and t-value print below the canvas). If the frame above is blank in your environment,
[open the demo directly](niivue_demo/emotionreg_ttest.html).*

This entire viewer — the NiiVue library, the image data, and the styling — is inlined into a
single self-contained `emotionreg_ttest.html` (~3.6&nbsp;MB). It was generated with one MATLAB
call (below).

---

## Quick start (MATLAB)

```matlab
dat = load_image_set('emotionreg');                 % 30 contrast images (fmri_data)
t   = threshold(ttest(dat), .005, 'unc', 'k', 10);  % statistic_image
canlab_niivue(t);                                   % write + open a standalone .html
```

`canlab_niivue` accepts any CANlab image-like object — `fmri_data`, `statistic_image`, `atlas`,
or a path to a `.nii`/`.nii.gz` — as the **overlay**, and supplies a default MNI anatomical
**underlay** automatically. For a `statistic_image` it writes the thresholded map by default, and
auto-derives the color limits from the suprathreshold values (so `cal_min` becomes the threshold
and below-threshold voxels are transparent).

```matlab
% Folder bundle instead of a single file (page + copied assets + data/*.nii.gz):
canlab_niivue(t, 'outdir', 'my_report', 'standalone', false);

% Fixed color limits, a different colormap, a custom underlay, no auto-open:
canlab_niivue(t, 'colormap', 'hot', 'cal_min', 3, 'cal_max', 8, ...
                 'underlay', 'mni152_1mm', 'noopen');
```

---

## How it gets from MATLAB to the browser

`canlab_niivue` writes the overlay (and underlay, if it is an object) to `.nii.gz`, then fills an
HTML template that boots the vendored NiiVue bundle and the CANlab viewer wrapper
(`canlabNiivue`). The object's data lives in `.dat`; NiiVue draws the slices on the GPU (WebGL2)
and handles the click-to-navigate crosshair, so MATLAB is only involved at export time — the page
is fully interactive on its own afterwards.

---

## Two output modes

| Mode | What you get | Use when |
|------|--------------|----------|
| **Standalone** (default) | One self-contained `.html` with the NiiVue library, the image data (base64), the viewer JS, and CSS all inlined. | Sharing a single file, emailing, or embedding one viewer in a report. A fine (≤2&nbsp;mm) underlay is auto-downsampled to ~2&nbsp;mm so the file stays small. |
| **Folder** (`'standalone', false`) | `index.html` + copies of `niivue.js`, `canlab_niivue_viewer.js`, `canlab_niivue.css`, and `data/*.nii.gz`. | Serving from a web folder, or many viewers that share one cached NiiVue bundle. Must be served over `http(s)` (ES-module imports are blocked on `file://`). |

---

## Embed the viewer in your own HTML page

Copy `niivue.js`, `canlab_niivue_viewer.js`, and `canlab_niivue.css` (from
`CanlabCore/Visualization_functions/canlab_niivue/assets/`) next to your page, then:

```html
<link rel="stylesheet" href="canlab_niivue.css">

<div id="canlab-niivue-controls"></div>
<div id="canlab-niivue-canvas-container"><canvas id="gl"></canvas></div>
<div id="canlab-niivue-readout"></div>

<script type="module">
  import { canlabNiivue } from "./canlab_niivue_viewer.js";
  canlabNiivue("gl", {
    underlay: "./underlay.nii.gz",   // optional anatomical
    overlay:  "./my_stat_map.nii.gz",// stat/mask
    colormap: "inferno",
    cal_min: 3, cal_max: 8           // optional; below cal_min is transparent
  });
</script>
```

The viewer builds the layout buttons, colormap dropdown, sliders, and readout into the three
`div`s above. Serve the folder over `http(s)`; opening via `file://` blocks the ES-module import.

---

## Embed inside an HTML report (with tables, stats, multiple viewers)

The most robust way to place a viewer next to other report content is an **`<iframe>`** around a
standalone page — it isolates NiiVue's canvas, element ids, and module scope, so several viewers
coexist on one report page with no collisions:

```html
<h2>Region X</h2>
<table> … your stats … </table>
<iframe src="region_x_viewer.html" width="720" height="480"
        style="border:1px solid #d4d8dd; border-radius:6px;"></iframe>
```

Generate each viewer while assembling the report:

```matlab
regions = {t_regionX, t_regionY};
for i = 1:numel(regions)
    canlab_niivue(regions{i}, 'outdir', reportdir, ...
        'fname', sprintf('region_%d.html', i), 'noopen');
    % then write <iframe src="region_i.html" …> into your report HTML
end
```

For a report with **many** viewers, prefer folder mode (one shared 2.7&nbsp;MB NiiVue bundle)
rather than repeating it inside every standalone file.

---

## Key MATLAB options (`canlab_niivue(obj, ...)`)

| Option | Default | Description |
|--------|---------|-------------|
| `'underlay'` | lab MNI template | Anatomical underlay: filename, keyword (e.g. `'mni152_1mm'`), or image object. `'none'` shows the overlay alone. |
| `'standalone'` | `true` | Single self-contained file vs. folder bundle. |
| `'colormap'` | `'inferno'` | Overlay positive colormap. |
| `'colormapNegative'` | `'winter'` | Overlay negative colormap (hot/cool split). |
| `'cal_min'` / `'cal_max'` | auto | Display range; `cal_min` is also the transparency cutoff. Auto-derived from suprathreshold magnitude. |
| `'opacity'` | `0.8` | Overlay opacity (also a live slider). |
| `'thresh'` | `true` for `statistic_image` | Write the thresholded map (zero out non-significant voxels). |
| `'underlay_2mm'` | `true` | In standalone mode, downsample a finer underlay to ~2&nbsp;mm to keep the file small. |
| `'outdir'` / `'fname'` | `./canlab_niivue_output` / `index.html` | Where to write. |
| `'open'` / `'noopen'` | open | Open the page in your browser. |

See `help canlab_niivue` for the full list.

## Key viewer options (`canlabNiivue("gl", config)`)

| Key | Default | Description |
|-----|---------|-------------|
| `underlay` / `overlay` | — | URL string or `{ base64, name }` payload. |
| `colormap` / `colormapNegative` | `"inferno"` / `"winter"` | Overlay positive / negative colormaps. |
| `cal_min` / `cal_max` | auto | Display range + transparency cutoff. |
| `opacity` | `0.8` | Overlay opacity. |
| `sliceType` | `"multiplanar"` | `multiplanar`/`axial`/`coronal`/`sagittal`/`render`. |
| `showRenderInMultiplanar` | `true` | Show the 3-D render tile beside the slices. |
| `controls` / `allowLoadOverlay` | `true` / `true` | Layout buttons + colormap dropdown + "Load overlay…" picker. |
| `showOpacity` / `showReadout` | `true` / `true` | Threshold + opacity sliders / crosshair readout. |
| `showOrientationLabels` | `false` | A/P/S/I/L/R letters (off avoids collisions at small sizes). |
| `crosshairWidth` / `textHeight` | `0.5` / `0.04` | Crosshair line width / colorbar font size. |

---

## Interactive controls in the page

- **Layout buttons** — Axial / Coronal / Sagittal / Multiplanar / 3D Render.
- **Colormap dropdown** — recolor the overlay live.
- **Load overlay…** — open a different `.nii`/`.nii.gz` from disk into the viewer.
- **Threshold slider** — raise/lower the display threshold (`cal_min`); below it is transparent.
- **Overlay / Underlay sliders** — fade each layer's opacity.
- **Crosshair** — click any slice; the MNI coordinate and overlay value print below the canvas.

---

## Files and updating NiiVue

The component lives in `CanlabCore/Visualization_functions/canlab_niivue/`:
`canlab_niivue.m` (the MATLAB bridge), `assets/` (the vendored `niivue.js`, the viewer JS, CSS,
and HTML template), and `sample/` (this emotionreg demo). NiiVue is pinned at **0.57.0**; to
update, re-download `dist/index.js` from unpkg into `assets/niivue.js` and re-test.

---

## See also

- [`fmri_data.montage`](individual_functions/fmri_data_montage.md) — static slice montage on a canonical underlay
- [`canlab_results_fmridisplay`](individual_functions/canlab_results_fmridisplay.md) — montage/surface scaffolds (`fmridisplay`)
- `surface`, `isosurface` (image_vector methods) — 3-D cortical-surface rendering in MATLAB
- `orthviews` / `canlab_orthviews` — SPM-based interactive orthviews in MATLAB
- [Object methods index](Object_methods.md)
