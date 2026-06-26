# canlab_niivue

A lightweight, dependency-free **web orthviews** for CANlab neuroimaging results, built on
[NiiVue](https://github.com/niivue/niivue). It renders a NIfTI underlay plus a statistic/mask
overlay in the browser — multiplanar / axial / coronal / sagittal / 3-D-render layouts,
transparent sub-threshold voxels (blobs over anatomy), colormap / threshold / opacity controls,
and a live MNI-coordinate + value readout. An atlas (attached by default) names the region under
the crosshair and highlights **just that region** — shown as an Outline (default), Shaded fill, or
Off via the "Atlas region" dropdown. No build step, no CDN, works offline.

## 📖 Full documentation

**See the guide (with a live embedded demo): [`docs/canlab_niivue_guide.md`](../../../docs/canlab_niivue_guide.md)**
— MATLAB usage, HTML embedding, report integration, and the full option reference.

## Quick start (MATLAB)

```matlab
t = threshold(ttest(load_image_set('emotionreg')), .005, 'unc', 'k', 10);
canlab_niivue(t);                                  % standalone .html (atlas on by default)
canlab_niivue(t, 'outdir', 'myfig', 'standalone', false);  % folder bundle (auto-served)
canlab_niivue(t, 'atlas', 'glasser');              % use a different atlas as the source
canlab_niivue(t, 'noatlas');                        % disable the atlas region readout
```

The atlas (`'atlas', obj | keyword | true`) is loaded as a hidden integer-label layer that names
the region under the crosshair and highlights only that region (Outline / Shaded / Off dropdown).
It is **on by default** in both modes — the label volume gzips to ~74 KB, so the standalone
`.html` cost is small. `'noatlas'` (or `'atlas','none'`) disables it.

> **Folder mode is auto-served.** Folder bundles use ES-module imports, which browsers block on
> `file://` (CORS) — opening one directly shows a blank page. `canlab_niivue` therefore starts a
> short-lived local `python3 -m http.server` and opens the `http://localhost` URL for you.
> Standalone mode has no such restriction (single self-contained file).

## Embed in your own HTML page (folder mode)

Copy `assets/niivue.js`, `assets/canlab_niivue_viewer.js`, `assets/canlab_niivue.css` next to your page:

```html
<link rel="stylesheet" href="canlab_niivue.css">
<div id="canlab-niivue-controls"></div>
<div id="canlab-niivue-canvas-container"><canvas id="gl"></canvas></div>
<div id="canlab-niivue-readout"></div>
<script type="module">
  import { canlabNiivue } from "./canlab_niivue_viewer.js";
  canlabNiivue("gl", { underlay: "./underlay.nii.gz", overlay: "./stat.nii.gz", colormap: "inferno" });
</script>
```

Serve over `http(s)` (ES-module imports are blocked on `file://`). The standalone mode produced by
`canlab_niivue(t)` has no such restriction — it is a single self-contained file.

## Files

- `canlab_niivue.m` — MATLAB bridge (object → viewer page).
- `assets/` — vendored `niivue.js` (pinned 0.57.0), `canlab_niivue_viewer.js`, `canlab_niivue.css`, HTML template.
- `sample/` — the emotionreg one-sample t-test demo (`emotionreg_ttest.html` standalone; `index.html` folder mode).

NiiVue is pinned at 0.57.0; to update, re-download `dist/index.js` from unpkg into `assets/niivue.js`.
