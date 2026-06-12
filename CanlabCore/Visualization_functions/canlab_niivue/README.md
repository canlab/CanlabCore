# canlab_niivue

A lightweight, dependency-free **web orthviews** for CANlab neuroimaging results, built on
[NiiVue](https://github.com/niivue/niivue). It renders a NIfTI underlay plus a statistic/mask
overlay in the browser — multiplanar / axial / coronal / sagittal / 3-D-render layouts,
transparent sub-threshold voxels (blobs over anatomy), colormap / threshold / opacity controls,
and a live MNI-coordinate + value readout. No build step, no CDN, works offline.

## 📖 Full documentation

**See the guide (with a live embedded demo): [`docs/canlab_niivue_guide.md`](../../../docs/canlab_niivue_guide.md)**
— MATLAB usage, HTML embedding, report integration, and the full option reference.

## Quick start (MATLAB)

```matlab
t = threshold(ttest(load_image_set('emotionreg')), .005, 'unc', 'k', 10);
canlab_niivue(t);                                  % standalone .html, opened in browser
canlab_niivue(t, 'outdir', 'myfig', 'standalone', false);  % folder bundle
```

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
