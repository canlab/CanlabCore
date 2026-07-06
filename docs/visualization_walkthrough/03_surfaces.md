# 3. Surfaces and 3‑D rendering

[← 2. Montages](02_montages.md) · [Back to index](index.md) · Next: [4. Colors and colormaps →](04_colormaps.md)

Surface renderings localize activity on the cortical sheet and in 3‑D subcortical structures, and make patterns comparable across studies. Most data objects (`fmri_data`, `statistic_image`, `region`) have a `surface()` method, and the underlying engine — `render_on_surface` — paints values onto any set of surface handles you build with `addbrain` or `isosurface`.

```matlab
obj = load_image_set('emotionreg', 'noverbose');
t   = threshold(ttest(obj), .05, 'unc');
```

---

## 3.1 The `surface()` method

Called with a layout keyword, `surface()` builds the surfaces and renders the map in one step. **`foursurfaces`** gives the standard lateral + medial view of both hemispheres:

```matlab
sh = surface(t, 'foursurfaces', 'noverbose');   % returns surface handles
```

![foursurfaces](figures/03_foursurfaces.png)

Positive values are warm, negative cool, with a split colorbar in t units. The returned handles `sh` can be recolored, made transparent, or rotated afterward.

Called with no layout, `surface(t)` produces a **cutaway** that reveals medial cortex and subcortical structures together:

```matlab
surface(t);
```

![surface cutaway](figures/03_surface_default.png)

---

## 3.2 Recoloring with `render_on_surface`

`render_on_surface` repaints existing surface handles — change the colormap, the color limits, or emphasize one sign. Here the same `foursurfaces` handles are recolored with a perceptual `summer` colormap:

```matlab
render_on_surface(t, sh, 'colormap', 'summer');
```

![foursurfaces recolored with summer](figures/03_foursurfaces_summer.png)

Key `render_on_surface` options:

| Option | Effect |
|--------|--------|
| `'colormap', name` | Any MATLAB colormap (`'hot'`, `'summer'`, `'winter'`, …) or an *n×3* matrix. |
| `'clim', [lo hi]` | Value→color limits. Use `[-3 3]` for a symmetric signed map, or `[2 6]` to emphasize positives. |
| `'splitcolor', {…}` | Custom warm/cool split colors. |
| `'srcdepth'`, `'sourcespace'`, `'targetsurface'` | Accurate volume→surface projection (see note below). |

> **Accurate projection.** MNI *volume* voxels don't map one‑to‑one onto standard *surface* vertices. For publication-grade cortical projection, pass `'sourcespace'` (e.g. `MNI152NLin2009cAsym`) and `'targetsurface'` (`fsLR_32k`, `fsaverage_164k`) so CANlab uses precomputed vertex mappings and samples across the cortical ribbon. See `help render_on_surface`.

---

## 3.3 Cutaways and slabs

`surface()` accepts many prepackaged 3‑D cutaways. **`coronal_slabs_4`** cuts the brain into four coronal slabs so you can see activation on interior cross-sections:

```matlab
surface(t, 'coronal_slabs_4');
```

![coronal slabs](figures/03_coronal_slabs_4.png)

Other cutaway keywords: `'left_cutaway'`, `'right_cutaway'`, `'left_insula_slab'`, `'right_insula_slab'`, `'accumbens_slab'`, `'coronal_slabs'`, `'coronal_slabs_5'`. Each can be passed to `surface()` (to render a map) or to `addbrain()` (to get bare surfaces).

---

## 3.4 Subcortical close-ups with `addbrain`

`addbrain` is the master surface builder: dozens of keywords return handles to cortical or subcortical anatomy that you compose and then color with `render_on_surface`. This builds a subcortical scene — thalamic nuclei, brainstem, amygdala — and renders the t‑map onto it:

```matlab
sh = addbrain('thalamus_group');
sh = [sh addbrain('brainstem')];
sh = [sh addbrain('amygdala')];
render_on_surface(t, sh, 'clim', [-3 3]);
set(sh, 'FaceAlpha', .9);
view(135, 12);  lightRestoreSingle;  camzoom(1.3);
```

![subcortical close-up](figures/03_subcortical_closeup.png)

`addbrain` knows individual structures (`'thalamus'`, `'amygdala'`, `'hippocampus'`, `'brainstem'`, `'pag'`, `'nucleus accumbens'`, `'caudate'`, `'putamen'`, `'cerebellum'`, …), composites (`'limbic'`, `'basal ganglia'`, `'brainstem_group'`, `'thalamus_group'`), and cortical surfaces (`'hires left'`, `'inflated left'`, …). Handles are standard MATLAB graphics: set `FaceAlpha`, rotate with `view()`, zoom with `camzoom`, and re-fix lighting with `lightRestoreSingle`. `addbrain('eraseblobs', sh)` resets a surface to gray so you can render a different contrast on the same scene.

---

## 3.5 Isosurfaces from atlases and images

`isosurface()` turns any image or atlas into 3‑D mesh blobs. Here every nucleus of a thalamic atlas becomes a surface; then we render the t‑map onto them:

```matlab
atl = load_atlas('thalamus');
sh  = isosurface(atl);
axis image vis3d off;  material dull;  view(210, 20);  lightRestoreSingle
render_on_surface(t, sh, 'colormap', 'hot');
```

Atlas nuclei as bare isosurfaces | t‑map rendered onto them
:---:|:---:
![thalamus isosurface](figures/03_isosurface_thalamus.png) | ![thalamus isosurface rendered](figures/03_isosurface_thalamus_rendered.png)

`isosurface` also builds custom cutaways from any anatomical image via `'thresh'`, `'nosmooth'`, and `'xlim'/'ylim'/'zlim'` bounds — useful for a bespoke publication surface. See `help image_vector.isosurface`.

---

[← 2. Montages](02_montages.md) · [Back to index](index.md) · Next: [4. Colors and colormaps →](04_colormaps.md)
