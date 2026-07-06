# 1. Getting started

[Back to index](index.md) · Next: [2. Montages and slices →](02_montages.md)

## Objects you'll visualize

CANlab visualization is organized around a few object classes. Knowing which one you hold tells you which display methods are available:

| Object | What it holds | Typical display |
|--------|---------------|-----------------|
| **`fmri_data`** | A set of brain images (data). | `plot()` QC overview, `mean()`, then montage/surface. |
| **`statistic_image`** | A statistic map with p‑values; supports `threshold()`. | `orthviews`, `montage`, `surface`. |
| **`region`** | Contiguous suprathreshold blobs, as analysis units. | `montage(r, 'regioncenters')`, `table()`. |
| **`atlas`** | A labeled parcellation (many named regions). | `montage`, `isosurface`, unique-color rendering. |
| **`fmridisplay`** | A display container: montages + surfaces + their handles. | `addblobs`, `removeblobs`, `controller`. |

Almost every method below exists on several of these classes, so `montage`, `surface`, and `orthviews` "just work" whatever you're holding.

## Load a dataset

The examples throughout this walkthrough use the bundled `emotionreg` sample dataset, so they run without downloading anything:

```matlab
obj = load_image_set('emotionreg', 'noverbose');   % 30 contrast images (fmri_data)
```

`load_image_set` also fetches larger published datasets by keyword (e.g. `'kragel18_alldata'`); see the [datasets tutorial](https://canlab.github.io/_pages/canlab_help_2c_loading_datasets/canlab_help_2c_loading_datasets.html).

## A quick quality-control look: `plot`

`plot(obj)` on an `fmri_data` object gives a six-panel QC overview (data matrix, covariance/correlation, histogram, global signal, and a Mahalanobis outlier plot) plus orthviews of the mean. Run it as a first check on any dataset:

```matlab
plot(obj);      % interactive QC figure + outlier report in the console
```

## Make a result to display

Most of the walkthrough displays a statistic map. Compute one with a voxelwise t‑test and threshold it:

```matlab
t = ttest(obj);                 % statistic_image
t = threshold(t, .05, 'unc');   % keep p < .05 uncorrected
```

## First look: `orthviews`

The fastest way to inspect a map interactively is `orthviews` — three planes with a movable crosshair. Click to move the crosshair through the volume:

```matlab
orthviews(t);
```

![orthviews of the thresholded map](figures/01_orthviews.png)

That's enough to start. From here the walkthrough builds up the full toolkit:

- **[2. Montages and slices](02_montages.md)** — publication slice arrays, per-region close-ups, customization.
- **[3. Surfaces and 3‑D rendering](03_surfaces.md)** — cortical surfaces, cutaways, subcortical structures, isosurfaces.
- **[4. Colors and colormaps](04_colormaps.md)** — the shared value→color pipeline; making maps comparable.
- **[5. The display controller](05_controller.md)** — the interactive panel, and its command-line twin.
- **[6. Atlases and regions](06_atlases_and_regions.md)** — parcellations in unique colors, on slices, surfaces, and in 3‑D.

---

[Back to index](index.md) · Next: [2. Montages and slices →](02_montages.md)
