# Loading and using image sets with `load_image_set` — how-to

`load_image_set` is the central registry that turns a short **keyword** into a
loaded brain-image object. It resolves the files on your MATLAB path (mostly from
the CANlab `Neuroimaging_Pattern_Masks` repo), loads them, and hands you back an
object ready to analyze:

- most keywords return an **`fmri_data`** object;
- surface / grayordinate (CIFTI) map sets return an **`fmri_surface_data`** object.

This page shows how to discover what is available, load it, and use it —
applying signatures, computing similarity to networks and topics, and dual
regression.

Prerequisites: CanlabCore **and** `Neuroimaging_Pattern_Masks` on your path
(`canlab_toolbox_setup`). SPM is needed for volume I/O.

---

## 1. Discover what you can load: `load_image_set('list')`

Start here. With no data, just run:

```matlab
load_image_set('list')
```

This prints a series of **categorized tables** to the screen, each under a
descriptive title:

| Table | What's in it | Returns |
|---|---|---|
| **Multivariate signatures** | Predictive patterns (NPS, SIIPS, PINES, VPS, …) with domain flags (pain, negemo, reward, …) | `fmri_data` |
| **Person-level datasets** | Subject-level sample data (`emotionreg`, `bmrk3`, `kragel270`, …) | `fmri_data` |
| **Network, ICA & topic maps** | `bucknerlab`, `neurosynth_topics_fi/ri`, `bgloops`, HCP group-ICA, … | `fmri_data` / `fmri_surface_data` |
| **Surface / grayordinate (CIFTI) map sets** | `hcp_ica15/25/50`, `spectral_bases` | `fmri_surface_data` |
| **Gradients & spatial-basis maps** | `transcriptomic_gradients`, `marg`, `mito_maps` | `fmri_data` |
| **Meta-analysis, receptor & curated-domain** | `emometa`, `pet` (Hansen), `kragelemotion`, `kragelschemas` | `fmri_data` |

`'list'` also **returns a struct** whose fields are those tables, so you can
query them in code:

```matlab
tmp = load_image_set('list');
tmp.surface                         % table of CIFTI map sets
tmp.signatures.keyword              % all signature keywords
tmp.signatures(tmp.signatures.pain==1, :)   % just the pain-related signatures
```

Load anything from the list by its keyword: `obj = load_image_set('<keyword>')`.

---

## 2. Load a person-level dataset

```matlab
data = load_image_set('emotionreg');   % ~30 subjects, reappraise-vs-look contrasts
descriptives(data);                     % quick QC
t = ttest(data);                        % one-sample t across subjects -> statistic_image
```

`data` is a standard `fmri_data` object — use it with `ttest`, `predict`,
`regress`, `apply_atlas`, `montage`, etc.

---

## 3. Apply a signature to your data

A **signature** is a weight map; applying it computes a per-image score (dot
product = "pattern expression", or a correlation). Load the signature with
`load_image_set` and apply it with `apply_mask`:

```matlab
data = load_image_set('emotionreg');
pines = load_image_set('pines');                 % negative-emotion signature

% Dot-product pattern expression (one value per image/subject):
pe  = apply_mask(data, pines, 'pattern_expression', 'ignore_missing');

% Or the spatial correlation with the pattern instead of the dot product:
r   = apply_mask(data, pines, 'pattern_expression', 'ignore_missing', 'correlation');
```

To apply **several** signatures and compare, load them as a set and loop, or use
`image_similarity_plot` (next section), which is built for map sets.

---

## 4. Compare your data to networks and topics (similarity)

`image_similarity_plot` computes the similarity (correlation or cosine) between
your image(s) and every map in a **map set**, and plots a wedge/polar summary. The
map set can be a `load_image_set` keyword or an `fmri_data` you pass in:

```matlab
t = ttest(load_image_set('emotionreg'));

% Similarity to the 7 Yeo/Buckner resting-state networks:
stats = image_similarity_plot(t, 'bucknerlab', 'cosine_similarity');

% Similarity to Neurosynth topic maps (reverse inference):
topics = load_image_set('neurosynth_topics_ri');
stats  = image_similarity_plot(t, 'mapset', topics, 'cosine_similarity', 'plotstyle', 'polar');
```

`stats.r` (or `.cosine_sim`) holds the per-map similarity values you can pull into
a table or figure of your own.

---

## 5. Dual regression with networks / ICA components

**Dual regression** takes a set of group spatial maps (e.g. ICA components or
networks) and, for each subject's **4-D time series**, (1) spatially regresses the
maps into the data to get component time courses, then (2) temporally regresses
those time courses back into the data to get subject-specific spatial maps.

```matlab
% Group maps: HCP resting-state ICA (surface) or a volumetric network set.
gmaps = load_image_set('bucknerlab');          % fmri_data (7 networks), or:
% gmaps = load_image_set('hcp_ica25');         % fmri_surface_data (25 components)

subj  = fmri_data('sub-01_task-rest_bold.nii.gz');   % your 4-D time series

[spatial_maps, timecourses, tmaps] = dual_regression(gmaps, subj);
% spatial_maps : subject-specific spatial map per component (fmri_data)
% timecourses  : component x timepoint matrix
% tmaps        : z-scored spatial maps
```

`dual_regression` resamples the data into the group-map space for you. For the
volumetric FSL-style engine directly, see `dual_regression_fsl`.

---

## 6. Surface / grayordinate (CIFTI) map sets

Surface map sets load natively as `fmri_surface_data` (no external toolbox):

```matlab
ica = load_image_set('hcp_ica25');   % 25 HCP group-ICA components, 91k grayordinates
ica                                  % fmri_surface_data, size(ica.dat) = [91282 x 25]

o2 = surface(get_wh_image(ica, 1));  % render component 1 on the cortical surface
o2 = set_colormap(o2, 'colormap', hot(256));

bases = load_image_set('spectral_bases');   % 200 spectral basis functions
```

Everything on the [`fmri_surface_data` how-to](fmri_surface_data_howto.md) applies:
`surface`, `montage` (subcortex on slices), `apply_parcellation`, `threshold`,
`vol2surf` / `surf2vol`, group analyses, etc.

---

## New keywords added

| Keyword(s) | Collection | Returns |
|---|---|---|
| `hcp_ica15` / `hcp_ica25` / `hcp_ica50` (aliases `hcp15`…) | HCP resting-state group-ICA components (15/25/50), 91k grayordinates | `fmri_surface_data` |
| `spectral_bases` (alias `hcp_bases`) | 200 spectral (Laplacian eigenmap) basis functions, 91k | `fmri_surface_data` |
| `mito_maps` (alias `mito`) | Mitochondrial energetic-capacity maps (CI, CII, CIV, MitoD, MRC, TRC) | `fmri_data` |

## See also

- [`fmri_data` methods](../fmri_data_methods.md), [`fmri_surface_data` methods](../fmri_surface_data_methods.md)
- [`load_atlas`](../atlases_regions_and_patterns.md) — the parallel registry for **parcellations** (run `load_atlas('list')` to see keywords)
- `image_similarity_plot`, `apply_mask`, `dual_regression`
