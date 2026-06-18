# Extracting and visualizing ROI / atlas data

A practical, copy-pasteable guide to pulling region-of-interest (ROI) and
atlas/parcel summaries out of CANlab `fmri_data` objects and plotting them.

All code uses the object-oriented methods that are the current API of record:
`extract_roi_averages`, `apply_parcellation` (and its atlas wrapper
`@atlas/extract_data`), plus the three plotting functions
`barplot_columns`, `lineplot_columns`, and `line_plot_multisubject`.

---

## Setup

Add CanlabCore (with subfolders) and make sure SPM is on the path. SPM is used
for all image I/O, so nothing below will run without it.

```matlab
addpath(genpath('/Users/f003vz1/Documents/GitHub/CanlabCore/CanlabCore'))

% SPM (SPM12 or later) must already be on your path, e.g.:
%   addpath('/path/to/spm')
% Verify:
assert(~isempty(which('spm_vol')), 'SPM is not on the MATLAB path.')
```

The sample atlases and datasets are resolved by keyword from the companion
`Neuroimaging_Pattern_Masks` repo; if `load_atlas`/`load_image_set` returns
empty, that repo is not on your path.

---

## Section A — One ROI mean per image, then a bar plot

We load the blessed test dataset (`'emotionreg'`: 30 single-subject contrast
images from Wager et al. 2008) and extract the mean within a single ROI. Here
the ROI is one region of a loaded atlas, but the same call works with any mask
filename, `region`, `atlas`, or `fmri_mask_image`.

```matlab
% 30 contrast images (one per subject), [voxels x 30]
imgs = load_image_set('emotionreg', 'noverbose');

% Pick a single ROI from an atlas. canlab2024 labels its regions by anatomical
% sub-part (e.g. amygdala sub-nuclei MTL_AMY_*), so match the prefix and
% 'flatten' the sub-regions into one ROI. (Friendly names like 'Amygdala'
% match nothing in this atlas -- always check atl.labels.)
atl = load_atlas('canlab2024', 'noverbose');
roi = select_atlas_subset(atl, {'MTL_AMY'}, 'flatten');   % one bilateral amygdala ROI

% Extract the mean within the ROI for every image.
% Returns a region-object array; with a single-region mask, one element.
cl = extract_roi_averages(imgs, roi, 'noverbose');

% cl(k).dat is [n_images x 1]: one ROI mean per image (subject).
roi_means = cl(1).dat;        % 30 x 1
size(roi_means)               % [30 1]
```

`cl(k).dat` has one row per image (30 here, matching `size(imgs.dat,2)`), and
`cl(k).all_data` is the full `[images x voxels]` matrix if you want the
voxelwise values.

Visualize the per-image values with `barplot_columns`. Its input is one column
per "bar" (here a single column = one ROI), rows are observations (subjects).
It plots the mean +/- SE, a violin, and individual points, and runs a t-test of
the column mean vs. zero.

```matlab
% One bar (the ROI), 30 subjects. statstable has Name/Mean_Value/Std_Error/T/P/Cohens_d.
[handles, dat, xdat, statstable] = barplot_columns(roi_means, ...
    'names', {'Amygdala'}, 'title', 'Emotion reg: ROI mean per subject', ...
    'nofig');

disp(statstable)
```

If you have several ROIs and want a bar per ROI, build a
`[subjects x ROIs]` matrix by horizontally concatenating each region's `.dat`.
Here we keep the amygdala and hippocampus sub-regions separate (no `'flatten'`),
giving one bar per sub-region:

```matlab
roi_multi = select_atlas_subset(atl, {'MTL_AMY', 'MTL_Hipp'});   % 14 sub-regions
cl_multi  = extract_roi_averages(imgs, roi_multi, 'noverbose');   % one element per region
M = cat(2, cl_multi.dat);                                         % [subjects x nROIs]
barplot_columns(M, 'names', roi_multi.labels, 'nofig', 'dolines');
```

---

## Section B — All parcel means from an atlas, then a line plot

To summarize an atlas in one shot, use `apply_parcellation`. It returns a full
`[images x parcels]` matrix of weighted means (one column per original parcel),
which is exactly the "atlas means" data matrix.

```matlab
imgs = load_image_set('emotionreg', 'noverbose');
atl  = load_atlas('canlab2024', 'noverbose');

% parcel_means is [n_images x n_parcels]. Pass the atlas (not an fmri_data) so
% integer labels resample correctly with nearest-neighbor interpolation.
parcel_means = apply_parcellation(imgs, atl);

size(parcel_means)            % [30  n_parcels]
n_parcels = size(parcel_means, 2);
labels    = atl.labels;       % 1 x n_parcels cell of region names
```

The atlas object also exposes the same computation as a method, which resamples
the data into the atlas space first:

```matlab
% Equivalent atlas-method entry point (delegates to apply_parcellation):
parcel_means2 = extract_data(atl, imgs);     % [n_images x n_parcels]
```

Visualize a mean trajectory across parcels (or across conditions) with
`lineplot_columns`. Its input is `[observations x variables]`; the x-axis is the
columns, and it draws a single line of the column means with error bars.

To plot the group mean across the first several parcels:

```matlab
nshow = min(20, n_parcels);
out = lineplot_columns(parcel_means(:, 1:nshow), ...
    'color', [.2 .2 .8], 'markerfacecolor', [.4 .4 1], 'within');

% out.m = column (parcel) means, out.ste = within-subject SE, out.CI95 = 95% CI
set(gca, 'XTick', 1:nshow, 'XTickLabel', labels(1:nshow));
xtickangle(45); ylabel('Mean contrast'); title('Parcel means across parcels');
```

`lineplot_columns` is the right tool when columns are *levels of one factor*
(e.g. ascending stimulus intensities). For that case, build a
`[subjects x conditions]` matrix where each column is a condition and pass it
directly.

---

## Section C — Multi-subject: ROI mean across temperature levels

The `bmrk3` dataset (`'bmrk3'` / `'pain'`) holds 33 participants x 6 heat
levels in a single `fmri_data` object (~198 images). Subject identity and the
within-subject condition (temperature) travel in `additional_info`, which is
exactly what `line_plot_multisubject` needs.

```matlab
% 33 subjects x 6 temperatures, all in one object [voxels x ~198]
[bmrk3, names] = load_image_set('bmrk3', 'noverbose');   % 'bmrk3' == 'pain'

subj  = bmrk3.additional_info.subject_id;     % grouping variable (per image)
temps = bmrk3.additional_info.temperatures;   % within-subject condition (per image)

% Extract one ROI mean per image (same call as Section A).
% Amygdala here; for a more pain-relevant region swap in an insular label,
% e.g. select_atlas_subset(atl, {'Ctx_AVI'}, 'flatten') (anterior insula).
atl = load_atlas('canlab2024', 'noverbose');
roi = select_atlas_subset(atl, {'MTL_AMY'}, 'flatten');   % one bilateral ROI
cl  = extract_roi_averages(bmrk3, roi, 'noverbose');
roi_vec = cl(1).dat;                              % [198 x 1], one mean per image
```

`line_plot_multisubject` requires **cell arrays `X` and `Y`, one cell per
subject**, where each cell is that subject's column vector of paired (x, y)
values. We group the per-image ROI means by subject:

```matlab
usubj = unique(subj);
nsubj = numel(usubj);

X = cell(1, nsubj);   % within-subject x (temperature)
Y = cell(1, nsubj);   % within-subject y (ROI mean)
for i = 1:nsubj
    wh   = subj == usubj(i);     % images belonging to this subject
    X{i} = temps(wh);            % column vector of this subject's temperatures
    Y{i} = roi_vec(wh);          % column vector of this subject's ROI means
end

% One regression line per subject; group t-test on the slopes.
[han, X, Y, slope_stats] = line_plot_multisubject(X, Y, 'center');

% slope_stats.b is [nsubj x 2] (intercept, slope); slope_stats.t/df/p are the
% group t-test on slopes; slope_stats.r_within is the mean within-subject r.
disp(slope_stats.t)
disp(slope_stats.p)
```

Equivalent vector form: instead of building cells, pass the concatenated
vectors plus a `'subjid'` integer vector and let the function split them.

```matlab
% X and Y are plain vectors here (NOT cells) when 'subjid' is given.
[han, X2, Y2, slope_stats2] = line_plot_multisubject(temps, roi_vec, ...
    'subjid', subj, 'center');
```

If you prefer a one-line-per-condition summary instead of per-subject
regression lines, reshape to `[subjects x temperatures]` and use
`lineplot_columns` or `barplot_columns` from Sections A/B:

```matlab
utemp = unique(temps);
Mcond = nan(nsubj, numel(utemp));
for i = 1:nsubj
    for j = 1:numel(utemp)
        Mcond(i, j) = mean(roi_vec(subj == usubj(i) & temps == utemp(j)));
    end
end
lineplot_columns(Mcond, 'color', 'r', 'within');   % mean ROI response per temp level
```

---

## Validation

The unit test
`CanlabCore/Unit_tests/data_extraction/canlab_test_extract_roi_methods.m`
checks that the different extraction paths above (`extract_roi_averages`,
`apply_mask`, `region/extract_data`, and `apply_parcellation`) agree on the same
data to within single-precision tolerance, recover analytically known synthetic
values, and stay highly correlated after on-the-fly resampling. Run it with:

```matlab
results = runtests('canlab_test_extract_roi_methods');
```

or run the whole suite (auto-discovers `canlab_test_*.m`) with
`canlab_run_all_tests`.

---

## Quick reference

| Goal | Method | Input | Output |
|---|---|---|---|
| One/few ROI means per image | `extract_roi_averages(obj, mask)` | `fmri_data` + mask/region/atlas/filename | `region` array; `cl(k).dat` = `[images x 1]` |
| All parcel means | `apply_parcellation(obj, atlas)` | `fmri_data` + `atlas` | `[images x parcels]` matrix |
| All parcel means (atlas method) | `extract_data(atlas, obj)` | `atlas` + `fmri_data` | `[images x parcels]` matrix |
| Bars/violins per condition | `barplot_columns(M)` | `[subjects x conditions]` | handles + stats `table` |
| One mean line across levels | `lineplot_columns(M)` | `[subjects x conditions]` | struct `out` (`.m`, `.ste`, `.CI95`) |
| One line per subject | `line_plot_multisubject(X, Y)` | cell arrays, one cell per subject | handles + `slope_stats` |
