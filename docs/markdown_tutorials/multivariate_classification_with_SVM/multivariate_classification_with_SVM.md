# Multivariate classification and cross-classification with CANlab tools (SVM)

A walkthrough of basic CANlab/SVM workflows for *decoding* fMRI condition
maps: training and evaluating a linear SVM, reading out cross-validated
scores, building ROC and confusion-matrix summaries, and laying the
groundwork for cross-classification across distinct affective domains.

> **Live Script version:** [`multivariate_classification_with_SVM_part1.mlx`](multivariate_classification_with_SVM_part1.mlx)
> &nbsp;·&nbsp; plain-text source (recommended for version control, opens in the
> Live Editor on R2025a+):
> [`multivariate_classification_with_SVM_part1.m`](multivariate_classification_with_SVM_part1.m)

## About the study (Woo et al., 2014)

This tutorial uses participant-level condition maps from the *DPSP*
romantic-rejection study (Woo et al., *Nat. Commun.* 2014). Sixty
participants who had recently experienced an unwanted romantic breakup
were scanned in two tasks on separate trials. In the **somatic-pain
task** they received painful **Heat (Hot)** or non-painful **Warmth**
on the volar forearm; in the **social-rejection task** they viewed
photographs of their **Ex-partner (Rejector)** or a **close Friend**,
recalling how they felt during the breakup or a positive shared
experience. The four single-subject contrast maps (Hot, Warm, Rejector,
Friend) let us pose two parallel **two-class classification** problems
(Hot vs. Warm; Rejector vs. Friend) and, by reusing the trained
classifiers across tasks, ask the **cross-classification** question that
motivated the paper: *do whole-brain or local multivariate patterns
generalize between physical pain and social rejection?* The headline
finding is that *both* whole-brain patterns and local patterns within
pain-processing regions (dACC, aINS, dpINS, S2) are *separately
modifiable* — each can be reliably decoded, but neither generalizes to
the other, suggesting distinct neural representations for pain and
rejection co-localized in similar gross anatomy.

- Local repo (Neuroimaging_Pattern_Masks):
  [`2015_Woo_NatureComms_Rejection`](https://github.com/canlab/Neuroimaging_Pattern_Masks/tree/master/Multivariate_signature_patterns/2015_Woo_NatureComms_Rejection)
- Paper (open access): [Woo et al., 2014, *Nat. Commun.* — doi:10.1038/ncomms6380](https://doi.org/10.1038/ncomms6380)
  · [local PDF](../../../Neuroimaging_Pattern_Masks/Multivariate_signature_patterns/2015_Woo_NatureComms_Rejection/Woo_2014_NatComms_dpSP_romantic_rejection.pdf)

## 1. Load the data

The sample data live in `CanlabCore/Sample_datasets/DPSP_pain_rejection_participant_maps/`
as four `fmri_data` objects — one per condition. Each object holds
first-level **condition maps** at the single-subject level (averaged
across runs). We define a condition map as **task − baseline**, where
the baseline is the *implicit* baseline captured by the first-level
GLM's intercept. So `single_subject_images_hot` stores a *Hot − Baseline*
map for each participant, and similarly for the other three conditions.
Because trials are jittered, the implicit baseline is estimable and
task − baseline images are stable. In many designs "baseline" effectively
means rest, but this study deliberately used an **18-s visuospatial
control task** instead — rest would have invited rumination about the
breakup, contaminating the comparison conditions.

```matlab
datadir = fullfile(fileparts(which('canlab_toolbox_setup')), ...
    'CanlabCore', 'Sample_datasets', 'DPSP_pain_rejection_participant_maps');

load(fullfile(datadir, 'DPSP_single_subject_images_hot.mat'));
load(fullfile(datadir, 'DPSP_single_subject_images_warm.mat'));
load(fullfile(datadir, 'DPSP_single_subject_images_rejector.mat'));
load(fullfile(datadir, 'DPSP_single_subject_images_friend.mat'));
```

Every object has **subjects in the same order** (one row per participant
in `.dat'`), so within-person operations like `Hot − Warm` reduce to
matrix subtraction on `.dat`. Each object also carries a
`metadata_table` field — a MATLAB `table` with one row per image. For
these data the columns are:

| column | meaning |
| --- | --- |
| `subj_id`         | participant identifier (one per row) |
| `Condition`       | condition label (`'Hot'`, `'Warm'`, `'Rejector'`, `'Friend'`) — handy when objects are concatenated for SVM |
| `Orig_partition`  | `'xval'` or `'test'` — assignment to the cross-validation set vs. final holdout set |
| `image_names`     | original filename (e.g., `con_0008.img`) |
| `orig_full_path_name` | original path on the lab filesystem |

```matlab
>> head(single_subject_images_hot.metadata_table)

    subj_id      Condition    Orig_partition
    _________    _________    ______________
    'dpsp002'    'Hot'        'xval'
    'dpsp003'    'Hot'        'xval'
    'dpsp004'    'Hot'        'xval'
    'dpsp005'    'Hot'        'test'
    'dpsp006'    'Hot'        'xval'
    ...
```

A **final-holdout** partition (here `'test'`, 18/59 participants) is set
aside so the model you ultimately ship can be evaluated *once* on data
that was never used to choose anything — model, threshold, or
hyperparameter. This is the cleanest form of out-of-sample accuracy you
can report, but it costs you participants in training; somewhere around
**n ≥ 60** is a reasonable floor for it to be meaningful. For the rest
of this tutorial we'll use all 59 participants in cross-validation, then
revisit the held-out partition in a later instalment when we add
cross-classification.

## 1a. Apply a gray-matter mask (optional)

Most CANlab analyses run on whole-brain images by default. Restricting
to gray matter is a reasonable preprocessing step before pattern
analyses — it cuts the feature count, removes voxels with little
expected task signal, and makes diagnostic plots easier to read.
That said, **non-gray compartments are useful too**: white-matter and
CSF time courses are diagnostic of motion and physiological artifact,
they are used by many denoising / normalization procedures (e.g., aCompCor),
and including them at the visualization stage can flag unexpected
artifacts (signal from ventricles, edges, etc.). So masking is a
*choice*, not an obligation.

```matlab
gm_mask = fmri_data(which('gray_matter_mask.nii'));
single_subject_images_hot      = apply_mask(single_subject_images_hot,      gm_mask);
single_subject_images_warm     = apply_mask(single_subject_images_warm,     gm_mask);
single_subject_images_rejector = apply_mask(single_subject_images_rejector, gm_mask);
single_subject_images_friend   = apply_mask(single_subject_images_friend,   gm_mask);
```

After masking each object contains ~195k in-mask voxels (the original
328k include non-brain and non-gray voxels).

## 2. Basic group analyses on the contrasts

Before any classifier, look at the **univariate** picture. We form two
within-person contrasts —

- *Hot − Warm* (pain vs. control)
- *Rejector − Friend* (rejection vs. control)

— by subtracting condition maps subject-by-subject. Because every object
has subjects in the same order in `.dat`, this is a one-liner with
`image_vector.image_math`:

```matlab
hot_vs_warm        = image_math(single_subject_images_hot,      single_subject_images_warm,   'minus');
rejector_vs_friend = image_math(single_subject_images_rejector, single_subject_images_friend, 'minus');
```

Now run a one-sample group t-test on each contrast and threshold at
*p* < .005 uncorrected, *k* ≥ 10 voxels:

```matlab
t_hw = ttest(hot_vs_warm);
t_rf = ttest(rejector_vs_friend);

t_hw = threshold(t_hw, .005, 'unc', 'k', 10);
t_rf = threshold(t_rf, .005, 'unc', 'k', 10);
```

> **Why `create_figure` before each `montage` / `surface` call.** CANlab
> montage and surface methods register an `fmridisplay` object against
> the *current* figure. Without an explicit `create_figure` (which opens
> a fresh, named window and clears it), a second call can overplot the
> previous figure instead of opening a new one — fine on the command
> line, but ugly when the same script runs again in a Live Script. The
> `axis off` keeps the placeholder axes from drawing a frame behind the
> brain slices.

### Hot − Warm (pain contrast)

```matlab
create_figure('Hot − Warm montage'); axis off
montage(t_hw);

create_figure('Hot − Warm surfaces'); axis off
surface(t_hw, 'foursurfaces_hcp');
```

![Hot vs Warm montage](pngs/hotvswarm_montage.png)
![Hot vs Warm surfaces](pngs/hotvswarm_foursurfaces_hcp.png)

Robust bilateral activation in dorsal posterior insula / S2, mid-insula,
thalamus, anterior cingulate, and somatosensory regions — the canonical
heat-evoked pain network.

### Rejector − Friend (rejection contrast)

```matlab
create_figure('Rejector − Friend montage'); axis off
montage(t_rf);

create_figure('Rejector − Friend surfaces'); axis off
surface(t_rf, 'foursurfaces_hcp');
```

![Rejector vs Friend montage](pngs/rejvsfriend_montage.png)
![Rejector vs Friend surfaces](pngs/rejvsfriend_foursurfaces_hcp.png)

The rejection contrast is sparser at this threshold and is dominated by
medial-prefrontal, posterior cingulate / precuneus, and right
temporoparietal regions — areas commonly engaged by mentalizing about
others and negative self-referential affect — together with some
anterior-cingulate / insular overlap with the pain map. The interesting
question, which we turn to now, is whether the *multivariate* patterns
underlying these two contrasts are the *same* or merely *anatomically
adjacent*.

## 3. Basic SVMs

**Decoding** is an umbrella for two related supervised problems:
**classification** (predicting a categorical label, e.g. *Hot vs Warm*)
and **regression** (predicting a continuous outcome, e.g. pain rating).
We'll focus on classification, using
[`xval_SVM`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/Statistics_tools/Cross_validated_Regression/xval_SVM.m).
For other algorithms or kernels, `fmri_data.predict` exposes a broader
menu (`cv_svm`, `cv_lassopcr`, `cv_pcr`, `cv_pls`, …) with a consistent
interface.

### What an SVM does (in one paragraph)

A **linear support vector machine** finds a hyperplane

$$ f(\mathbf{x}) = \mathbf{w}^\top \mathbf{x} + b $$

that separates the two classes (here `Y = +1` for Hot, `Y = −1` for
Warm) while *maximising the margin* between the closest points on each
side. With soft-margin slack variables $\xi_i \ge 0$ (allowing some
training-set misclassification) the linear SVM solves

$$ \min_{\mathbf{w},b,\boldsymbol{\xi}} \; \tfrac{1}{2}\lVert\mathbf{w}\rVert^2 + C\sum_i \xi_i
\quad\text{s.t.}\quad y_i\,(\mathbf{w}^\top \mathbf{x}_i + b) \ge 1 - \xi_i,\; \xi_i \ge 0. $$

Predictions are made from the sign of $f(\mathbf{x})$; the signed value
$f(\mathbf{x})$ — the **distance from the hyperplane** — is a useful
*continuous* score (more on this below). The hyperparameter $C$ trades
off margin width vs. training-set errors; the default in `xval_SVM` is
the MATLAB `fitcsvm` / `fitclinear` default, optionally tuned via nested
cross-validation.

### Building the input: stacked Hot + Warm

To classify Hot vs Warm we concatenate the two condition objects, set
`.Y` to `+1` / `-1`, and build a **grouping vector** of subject IDs so
that **both** maps from the same participant land in the same fold of
cross-validation. Mixing within-person observations across train and
test sets is the most common source of data leakage in brain decoding
and yields wildly optimistic accuracy. CANlab provides
[`xval_stratified_holdout_leave_whole_subject_out`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/Statistics_tools/Cross_validated_Regression/xval_stratified_holdout_leave_whole_subject_out.m)
for this — it builds k-fold partitions that (a) keep all images from
the same `id` together and (b) stratify each fold on class membership
(or, for continuous outcomes, on quartiles of `Y`). `xval_SVM` calls
this routine under the hood whenever you pass an `id` vector.

```matlab
% Stack hot and warm
hw_obj = cat(single_subject_images_hot, single_subject_images_warm);
hw_obj = remove_empty(hw_obj);    % drop the all-zero voxels reintroduced by cat()

% Effects-coded class labels
n = size(single_subject_images_hot.dat, 2);
hw_obj.Y = [ones(n,1); -ones(n,1)];

% Grouping vector: integer subject id, same for both maps of each participant
hw_id_strings = [single_subject_images_hot.metadata_table.subj_id; ...
                 single_subject_images_warm.metadata_table.subj_id];
[~, ~, hw_id] = unique(hw_id_strings, 'stable');
```

### Train and cross-validate

`xval_SVM` expects an `[observations × features]` matrix `X`, a
`[observations × 1]` vector `Y` of `±1` outcomes, and an
`[observations × 1]` grouping vector `id`. With ~195k features per
image, we pass `'highdimensional', true` so it dispatches to
`fitclinear` (much faster than `fitcsvm` for wide data). For a quick
first pass we turn off hyperparameter optimization, repeats, and
bootstrap:

```matlab
X  = double(hw_obj.dat');   % images × voxels
Y  = hw_obj.Y;              % +1 / -1
id = hw_id;                 % subject grouping (1..59)

rng(2026);                  % reproducibility
S_hw = xval_SVM(X, Y, id, ...
    'highdimensional', true, ...
    'nooptimize', 'norepeats', 'nobootstrap');
```

`S_hw` is a `predictive_model` object. The fields you'll most often use:

| field | what it is |
| --- | --- |
| `Y` / `yfit`                  | true and **cross-validated predicted** class labels (`±1`) |
| `dist_from_hyperplane_xval`   | cross-validated continuous SVM scores (signed distance to the hyperplane) — use these as your "what the brain said" measure |
| `class_probability_xval`      | cross-validated class-membership probabilities (Platt scaling) |
| `crossval_accuracy`           | single-interval cross-validated accuracy (%) |
| `classification_d_singleinterval`, `d_within` | classification effect sizes (Cohen's *d*) — see below |
| `w`                           | model weights (one per voxel) from the final model fit to *all* data |
| `trIdx` / `teIdx`             | training / test indices per fold |
| `ClassificationModel`         | the underlying `ClassificationLinear` / `ClassificationSVM` object |

For this run the headline numbers are

| metric | value |
| --- | --- |
| cross-validated accuracy             | **77.1 %** (chance = 50 %; *P* < 10⁻⁸) |
| classification *d* (single-interval) | **1.05** |
| classification *d* (within-person)   | **0.96** |

### Reading the ROC plot

```matlab
create_figure('ROC');
ROC = roc_plot(S_hw.dist_from_hyperplane_xval, S_hw.Y > 0, 'threshold', 0);
set(gca, 'FontSize', 16);
```

![ROC: Hot vs Warm](pngs/svm_hotvswarm_roc.png)

The ROC sweeps the decision threshold across the cross-validated SVM
scores and traces sensitivity (true-positive rate for "Hot") against
1 − specificity (false-positive rate). The diagonal is chance. With
`'threshold', 0` we evaluate accuracy at the SVM's natural decision
boundary (the hyperplane itself); `roc_plot` also reports the
threshold-free **AUC** (here 0.80) and a Gaussian-model effect size
*d_a* (here 1.05). At the *a priori* threshold of 0 we get
*sensitivity* = 81 % (95 % CI 70–91 %), *specificity* = 73 %
(62–84 %), and *PPV* = 75 %.

### Continuous scores → Cohen's *d* (often more sensitive than accuracy)

`dist_from_hyperplane_xval` is a continuous summary of how strongly
the trained brain pattern "votes" for class +1 on each held-out image.
Two useful effect sizes fall out of these scores:

- **Single-interval *d*** — treat the scores as a continuous response
  variable and compute the standardised mean difference between classes.
  This is `S_hw.classification_d_singleinterval` (and equals the Gaussian
  *d_a* reported by `roc_plot`).
- **Within-person *d*** — when you have paired observations per
  subject (as here), compute the *paired* score difference
  (`score_Hot − score_Warm`) for each participant and standardise it.
  This is `S_hw.d_within`, and uses the within-subject standard
  deviation, which is typically smaller than the between-subject SD.

Both are usually more sensitive than thresholded accuracy because they
exploit the continuous information in the scores rather than collapsing
each prediction into a 0/1 hit. Accuracy can flatline at chance for a
classifier that *would* discriminate well with a better threshold;
the *d*s will still register signal.

### Confusion matrix

`confusion_matrix(S_hw)` will plot a chart with the underlying ±1 codes
as axis labels. To get nicer Hot/Warm labels (plus row and column
percentage summaries) we pull the raw counts with `'noplot'` and pass
them to MATLAB's `confusionchart` ourselves:

```matlab
[rawConf, normConf] = confusion_matrix(S_hw, 'noplot');
fig_cm = create_figure('Confusion matrix'); clf(fig_cm);
cm = confusionchart(fig_cm, rawConf, {'Warm','Hot'}, ...
    'Title', 'Hot vs Warm — cross-validated', ...
    'RowSummary', 'row-normalized', ...
    'ColumnSummary', 'column-normalized', ...
    'FontSize', 14);
```

The `clf(fig_cm)` removes the placeholder axes that `create_figure`
opened — `confusionchart` can't share an axes container, so without
this the call errors with *"Adding ConfusionMatrixChart to axes is not
supported. Turn hold off."*

![Confusion matrix: Hot vs Warm](pngs/svm_hotvswarm_confusion.png)

Rows are *true* labels, columns are *predicted* labels, and the cells
along each row and column are summarised as row- and column-normalised
percentages. The classifier correctly classifies 48 / 59 Hot maps and
43 / 59 Warm maps, i.e., it is slightly more sensitive to *Hot* than to *Warm*
at this decision threshold. Raw counts are in `rawConf`.

---

## What's next

This is the end of *Part 1*. We have

- loaded and masked the four DPSP condition objects,
- computed two univariate group contrasts (*Hot − Warm*, *Rejector − Friend*) and visualised them at *p* < .005 uncorrected,
- trained and cross-validated a whole-brain linear SVM for *Hot vs Warm*,
- read out accuracy, AUC, classification *d*, and the confusion matrix.

In Part 2 we'll add the parallel rejection classifier, run the
**cross-classification** test (does the Hot/Warm pattern discriminate
Rejector/Friend, and vice versa?), repeat the analysis **within
specific regions** (dACC, aINS, …) using a searchlight or atlas-based
parcellation, evaluate the final model on the held-out `'test'`
partition, and expose more of the methods that live on the
`predictive_model` object.

## References

- Woo C-W, Koban L, Kross E, Lindquist MA, Banich MT, Ruzic L, Andrews-Hanna JR, Wager TD (2014). **Separate neural representations for physical pain and social rejection.** *Nature Communications* 5:5380. [doi:10.1038/ncomms6380](https://doi.org/10.1038/ncomms6380)
- CANlab tutorials and walkthroughs: <https://github.com/canlab/CANlab_help_examples>
- CanlabCore function reference: <https://canlabcore.readthedocs.org/en/latest/>
