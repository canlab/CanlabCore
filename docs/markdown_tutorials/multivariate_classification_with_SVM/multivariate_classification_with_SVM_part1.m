%[text] # Multivariate classification with SVM — Part 1
%[text] A walkthrough of basic CANlab/SVM workflows for decoding fMRI condition maps from the **DPSP** romantic-rejection study (Woo et al., *Nat. Commun.* 2014). We train and cross-validate a whole-brain linear SVM on Hot vs. Warm, read out accuracy, AUC, Cohen's *d*, and a confusion matrix, after first inspecting the parallel univariate contrasts.
%[text] ## About the study
%[text] Sixty participants who had recently experienced an unwanted romantic breakup were scanned in two tasks. In the **somatic-pain task** they received painful **Hot** or non-painful **Warm** thermal stimuli; in the **social-rejection task** they viewed photos of their **Ex-partner** (Rejector) or a close **Friend** while recalling breakup-related or positive shared experiences. The four single-subject contrast maps (Hot, Warm, Rejector, Friend) let us pose two parallel two-class classification problems and, by re-using the trained classifiers across tasks, ask whether whole-brain or local multivariate patterns generalize between physical pain and social rejection. The paper's headline finding is that both whole-brain patterns and local patterns within pain-processing regions (dACC, aINS, dpINS, S2) are *separately modifiable*: each can be decoded reliably, but neither generalizes to the other.
%[text] - Paper (open access): Woo et al., 2014, *Nat. Commun.* 5:5380, doi:10.1038/ncomms6380
%[text] - Local pattern weights: |Neuroimaging_Pattern_Masks/Multivariate_signature_patterns/2015_Woo_NatureComms_Rejection| \
%%
%[text] ## 1. Load the data
%[text] The sample data live in |CanlabCore/Sample_datasets/DPSP_pain_rejection_participant_maps| as four |fmri_data| objects — one per condition. Each object holds first-level **condition maps** at the single-subject level, defined as |task − baseline| where baseline is the *implicit* baseline captured by the first-level GLM's intercept. Trials are jittered, so the implicit baseline is estimable. The study used an 18-s visuospatial control task instead of true rest, to avoid rumination contamination.
canlabcore_dir = fileparts(fileparts(which('fmri_data')));
sample_dir = fullfile(canlabcore_dir, 'Sample_datasets', 'DPSP_pain_rejection_participant_maps');
load(fullfile(sample_dir, 'DPSP_single_subject_images_hot.mat'));
load(fullfile(sample_dir, 'DPSP_single_subject_images_warm.mat'));
load(fullfile(sample_dir, 'DPSP_single_subject_images_rejector.mat'));
load(fullfile(sample_dir, 'DPSP_single_subject_images_friend.mat'));
%[text] All four objects have subjects in the same order, so within-person operations like Hot − Warm reduce to subtraction on |.dat|. Each object also carries a |metadata_table| field — a MATLAB |table| with one row per image, with columns |subj_id| (participant ID), |Condition| (label — handy when objects are concatenated for SVM), |Orig_partition| (|xval| or |test|, splitting subjects into a cross-validation set and a final-holdout set), and the original image filename and path.
head(single_subject_images_hot.metadata_table)
%[text] The |Orig_partition| column splits subjects into |xval| (used here for cross-validation, 41 of 59) and |test| (a final-holdout group, 18 of 59 participants, to be evaluated once on the model we ship in a later part).
%%
%[text] ## 1a. Apply a gray-matter mask (optional)
%[text] Restricting to gray matter is a reasonable preprocessing step before pattern analyses — it cuts the feature count and removes voxels with little expected task signal. Non-gray compartments (white matter, CSF) are still useful for denoising and diagnostics (e.g., aCompCor, motion QC), so masking is a *choice*, not an obligation.
gm_mask = fmri_data(which('gray_matter_mask.nii'));
single_subject_images_hot      = apply_mask(single_subject_images_hot,      gm_mask);
single_subject_images_warm     = apply_mask(single_subject_images_warm,     gm_mask);
single_subject_images_rejector = apply_mask(single_subject_images_rejector, gm_mask);
single_subject_images_friend   = apply_mask(single_subject_images_friend,   gm_mask);
n_in_mask_voxels = size(single_subject_images_hot.dat, 1)
%%
%[text] ## 2. Group analyses on within-person contrasts
%[text] Before any classifier, look at the **univariate** picture. We form two within-person contrasts — Hot − Warm and Rejector − Friend — by subject-wise subtraction with |image_vector.image_math|, then run one-sample group t-tests and threshold at *p* < .005 uncorrected, *k* ≥ 10 voxels.
hot_vs_warm        = image_math(single_subject_images_hot,      single_subject_images_warm,   'minus');
rejector_vs_friend = image_math(single_subject_images_rejector, single_subject_images_friend, 'minus');
t_hw = ttest(hot_vs_warm);
t_rf = ttest(rejector_vs_friend);
t_hw = threshold(t_hw, .005, 'unc', 'k', 10);
t_rf = threshold(t_rf, .005, 'unc', 'k', 10);
%%
%[text] ### Hot − Warm montage
%[text] Robust bilateral activation in dorsal posterior insula / S2, mid-insula, thalamus, anterior cingulate, and somatosensory regions — the canonical heat-evoked pain network. The |create_figure| call opens a fresh, named window so the |montage| does not overplot the previous section's figure.
create_figure('Hot − Warm montage'); axis off
montage(t_hw);
%%
%[text] ### Hot − Warm cortical surfaces (foursurfaces_hcp)
create_figure('Hot − Warm surfaces'); axis off
surface(t_hw, 'foursurfaces_hcp');
%%
%[text] ### Rejector − Friend montage
%[text] Sparser at this threshold and dominated by medial-prefrontal, posterior cingulate / precuneus, and right temporoparietal regions — areas commonly engaged by mentalizing about others and negative self-referential affect — together with some anterior-cingulate / insular overlap with the pain map. The interesting question, which we turn to next, is whether the *multivariate* patterns underlying these two contrasts are the *same* or merely *anatomically adjacent*.
create_figure('Rejector − Friend montage'); axis off
montage(t_rf);
%%
%[text] ### Rejector − Friend cortical surfaces
create_figure('Rejector − Friend surfaces'); axis off
surface(t_rf, 'foursurfaces_hcp');
%%
%[text] ## 3. Basic SVMs
%[text] **Decoding** is an umbrella for two related supervised problems: **classification** (predicting a categorical label, e.g. Hot vs Warm) and **regression** (predicting a continuous outcome, e.g. pain rating). We'll focus on classification, using |xval_SVM|. For other algorithms or kernels, |fmri_data.predict| exposes a broader menu (|cv_svm|, |cv_lassopcr|, |cv_pcr|, |cv_pls|, …) with a consistent interface.
%[text] A **linear support vector machine** finds a hyperplane $ f(\\mathbf{x}) = \\mathbf{w}^\\top \\mathbf{x} + b $ that separates the two classes (here $ Y = +1 $ for Hot, $ Y = -1 $ for Warm) while *maximising the margin* between the closest points on each side. With soft-margin slack variables $ \\xi_i \\ge 0 $ (allowing some training-set misclassification) the linear SVM solves $ \\min_{\\mathbf{w},b,\\boldsymbol{\\xi}} \\; \\tfrac{1}{2}\\lVert\\mathbf{w}\\rVert^2 + C\\sum_i \\xi_i $ subject to $ y_i(\\mathbf{w}^\\top \\mathbf{x}_i + b) \\ge 1 - \\xi_i $. Predictions are made from the sign of $ f(\\mathbf{x}) $; the signed value $ f(\\mathbf{x}) $ — the **distance from the hyperplane** — is a useful continuous score.
%[text] To classify Hot vs Warm we **concatenate** the two condition objects, set |.Y| to ±1, and build a **grouping vector** of subject IDs so that both maps from the same participant land in the same fold of cross-validation. Mixing within-person observations across train and test sets is the most common source of data leakage in brain decoding. CANlab provides |xval_stratified_holdout_leave_whole_subject_out| for this — it builds k-fold partitions that keep all images from the same |id| together and stratify each fold on class membership. |xval_SVM| calls it under the hood whenever you pass an |id| vector.
hw_obj = cat(single_subject_images_hot, single_subject_images_warm);
hw_obj = remove_empty(hw_obj);
n = size(single_subject_images_hot.dat, 2);
hw_obj.Y = [ones(n,1); -ones(n,1)];
hw_id_strings = [single_subject_images_hot.metadata_table.subj_id; ...
                 single_subject_images_warm.metadata_table.subj_id];
[~, ~, hw_id] = unique(hw_id_strings, 'stable');
%%
%[text] ## Train and cross-validate
%[text] |xval_SVM| expects an [observations × features] matrix |X|, a [observations × 1] vector |Y| of ±1 outcomes, and a grouping vector |id|. With ~195k features per image, |'highdimensional', true| dispatches to |fitclinear| (faster than |fitcsvm| for wide data). For a quick first pass we turn off hyperparameter optimization, repeats, and bootstrap.
X  = double(hw_obj.dat');
Y  = hw_obj.Y;
id = hw_id;
rng(2026)
S_hw = xval_SVM(X, Y, id, ...
    'highdimensional', true, ...
    'nooptimize', 'norepeats', 'nobootstrap', 'noverbose', 'noplot')
%[text] Useful fields on the returned |predictive_model| object: |dist_from_hyperplane_xval| (cross-validated continuous SVM scores), |yfit| (cross-validated predictions), |crossval_accuracy| (% correct, chance = 50%), |classification_d_singleinterval| and |d_within| (effect sizes — see below), |w| (model weights from the final fit on all data), and |trIdx|/|teIdx| (per-fold splits).
%%
%[text] ## ROC plot
%[text] The ROC sweeps the decision threshold across the cross-validated SVM scores and traces sensitivity (true-positive rate for Hot) against 1 − specificity (false-positive rate). The diagonal is chance. With |'threshold', 0| we evaluate accuracy at the SVM's natural decision boundary (the hyperplane itself); |roc_plot| also reports the threshold-free **AUC** and a Gaussian-model effect size $ d_a $.
create_figure('ROC');
ROC = roc_plot(S_hw.dist_from_hyperplane_xval, S_hw.Y > 0, 'threshold', 0);
set(gca, 'FontSize', 16);
%%
%[text] ## Continuous scores and Cohen's *d*
%[text] |dist_from_hyperplane_xval| is a continuous summary of how strongly the trained brain pattern "votes" for class +1 on each held-out image. Two useful effect sizes fall out of these scores. The **single-interval** *d* treats the scores as a continuous response and computes the standardised mean difference between classes. The **within-person** *d* uses the *paired* score difference (|score_Hot − score_Warm|) for each participant and standardises with the within-subject SD, which is typically smaller than the between-subject SD. Both are usually more sensitive than thresholded accuracy because they exploit the continuous information in the scores rather than collapsing each prediction into a 0/1 hit — accuracy can flatline at chance for a classifier that *would* discriminate well with a better threshold, but the *d* values will still register signal.
classification_d_singleinterval = S_hw.classification_d_singleinterval
classification_d_within         = S_hw.d_within
%%
%[text] ## Confusion matrix
%[text] Rows are *true* labels, columns are *predicted* labels, and cells are percentages of each true class (row-normalised). For this run the classifier is slightly more sensitive to Hot than to Warm at this decision threshold. We pull the raw counts via |confusion_matrix(..., 'noplot')| and then build our own |confusionchart| with text labels (Warm, Hot) rather than the underlying ±1 codes.
[rawConf, normConf] = confusion_matrix(S_hw, 'noplot');
fig_cm = create_figure('Confusion matrix'); clf(fig_cm);
cm = confusionchart(fig_cm, rawConf, {'Warm','Hot'}, ...
    'Title', 'Hot vs Warm — cross-validated', ...
    'RowSummary', 'row-normalized', ...
    'ColumnSummary', 'column-normalized', ...
    'FontSize', 14);
%%
%[text] ## Summary
%[text] We loaded and masked the four DPSP condition objects, computed two univariate group contrasts (Hot − Warm, Rejector − Friend) and visualised them at *p* < .005 uncorrected, then trained and cross-validated a whole-brain linear SVM for Hot vs Warm and read out accuracy, AUC, classification *d*, and the confusion matrix.
%[text] In Part 2 we'll add the parallel rejection classifier, run the **cross-classification** test (does the Hot/Warm pattern discriminate Rejector/Friend, and vice versa?), repeat the analysis within specific regions (dACC, aINS, …) using a searchlight or atlas-based parcellation, evaluate the final model on the held-out |test| partition, and expose more of the methods that live on the |predictive_model| object.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
