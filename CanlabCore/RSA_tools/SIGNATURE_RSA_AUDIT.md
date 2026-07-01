# Cross-subject bodysite SVM-signature RSA — audit & corrected pipeline

Audit of the WASABI **AcceptMap** and **DistractMap** "Bodysite Pain SVM
Signature Visualization Reports" Live Scripts, and the corrected, publication-
defensible reimplementation now living in `RSA_tools/`.

The two `.mlx` files are near-identical; everything below applies to both.

---

## 0. The central conceptual issue (read this first)

The scripts ask a bodysite-**specificity** question but compute it on the wrong
statistical object, and this — not a coding bug — is the main reason the results
are puzzling.

There are **two distinct RSAs** being conflated:

| | What it correlates | What it answers |
|---|---|---|
| **Data-level RSA** (the conventional RSA that *worked* in BodyMap) | the **evoked response patterns** for each bodysite | "is leftarm-pain represented differently from leftleg-pain?" → somatotopy |
| **Signature-level RSA** (what these scripts compute) | the learned **SVM weight maps** | "do two within-site *hot-vs-warm decoders* point the same way?" |

The SVM `weight_obj` is a **backward-model discriminative** map trained to
separate hot from warm **within one bodysite**. Correlating bodysite A's
hot-vs-warm decoder with bodysite B's hot-vs-warm decoder is a question about
**generalizability of pain decoding across sites**, *not* about somatotopic
bodysite representation. So:

- It is expected — not a failure — that signature-RSA shows weak/dirty
  somatotopy while data-level RSA shows clean somatotopy. They measure different
  things.
- The signatures are dominated by a **shared nociceptive (hot>warm) component**
  common to every site (this inflates dot-product / pattern-expression
  off-diagonals) plus **noisy backward-model rotation** (correlated voxels get
  arbitrary, subject-specific weights — this deflates weight-map correlations).
  Neither cleanly recovers somatotopy. The scripts themselves drift toward this
  conclusion ("maps likely contain a shared nociceptive component plus a smaller
  site-specific component … pattern expression is dominated by the shared
  component; correlation isolates the spatial differences").
- The **failure to generalize to AcceptMap/DistractMap, even same-subject**, is
  consistent with this: a hot-vs-warm decoder trained on BodyMap need not express
  cleanly on a different paradigm whose pain manipulation, baseline, and noise
  structure differ — and run-level test maps are low-SNR.

**Recommendation.** Keep the signature cross-subject RSA as a *secondary*
analysis ("are the learned decoders reproducible across people?"). Answer the
**bodysite-specificity** question with **data-level RSA on the evoked patterns**
via the existing `compute_rsm` / `rsa_parcelwise` route (see
`examples/rsa_acceptmap_pipeline.m`), and/or run the signature RSA on
**Haufe / structure-coefficient maps** rather than raw weights.

---

## 1. Audit findings (code-level)

### 1a. Tautological / trivial analyses
- **Self-correlations on the diagonal.** Within-subject 8×8 `corr_*{s}` have
  `diag == 1` by construction; `atanh(1) = Inf`. Any "within" defined as the
  diagonal is trivial. ✅ You already flagged this; the corrected pipeline never
  uses the diagonal for inference (`build_crosssubject_signature_rsa` flags it;
  all tests use `RSA.cross_subject_mask`).

### 1b. Non-independence (the dominant statistical error)
- `get_crosssubject_site_effect` pools **every cross-subject pairwise cell** into
  `same_r` / `diff_r`, then `ttest2(same_z, diff_z)`. With 8 subjects × 8 sites
  there are 8·C(8,2)=**224** same-site cells and **1568** different-site cells,
  but only **8 independent subjects**. Each subject's maps appear in ~14 pairs,
  so the cells are strongly dependent. `ttest2` uses df≈1790 instead of ~7 →
  the t (≈4–6) and p (1e-5…1e-7) are **massively anticonservative**.
- Same problem in the sitewise `ttest2(dpIns_z, S1_z)` (t≈7.4) and in the
  "dpIns > S1" comparison.
- **Fix:** collapse to **one observation per subject** before testing
  (`subject_level_crosssubject_effect`, df = nSub−1 = 7), and/or a
  **within-subject label permutation** test
  (`permutation_test_site_specificity`).

### 1c. The "extend to whole brain" subject loop is degenerate
```matlab
for s = 1:nSub
    idx = RSA_s1_corr.subject_idx ~= s;     % computed but never used
    s1_sub(s) = mean(cellfun(@mean, sim_s1)); % <-- does NOT depend on s
    dp_sub(s) = mean(cellfun(@mean, sim_dp));
end
[~, p, ~, stats] = ttest(dp_sub, s1_sub)
```
`s1_sub`/`dp_sub` are **constant across `s`** (the body ignores `s`), so this
"subject-level" t-test is meaningless. `subject_level_crosssubject_effect`
replaces it with a real per-anchor aggregation.

### 1d. Fisher-z applied indiscriminately
`get_crosssubject_site_effect` does `atanh(M)` regardless of metric. For the
**dot-product / pattern-expression** RSA, `M` is unbounded (values ≫ 1), so
`atanh` is complex/nonsense. **Fix:** `rsa_fisher_z` is only applied to the
correlation metric (clamped to (−1,1)); the helpers return empty z and warn for
other metrics.

### 1e. `apply_mask` argument-order asymmetry (latent, not currently biting you)
`apply_mask(map1, map2, 'pattern_expression', metric)` →
`canlab_pattern_similarity(map1.dat, map2.dat, …)` where `badvals` (excluded
voxels) is defined **from the first argument only**. If two maps had *different*
voxel support, `M(i,j) ≠ M(j,i)`. It is symmetric in your data **only because
all maps within one ROI share the ROI mask**. This is fragile. **Fix:** the
builder stacks all maps on a **common voxel support** and computes
`corr`/cosine/dot in one shot, so symmetry is guaranteed and order-independent.
(Also faster: one `corr(W)` vs an O(K²) `apply_mask` loop.)

### 1f. Metric token confusion
`apply_mask(a, b, 'pattern_expression', 'pattern_expression')` — the second
`'pattern_expression'` is **not** a recognized metric in
`canlab_pattern_similarity`, so it silently falls back to **dot product**. Works,
but is accidental. The corrected builder takes an explicit
`'correlation'|'cosine'|'dotproduct'`.

### 1g. Index fragility
`expr_matrix{s} = svm_stats_cell{s+1}{b}{1}` for `s==8` (an ad-hoc subject-skip)
is a latent off-by-one hazard; the corrected code drives everything off the
explicit `subjects`/`bodysite_names` lists.

---

## Answers to the specific audit questions

- **Fisher z — needed where?** Only for the **correlation** metric, and only when
  averaging or testing. Compute stats in z, report descriptives in r via `tanh`.
  Never `atanh` dot-product/cosine values. Never `atanh` the diagonal.
- **pattern_expression appropriate for RSA?** As the *primary* RSA measure, no —
  it is magnitude-/coverage-driven and dominated by the shared nociceptive
  component. Keep **correlation** primary; pattern-expression/cosine secondary.
- **Correlate SVM weights, Haufe maps, or both?** Report **both**, but lead with
  **Haufe / structure-coefficient** maps for *representational* claims: raw
  backward weights are spatially unreliable (correlated-voxel rotation), which
  artificially deflates cross-subject correlation. Haufe maps are the
  interpretable forward model. Raw weights remain the object for *decoding*
  claims.
- **Full-data weights vs fold-averaged?** For a *single representative signature*,
  the full-data `weight_obj` is the canonical choice and is appropriate for RSA.
  Fold-averaged maps are for **stability/QC** (the report's fold-correlation
  panel), not for inference across folds (folds are not independent).
- **Is full-data training appropriate for representational analysis?** Yes —
  cross-subject comparisons are unaffected by within-subject CV because subjects
  are independent (no cross-subject leakage). Leave-one-session-out CV affects
  *performance* interpretation, not the cross-subject RSA of the weight maps.
- **Are `S.weight_obj.dat` dims compatible across subjects/parcels?** Within one
  ROI, yes (shared mask). Across ROIs/parcels they differ — never mix ROIs in one
  matrix. The builder enforces equal voxel counts and errors with guidance
  otherwise.
- **Missing maps / NaNs?** The builder isolates a common finite, non-zero voxel
  support and flags unextractable maps (`RSA.qc.dropped`) instead of crashing.

---

## 2–6. What the corrected code provides

| Requested | File | Notes |
|---|---|---|
| A `build_crosssubject_signature_rsa` | `build_crosssubject_signature_rsa.m` | common-support matrix, explicit metric, masks, QC |
| B `get_crosssubject_site_effect` | `get_crosssubject_site_effect.m` | **descriptive only**; metric-aware Fisher-z |
| C `get_sitewise_crosssubject_similarity` | `get_sitewise_crosssubject_similarity.m` | descriptive; r and z pools |
| D `collapse_rsa_to_bodysite_matrix` | `collapse_rsa_to_bodysite_matrix.m` | z-space averaging, cross-subject only |
| E `plot_same_vs_different_site` | `plot_same_vs_different_site.m` | **subject-level** error bars + test |
| F `plot_sitewise_crosssubject_similarity` | `plot_sitewise_crosssubject_similarity.m` | subject-level error bars |
| G `subject_level_crosssubject_effect` | `subject_level_crosssubject_effect.m` | leave-one-anchor; paired t, df=nSub−1 |
| H `permutation_test_site_specificity` | `permutation_test_site_specificity.m` | within-subject label permutation |
| (3) Fisher z | `rsa_fisher_z.m` | clamped `atanh` |
| (5) parcelwise | `signature_specificity_parcelwise.m` | per-parcel RSA → brain maps + BH-FDR |
| (6) matrix plots | `plot_signature_rsa_matrix.m`, `add_subject_boundaries.m` | full + collapsed |
| pipeline | `examples/rsa_signature_crosssubject_pipeline.m` | end-to-end |

**Recommended inference (item 3).** Use **both** the subject-level paired test
(primary; interpretable effect size = Cohen's dz) and the within-subject
**permutation** test (confirmatory; distribution-free). A mixed-effects model
(`similarity ~ same_site + (1|subject_i) + (1|subject_j)`) with crossed random
effects is the most complete option, but is overkill at N=8 and harder to defend
than the simple subject-level test; we provide the subject-level + permutation
pair as the practical, robust choice.

Validated on synthetic data: with an injected site component the subject-level
test gives df=nSub−1 and the permutation test p<.001; under a no-site null both
are non-significant (no false positive).

---

## 7. Suggested reporting language (cautious)

> Learned per-subject bodysite signatures (SVM weight maps) were compared across
> subjects with representational similarity analysis. To avoid the
> non-independence of pairwise similarities, similarity was Fisher-z transformed
> and averaged to a single same-site and different-site value per subject; a
> paired t-test across the 8 subjects, confirmed by a within-subject label
> permutation test (5,000 permutations), tested same- vs different-bodysite
> similarity. Cross-subject same-bodysite similarity modestly exceeded
> different-bodysite similarity [report dz, t(7), permutation p], indicating that
> the *discriminative pattern learned for pain at a given site* is partially
> reproducible across individuals. Effect sizes were small (Δr on the order of
> 0.0x), so this reflects weak shared structure rather than a clean somatotopic
> code. dpIns showed [greater/comparable] cross-subject consistency than S1
> [paired t(7), p]. Because these are backward-model discriminative weights —
> dominated by a shared nociceptive component and sensitive to correlated-voxel
> rotation — they are not a direct readout of somatotopic representation; the
> somatotopy question is addressed by representational analysis of the evoked
> response patterns (data-level RSA).

Frame the S1-vs-dpIns dissociation as a **hypothesis** supported by a small
effect, not an established result, and only if it survives the subject-level /
permutation tests. Avoid asserting "dpIns encodes shared interoceptive coding,
S1 idiosyncratic somatotopy" on the basis of Δr≈0.015–0.030.

---

## 8. Threats to validity — flags

1. **Wrong object for the question** (§0): weight-RSA ≠ somatotopy. *Highest
   priority.*
2. **Non-independence** of pairwise cells → use subject-level/permutation. *Fixed.*
3. **Raw weights vs Haufe**: raw backward weights deflate cross-subject
   correlation; rerun on Haufe maps before drawing representational conclusions.
4. **Shared nociceptive component** inflates pattern-expression off-diagonals →
   keep correlation primary; consider partialling out the across-site mean map.
5. **Scale/norm confounds** across subjects → correlation (scale-free) primary;
   never compare dot products across subjects without normalization.
6. **Parcel/ROI mask-size differences** change voxel counts and similarity
   variance → never mix ROIs in one matrix; report per-parcel `n` voxels.
7. **Smoothing / registration**: cross-subject similarity is inflated by shared
   smoothing/template structure and deflated by misregistration; both are
   spatially confounded with the result.
8. **Class imbalance / fold structure**: leave-one-session-out CV with few
   sessions makes fold weights noisy; use full-data weight for the signature.
9. **Self-correlations** (atanh(1)) — excluded by construction. *Fixed.*
10. **Over-interpreting small r** (Δr≈0.0x): statistically detectable ≠
    representationally meaningful. Report effect sizes, not just p.
11. **Low-SNR run-level test maps** in Accept/Distract: the generalization
    "failure" may be a measurement-quality issue, not absence of signal. Try
    subject-averaged single-trial maps or first-level betas before concluding.
12. **Leakage**: none across subjects; within-subject CV does not contaminate
    cross-subject comparisons.

— 2026 Michael Sun. GPL v3.
