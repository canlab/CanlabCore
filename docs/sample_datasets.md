# Sample datasets and signature patterns in CanlabCore

CanlabCore ships with a small collection of test datasets and exposes a much
larger registry of group-level images, parcellations, meta-analytic maps, and
multivariate "signature" patterns through the `load_image_set` function.
Together, these are the canonical inputs for tutorials, walkthroughs, and the
example scripts in [CANlab_help_examples](https://github.com/canlab/CANlab_help_examples).

Most of the *built-in* test data lives in `CanlabCore/Sample_datasets/` and
ships with the toolbox itself. The remaining datasets and patterns are loaded
on demand from the CANlab `Neuroimaging_Pattern_Masks` repository, the
`MasksPrivate` repository (lab-internal weight maps), the
[canlab_single_trials](https://github.com/canlab/canlab_single_trials)
repository (single-trial pain datasets), or in some cases by automatic
download from Neurovault.

The central entry point is:

```matlab
[image_obj, networknames, imagenames] = load_image_set('keyword');
```

`load_image_set` resolves a keyword (e.g. `'emotionreg'`, `'nps'`, `'kragel18'`,
`'bucknerlab'`) to the matching image set, returns an `fmri_data` object with
the maps loaded, a cell array of formatted network names suitable for plot
labels, and a cell array of full image filenames. If the input is not a
recognized keyword, it is treated as a list of NIfTI filenames and loaded
directly.

For a longer narrative walkthrough see the original CANlab documentation page
at <https://canlab.github.io/_pages/canlab_help_2c_loading_datasets/canlab_help_2c_loading_datasets.html>.
The keyword inventory has grown since that page was written, so the table
below reflects the *current* state of the registry from the live
`load_image_set.m` source. A few datasets — particularly newer signatures
such as `pifonem`, `transcriptomic_gradients`, and the Hansen22 PET maps —
are not described in the legacy help page.

## Built-in `Sample_datasets/` files

The following test data ships with the CanlabCore repository under
`CanlabCore/Sample_datasets/`:

| File / folder | Description | Citation |
|---|---|---|
| `Atlas_2012_REMI_behavioral_data.mat` | Reappraisal-of-emotional-images (REMI) behavioral data accompanying Atlas et al. 2012. | [Atlas, L. Y., Whittington, R. A., Lindquist, M. A., Wielgosz, J., Sonty, N., & Wager, T. D. (2012). Dissociable influences of opiates and expectations on pain. J. Neurosci., 32(23), 8053-8064.](https://doi.org/10.1523/JNEUROSCI.0383-12.2012) <!-- TODO: confirm REMI citation matches this paper --> |
| `emotion regulation_pAgF_z_FDR_0.01_8_14_2015.nii` | Neurosynth-derived "emotion regulation" reverse-inference map (z, FDR 0.01), generated 14-Aug-2015. | [Yarkoni, T., Poldrack, R. A., Nichols, T. E., Van Essen, D. C., & Wager, T. D. (2011). Large-scale automated synthesis of human functional neuroimaging data. Nature Methods, 8(8), 665-670.](https://doi.org/10.1038/nmeth.1635) |
| `Jepma_IE2_single_trial_canlab_dataset.mat` | Single-trial dataset from Jepma IE2 (instructed expectancy 2). Stored as a `canlab_dataset`. | [Jepma, M., Koban, L., van Doorn, J., Jones, M., & Wager, T. D. (2018). Behavioural and neural evidence for self-reinforcing expectancy effects on pain. Nature Human Behaviour, 2(11), 838-855.](https://doi.org/10.1038/s41562-018-0455-8) <!-- TODO: confirm IE2 maps to this paper --> |
| `Pinel_sample_fMRI_time_series/` | Single-subject Pinel localizer raw fMRI timeseries with events file and confounds. | [Pinel, P., Thirion, B., Meriaux, S., Jobert, A., Serres, J., Le Bihan, D., Poline, J.-B., & Dehaene, S. (2007). Fast reproducible identification and large-scale databasing of individual functional cognitive networks. BMC Neuroscience, 8, 91.](https://doi.org/10.1186/1471-2202-8-91) |
| `Wager_et_al_2008_Neuron_EmotionReg/` | 30 first-level contrast images for [reappraise negative vs. look negative], plus metadata. Loaded by `load_image_set('emotionreg')`. | [Wager, T. D., Davidson, M. L., Hughes, B. L., Lindquist, M. A., & Ochsner, K. N. (2008). Prefrontal-subcortical pathways mediating successful emotion regulation. Neuron, 59(6), 1037-1050.](https://doi.org/10.1016/j.neuron.2008.09.006) |
| `Woo_2015_PlosBio_BMRK3_pain_6levels/` | BMRK3 pain stimulation 6-levels dataset: 33 participants, brain responses to six levels of heat. Loaded (when full version is on path) by `load_image_set('bmrk3')` or `'pain'`. | [Woo, C.-W., Roy, M., Buhle, J. T., & Wager, T. D. (2015). Distinct brain systems mediate the effects of nociceptive input and self-regulation on pain. PLoS Biology, 13(1), e1002036.](https://doi.org/10.1371/journal.pbio.1002036) |

## `load_image_set` keyword registry

The registry is organized into four families, mirroring the `:Available
Keywords:` block in the `load_image_set.m` header. Type the keyword (or any
listed alias) as the first argument to `load_image_set`. Keywords are matched
case-insensitively. Some entries require additional CANlab repositories or
private data on the MATLAB path (`MasksPrivate`,
`Neuroimaging_Pattern_Masks`, `canlab_single_trials`); a couple are
auto-downloaded from Neurovault.

Special meta-keywords:

- `'list'` — returns a table of registered signatures (as the first output) and prints it.
- `'all'` — loads every signature in the registry (large object).

### Sample test datasets — one image per subject

| Keyword(s) | Description | Type | Citation |
|---|---|---|---|
| `emotionreg`, `emotionregulation` | N=30 contrast images for [reappraise negative vs. look negative]. | group-level images | [Wager, T. D., Davidson, M. L., Hughes, B. L., Lindquist, M. A., & Ochsner, K. N. (2008). Prefrontal-subcortical pathways mediating successful emotion regulation. Neuron, 59(6), 1037-1050.](https://doi.org/10.1016/j.neuron.2008.09.006) |
| `bmrk3`, `pain` | 33 participants x 6 levels of noxious heat (BMRK3). Requires `bmrk3_6levels_pain_dataset.mat` on path. | group-level images | [Woo, C.-W., Roy, M., Buhle, J. T., & Wager, T. D. (2015). Distinct brain systems mediate the effects of nociceptive input and self-regulation on pain. PLoS Biology, 13(1), e1002036.](https://doi.org/10.1371/journal.pbio.1002036) |
| `kragel270`, `kragel18_alldata`, `kragel2018_alldata`, `kragel18_testdata` | 270 single-subject maps from Kragel 2018 (combined across pain, cognitive control, and negative emotion studies). Auto-downloaded from Neurovault if missing. | group-level images | [Kragel, P. A., Kano, M., Van Oudenhove, L., Ly, H. G., Dupont, P., Rubio, A., Delon-Martin, C., Bonaz, B. L., Manuck, S. B., Gianaros, P. J., Ceko, M., Reynolds Losin, E. A., Woo, C.-W., Nichols, T. E., & Wager, T. D. (2018). Generalizable representations of pain, cognitive control, and negative emotion in medial frontal cortex. Nature Neuroscience, 21(2), 283-289.](https://doi.org/10.1038/s41593-017-0051-7) |

### Sample test datasets — one image per trial (single-trial datasets)

These keywords dispatch to `load_<keyword>.m` in the
[canlab_single_trials](https://github.com/canlab/canlab_single_trials)
repository (compiled and maintained by Bogdan Petre). Each loads as an
`fmri_data` object whose `.metadata_table` field stores per-trial
information. The repository must be on your MATLAB path. The keyword
`all_single_trials` loads all of them.

| Keyword(s) | Description | Type | Citation |
|---|---|---|---|
| `nsf` | NSF heat-pain study (early Wager-lab pain dataset). | single-trial timeseries | [Wager, T. D., Atlas, L. Y., Lindquist, M. A., Roy, M., Woo, C.-W., & Kross, E. (2013). An fMRI-based neurologic signature of physical pain. New England Journal of Medicine, 368(15), 1388-1397.](https://doi.org/10.1056/NEJMoa1204471) <!-- TODO: confirm primary publication for NSF single-trial dataset --> |
| `bmrk3pain` | BMRK3 painful-heat trials (single-trial version of the BMRK3 dataset). | single-trial timeseries | [Woo, C.-W., Roy, M., Buhle, J. T., & Wager, T. D. (2015). Distinct brain systems mediate the effects of nociceptive input and self-regulation on pain. PLoS Biology, 13(1), e1002036.](https://doi.org/10.1371/journal.pbio.1002036) |
| `bmrk3warm` | BMRK3 warm (non-painful) trials. | single-trial timeseries | [Woo, C.-W., Roy, M., Buhle, J. T., & Wager, T. D. (2015). Distinct brain systems mediate the effects of nociceptive input and self-regulation on pain. PLoS Biology, 13(1), e1002036.](https://doi.org/10.1371/journal.pbio.1002036) |
| `bmrk4` | BMRK4 capsaicin/heat dataset; basis for the VPS. | single-trial timeseries | [Krishnan, A., Woo, C.-W., Chang, L. J., Ruzic, L., Gu, X., Lopez-Sola, M., Jackson, P. L., Pujol, J., Fan, J., & Wager, T. D. (2016). Somatic and vicarious pain are represented by dissociable multivariate brain patterns. eLife, 5, e15166.](https://doi.org/10.7554/eLife.15166) <!-- TODO: confirm BMRK4 mapping --> |
| `exp` | Expectancy / placebo single-trial pain dataset. | single-trial timeseries | [Atlas, L. Y., Bolger, N., Lindquist, M. A., & Wager, T. D. (2010). Brain mediators of predictive cue effects on perceived pain. Journal of Neuroscience, 30(39), 12964-12977.](https://doi.org/10.1523/JNEUROSCI.0057-10.2010) <!-- TODO: confirm 'exp' dataset citation --> |
| `ie` | Instructed expectancy (IE) pain dataset. | single-trial timeseries | [Jepma, M., & Wager, T. D. (2015). Conceptual conditioning: Mechanisms mediating conditioning effects on pain. Psychological Science, 26(11), 1728-1739.](https://doi.org/10.1177/0956797615597658) <!-- TODO: confirm IE primary citation --> |
| `ie2` | Instructed expectancy 2 (Jepma IE2). | single-trial timeseries | [Jepma, M., Koban, L., van Doorn, J., Jones, M., & Wager, T. D. (2018). Behavioural and neural evidence for self-reinforcing expectancy effects on pain. Nature Human Behaviour, 2(11), 838-855.](https://doi.org/10.1038/s41562-018-0455-8) |
| `ilcp` | ILCP chronic-pain single-trial dataset. | single-trial timeseries | <!-- TODO: find DOI for ILCP single-trial dataset --> |
| `romantic` | Romantic-rejection / partner-feedback pain study (basis for the rejection signature). | single-trial timeseries | [Woo, C.-W., Koban, L., Kross, E., Lindquist, M. A., Banich, M. T., Ruzic, L., Andrews-Hanna, J. R., & Wager, T. D. (2014). Separate neural representations for physical pain and social rejection. Nature Communications, 5, 5380.](https://doi.org/10.1038/ncomms6380) |
| `scebl` | Social / emotional placebo / social-cue pain study (SCEBL). | single-trial timeseries | [Koban, L., & Wager, T. D. (2016). Beyond conformity: Social influences on pain reports and physiology. Emotion, 16(1), 24-32.](https://doi.org/10.1037/emo0000087) <!-- TODO: confirm SCEBL primary citation --> |
| `stephan` | Stephan et al. pain dataset. | single-trial timeseries | <!-- TODO: find DOI for stephan single-trial dataset --> |
| `all_single_trials` | Loads every single-trial study above into one combined object. | single-trial timeseries | (see individual studies) |

### Parcellations and large-scale networks / patterns

| Keyword(s) | Description | Type | Citation |
|---|---|---|---|
| `bucknerlab` | 7-network cortical parcellation (Yeo/Buckner 2011). Cortex only. | parcellation | [Yeo, B. T. T., Krienen, F. M., Sepulcre, J., Sabuncu, M. R., Lashkari, D., Hollinshead, M., Roffman, J. L., Smoller, J. W., Zollei, L., Polimeni, J. R., Fischl, B., Liu, H., & Buckner, R. L. (2011). The organization of the human cerebral cortex estimated by intrinsic functional connectivity. Journal of Neurophysiology, 106(3), 1125-1165.](https://doi.org/10.1152/jn.00338.2011) |
| `bucknerlab_wholebrain` | 7 networks extended to cortex + basal ganglia + cerebellum. | parcellation | [Yeo, B. T. T., et al. (2011). The organization of the human cerebral cortex estimated by intrinsic functional connectivity. Journal of Neurophysiology, 106(3), 1125-1165.](https://doi.org/10.1152/jn.00338.2011) |
| `bucknerlab_wholebrain_plus` | 7-network parcellation + SPM Anatomy Toolbox regions + brainstem. | parcellation | [Yeo, B. T. T., et al. (2011). The organization of the human cerebral cortex estimated by intrinsic functional connectivity. Journal of Neurophysiology, 106(3), 1125-1165.](https://doi.org/10.1152/jn.00338.2011) |
| `allengenetics` | Five maps of human gene expression compiled by Luke Chang from the Allen Human Brain Atlas. | meta-analytic map | [Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen, E. H., Ng, L., Miller, J. A., et al. (2012). An anatomically comprehensive atlas of the adult human brain transcriptome. Nature, 489(7416), 391-399.](https://doi.org/10.1038/nature11405) |
| `bgloops`, `pauli` | 5 basal-ganglia parcels with associated cortical networks (Pauli 2016). | parcellation | [Pauli, W. M., O'Reilly, R. C., Yarkoni, T., & Wager, T. D. (2016). Regional specialization within the human striatum for diverse psychological functions. PNAS, 113(7), 1907-1912.](https://doi.org/10.1073/pnas.1507610113) |
| `bgloops17`, `pauli17` | 17-parcel striatal regions (Pauli 2016). | parcellation | [Pauli, W. M., O'Reilly, R. C., Yarkoni, T., & Wager, T. D. (2016). Regional specialization within the human striatum for diverse psychological functions. PNAS, 113(7), 1907-1912.](https://doi.org/10.1073/pnas.1507610113) |
| `bgloops_cortex`, `pauli_cortex` | Cortical regions most strongly coupled to the Pauli 5-region striatal clusters. | parcellation | [Pauli, W. M., O'Reilly, R. C., Yarkoni, T., & Wager, T. D. (2016). PNAS, 113(7), 1907-1912.](https://doi.org/10.1073/pnas.1507610113) |
| `pauli_subcortical` | Pauli probabilistic subcortical atlas (CIT168). | parcellation | [Pauli, W. M., Nili, A. N., & Tyszka, J. M. (2018). A high-resolution probabilistic in vivo atlas of human subcortical brain nuclei. Scientific Data, 5, 180063.](https://doi.org/10.1038/sdata.2018.63) |
| `pet_nr_map`, `hansen22`, `pet`, `receptorbinding` | 30 PET tracer / neurotransmitter-receptor maps (Hansen 2022). | meta-analytic map | [Hansen, J. Y., Shafiei, G., Markello, R. D., Smart, K., Cox, S. M. L., Norgaard, M., Beliveau, V., Wu, Y., Gallezot, J.-D., Aumont, E., Servaes, S., Scala, S. G., DuBois, J. M., Wainstein, G., Bezgin, G., Funck, T., Schmitz, T. W., Spreng, R. N., Galovic, M., Koepp, M. J., Duncan, J. S., Coles, J. P., Fryer, T. D., Aigbirhio, F. I., McGinnity, C. J., Hammers, A., Soucy, J.-P., Baillet, S., Guimond, S., Hietala, J., Bedard, M.-A., Leyton, M., Kobayashi, E., Rosa-Neto, P., Ganz, M., Knudsen, G. M., Palomero-Gallagher, N., Shine, J. M., Carson, R. E., Tuominen, L., Dagher, A., & Misic, B. (2022). Mapping neurotransmitter systems to the structural and functional organization of the human neocortex. Nature Neuroscience, 25(11), 1569-1581.](https://doi.org/10.1038/s41593-022-01186-3) |
| `emometa`, `emotionmeta`, `2015emotionmeta` | Meta-analytic maps for 5 basic emotion categories (Anger, Disgust, Fear, Happy, Sad). | meta-analytic map | [Wager, T. D., Kang, J., Johnson, T. D., Nichols, T. E., Satpute, A. B., & Barrett, L. F. (2015). A Bayesian model of category-specific emotional brain responses. PLOS Computational Biology, 11(4), e1004066.](https://doi.org/10.1371/journal.pcbi.1004066) |
| `marg`, `transmodal`, `principalgradient` | Margulies et al. 2016 first principal connectivity gradient (unimodal-to-transmodal). MNI152NLin2009cAsym version. | meta-analytic map | [Margulies, D. S., Ghosh, S. S., Goulas, A., Falkiewicz, M., Huntenburg, J. M., Langs, G., Bezgin, G., Eickhoff, S. B., Castellanos, F. X., Petrides, M., Jefferies, E., & Smallwood, J. (2016). Situating the default-mode network along a principal gradient of macroscale cortical organization. PNAS, 113(44), 12574-12579.](https://doi.org/10.1073/pnas.1608282113) |
| `margfsl` | Same as `marg` in MNI152NLin6Asym (FSL) space. | meta-analytic map | [Margulies, D. S., et al. (2016). PNAS, 113(44), 12574-12579.](https://doi.org/10.1073/pnas.1608282113) |
| `transcriptomic_gradients` | Principal transcriptomic gradients of the human brain. | meta-analytic map | [Hawrylycz, M. J., et al. (2012). An anatomically comprehensive atlas of the adult human brain transcriptome. Nature, 489(7416), 391-399.](https://doi.org/10.1038/nature11405); [Vogel, J. W., Alexander-Bloch, A., Wagstyl, K., Bertolero, M. A., Markello, R. D., Pines, A., Sydnor, V. J., Diaz-Papkovich, A., Hansen, J. Y., Evans, A. C., Bernhardt, B., Misic, B., Satterthwaite, T. D., & Seidlitz, J. (2024). Deciphering the functional specialization of whole-brain spatiomolecular gradients in the adult brain. PNAS, 121(25), e2219379121.](https://doi.org/10.1073/pnas.2219379121) |

### Signature patterns and predictive models

Multivariate "signature" patterns are weight maps from peer-reviewed
predictive-modeling studies. Most are loaded by a single keyword and stored
either in the public `Neuroimaging_Pattern_Masks` repository or in the
private `MasksPrivate` repository (lab-internal). Convenience meta-keywords:

- `npsplus` — NPS (NPSpos / NPSneg), SIIPS, PINES, Romantic Rejection, VPS, etc.
- `painsig` — NPS (NPSpos / NPSneg) and SIIPS only.
- `fibromyalgia`, `fibro`, `fm` — NPSp + FM-pain + FM-multisensory.
- `all` — load all registered signatures.

| Keyword(s) | Description | Type | Citation |
|---|---|---|---|
| `nps` | Neurologic Pain Signature (NPS). Predicts physical heat pain. | signature pattern | [Wager, T. D., Atlas, L. Y., Lindquist, M. A., Roy, M., Woo, C.-W., & Kross, E. (2013). An fMRI-based neurologic signature of physical pain. New England Journal of Medicine, 368(15), 1388-1397.](https://doi.org/10.1056/NEJMoa1204471) |
| `vps` | Vicarious Pain Signature (VPS). Predicts pain observed in others. | signature pattern | [Krishnan, A., Woo, C.-W., Chang, L. J., Ruzic, L., Gu, X., Lopez-Sola, M., Jackson, P. L., Pujol, J., Fan, J., & Wager, T. D. (2016). Somatic and vicarious pain are represented by dissociable multivariate brain patterns. eLife, 5, e15166.](https://doi.org/10.7554/eLife.15166) |
| `rejection` | Romantic-rejection signature. | signature pattern | [Woo, C.-W., Koban, L., Kross, E., Lindquist, M. A., Banich, M. T., Ruzic, L., Andrews-Hanna, J. R., & Wager, T. D. (2014). Separate neural representations for physical pain and social rejection. Nature Communications, 5, 5380.](https://doi.org/10.1038/ncomms6380) |
| `siips` | Stimulus-Intensity-Independent Pain Signature. | signature pattern | [Woo, C.-W., Schmidt, L., Krishnan, A., Jepma, M., Roy, M., Lindquist, M. A., Atlas, L. Y., & Wager, T. D. (2017). Quantifying cerebral contributions to pain beyond nociception. Nature Communications, 8, 14211.](https://doi.org/10.1038/ncomms14211) |
| `pines` | Picture-Induced Negative Emotion Signature. | signature pattern | [Chang, L. J., Gianaros, P. J., Manuck, S. B., Krishnan, A., & Wager, T. D. (2015). A sensitive and specific neural signature for picture-induced negative affect. PLOS Biology, 13(6), e1002180.](https://doi.org/10.1371/journal.pbio.1002180) |
| `gsr` | Stress-induced skin-conductance signature. | signature pattern | [Eisenbarth, H., Chang, L. J., & Wager, T. D. (2016). Multivariate brain prediction of heart rate and skin conductance responses to social threat. Journal of Neuroscience, 36(47), 11987-11998.](https://doi.org/10.1523/JNEUROSCI.3672-15.2016) |
| `hr` | Stress-induced heart-rate signature. | signature pattern | [Eisenbarth, H., Chang, L. J., & Wager, T. D. (2016). Multivariate brain prediction of heart rate and skin conductance responses to social threat. Journal of Neuroscience, 36(47), 11987-11998.](https://doi.org/10.1523/JNEUROSCI.3672-15.2016) |
| `multisensory` | Fibromyalgia multisensory pattern. | signature pattern | [Lopez-Sola, M., Woo, C.-W., Pujol, J., Deus, J., Harrison, B. J., Monfort, J., & Wager, T. D. (2017). Towards a neurophysiological signature for fibromyalgia. Pain, 158(1), 34-47.](https://doi.org/10.1097/j.pain.0000000000000707) |
| `fmpain` | Fibromyalgia pain-period pattern. | signature pattern | [Lopez-Sola, M., et al. (2017). Towards a neurophysiological signature for fibromyalgia. Pain, 158(1), 34-47.](https://doi.org/10.1097/j.pain.0000000000000707) |
| `plspain` | PLS pain-related pattern. | signature pattern | [Kragel, P. A., Kano, M., Van Oudenhove, L., Ly, H. G., Dupont, P., Rubio, A., Delon-Martin, C., Bonaz, B. L., Manuck, S. B., Gianaros, P. J., Ceko, M., Reynolds Losin, E. A., Woo, C.-W., Nichols, T. E., & Wager, T. D. (2018). Generalizable representations of pain, cognitive control, and negative emotion in medial frontal cortex. Nature Neuroscience, 21(2), 283-289.](https://doi.org/10.1038/s41593-017-0051-7) |
| `cpdm` | Combined PDM (multivariate-mediation pain pattern). | signature pattern | [Geuter, S., Reynolds Losin, E. A., Roy, M., Atlas, L. Y., Schmidt, L., Krishnan, A., Koban, L., Wager, T. D., & Lindquist, M. A. (2020). Multiple brain networks mediating stimulus-pain relationships in humans. Cerebral Cortex, 30(7), 4204-4219.](https://doi.org/10.1093/cercor/bhaa048) |
| `pain_pdm`, `pdm` | 10 individual PDM maps and a combined PDM weighting. | signature pattern | [Geuter, S., et al. (2020). Multiple brain networks mediating stimulus-pain relationships in humans. Cerebral Cortex, 30(7), 4204-4219.](https://doi.org/10.1093/cercor/bhaa048) |
| `npsplus` | NPS (incl. NPSpos / NPSneg), SIIPS, PINES, Rejection, VPS and more in one object. | signature pattern | (see component signatures) |
| `painsig` | NPS (incl. NPSpos / NPSneg) and SIIPS only. | signature pattern | (see component signatures) |
| `fibromyalgia`, `fibro`, `fm` | NPSp, FM-pain, FM-multisensory bundle. | signature pattern | [Lopez-Sola, M., et al. (2017). Pain, 158(1), 34-47.](https://doi.org/10.1097/j.pain.0000000000000707) |
| `guilt`, `guilt_behavior` | Guilt-behavior SVM pattern. | signature pattern | [Yu, H., Koban, L., Chang, L. J., Wagner, U., Krishnan, A., Vuilleumier, P., Zhou, X., & Wager, T. D. (2020). A generalizable multivariate brain pattern for interpersonal guilt. Cerebral Cortex, 30(6), 3558-3572.](https://doi.org/10.1093/cercor/bhz326) |
| `neurosynth`, `neurosynth_featureset1` | 525 reverse-inference z-score maps from Neurosynth (2013). | meta-analytic map | [Yarkoni, T., Poldrack, R. A., Nichols, T. E., Van Essen, D. C., & Wager, T. D. (2011). Large-scale automated synthesis of human functional neuroimaging data. Nature Methods, 8(8), 665-670.](https://doi.org/10.1038/nmeth.1635) |
| `neurosynth_topics_forwardinference`, `neurosynth_topics_fi` | 54 forward-inference topic maps from Yarkoni & Poldrack (2014) topic-modeling, with ChatGPT-summarized topic labels from Ke et al. 2024. | meta-analytic map | [Poldrack, R. A., Mumford, J. A., Schonberg, T., Kalar, D., Barman, B., & Yarkoni, T. (2012). Discovering relations between mind, brain, and mental disorders using topic mapping. PLOS Computational Biology, 8(10), e1002707.](https://doi.org/10.1371/journal.pcbi.1002707); [Ke, J., Kang, J., Aerts, H. J. W. L., Wager, T. D., & Lindquist, M. A. (2024). A unified, scalable framework for neural population decoding. Nature Neuroscience, 27, 2310-2321.](https://doi.org/10.1038/s41593-024-01807-z) <!-- TODO: confirm Ke 2024 citation matches the topic-label paper --> |
| `neurosynth_topics_reverseinference`, `neurosynth_topics_ri` | 54 reverse-inference topic maps; same source as above. | meta-analytic map | [Poldrack, R. A., et al. (2012). PLOS Computational Biology, 8(10), e1002707.](https://doi.org/10.1371/journal.pcbi.1002707) |
| `pain_cog_emo`, `kragel18` | PLS maps for generalizable representations of pain, cognitive control, and emotion (24 maps). | signature pattern | [Kragel, P. A., Kano, M., Van Oudenhove, L., Ly, H. G., Dupont, P., Rubio, A., Delon-Martin, C., Bonaz, B. L., Manuck, S. B., Gianaros, P. J., Ceko, M., Reynolds Losin, E. A., Woo, C.-W., Nichols, T. E., & Wager, T. D. (2018). Generalizable representations of pain, cognitive control, and negative emotion in medial frontal cortex. Nature Neuroscience, 21(2), 283-289.](https://doi.org/10.1038/s41593-017-0051-7) |
| `kragelemotion` | 7 emotion-predictive models (Kragel & LaBar 2015). | signature pattern | [Kragel, P. A., & LaBar, K. S. (2015). Multivariate neural biomarkers of emotional states are categorically distinct. Social Cognitive and Affective Neuroscience, 10(11), 1437-1448.](https://doi.org/10.1093/scan/nsv032) |
| `kragelschemas` | 20 visual emotion-schema patterns. | signature pattern | [Kragel, P. A., Reddan, M. C., LaBar, K. S., & Wager, T. D. (2019). Emotion schemas are embedded in the human visual system. Science Advances, 5(7), eaaw4358.](https://doi.org/10.1126/sciadv.aaw4358) |
| `reddanCSplus`, `threat` | CS+ vs. CS- threat-conditioning classifier. | signature pattern | [Reddan, M. C., Wager, T. D., & Schiller, D. (2018). Attenuating neural threat expression with imagination. Neuron, 100(4), 994-1005.e4.](https://doi.org/10.1016/j.neuron.2018.10.047) |
| `zhouvps` | Generalized vicarious-pain signature (Zhou 2020). | signature pattern | [Zhou, F., Li, J., Zhao, W., Xu, L., Zheng, X., Fu, M., Yao, S., Kendrick, K. M., Wager, T. D., & Becker, B. (2020). Empathic pain evoked by sensory and emotional-communicative cues share common and process-specific neural representations. eLife, 9, e56929.](https://doi.org/10.7554/eLife.56929) |
| `multiaversive`, `mpa2` | Multiple Predictive patterns for Aversive experience (MPA2): General, Mechanical, Sounds, Thermal, Visual aversive. | signature pattern | [Ceko, M., Kragel, P. A., Woo, C.-W., Lopez-Sola, M., & Wager, T. D. (2022). Common and stimulus-type-specific brain representations of negative affect. Nature Neuroscience, 25(6), 760-770.](https://doi.org/10.1038/s41593-022-01082-w) |
| `stroop` | Stroop-demand SVM pattern. | signature pattern | [Silvestrini, N., Chen, J.-I., Piche, M., Roy, M., Vachon-Presseau, E., Woo, C.-W., Wager, T. D., & Rainville, P. (2020). Distinct fMRI patterns colocalized in the cingulate cortex underlie the after-effects of cognitive control on pain. NeuroImage, 217, 116898.](https://doi.org/10.1016/j.neuroimage.2020.116898) <!-- TODO: confirm Silvestrini 2020 citation is the Stroop-demand SVM paper --> |
| `ncs` | Neurobiological Craving Signature: combined drug + food, drug-only, food-only weight maps. | signature pattern | [Koban, L., Wager, T. D., & Kober, H. (2023). A neuromarker for drug and food craving distinguishes drug users from non-users. Nature Neuroscience, 26(2), 316-325.](https://doi.org/10.1038/s41593-022-01228-w) |
| `pifonem` | Picture-Induced Fear of Neck Movement (PiFoneM); predicts fear of neck movement in acute and chronic whiplash. | signature pattern | Murillo, C., et al. (2026). PiFoneM: a brain pattern for fear of neck movement in whiplash. Journal of Pain. <!-- TODO: find DOI for Murillo .. Ashar 2026 J Pain PiFoneM paper --> |
| `vifs` | Vicarious / instrument fear pattern (registered via the signature table; see `'list'` for current entry). | signature pattern | <!-- TODO: confirm vifs citation; entry is in the signature-table fallback path --> |
| `list` | Print and return the live signature registry table. | (registry helper) | — |
| `all` | Load every registered signature in one object. | signature pattern | (see component signatures) |

## Worked examples

These examples are adapted from the `:Examples:` block in `load_image_set.m`.

### Loading the NPS plus several other signatures by name

```matlab
imagenames = {'weights_NSF_grouppred_cvpcr.img' ...   % NPS
              'Rating_Weights_LOSO_2.nii'  ...        % PINES
              'dpsp_rejection_vs_others_weights_final.nii' ... % rejection
              'bmrk4_VPS_unthresholded.nii'};         % VPS

[obj, netnames, imgnames] = load_image_set(imagenames);

% Equivalent (and richer) one-liner using the npsplus keyword:
[obj, netnames, imgnames] = load_image_set('npsplus');
```

### Applying the Kragel 2018 PLS signatures to the emotion-regulation dataset

```matlab
% Load PLS signatures from Kragel et al. 2018
[obj, names] = load_image_set('pain_cog_emo');
bpls_wholebrain   = get_wh_image(obj, [8 16 24]);
names_wholebrain  = names([8 16 24]);
bpls_subregions   = get_wh_image(obj, [1:6 9:14 17:22]);
names_subregions  = names([1:6 9:14 17:22]);

% Load test data: emotion regulation contrasts (Wager et al. 2008)
test_data_obj = load_image_set('emotionreg');

% Compare the test data to each Kragel pattern
create_figure('Kragel Pain-Cog-Emo maps', 1, 2);
stats = image_similarity_plot(test_data_obj, 'average', 'mapset', ...
    bpls_wholebrain, 'networknames', names_wholebrain, 'nofigure');
subplot(1, 2, 2)
stats = image_similarity_plot(test_data_obj, 'average', 'mapset', ...
    bpls_subregions, 'networknames', names_subregions, 'nofigure');
```

### Browsing the live signature registry

```matlab
% Print and return the registry as a MATLAB table
sig_table = load_image_set('list');

% Load every registered signature into one fmri_data object
[obj, names] = load_image_set('all');
```

## See also

- The original CANlab loading-datasets walkthrough: <https://canlab.github.io/_pages/canlab_help_2c_loading_datasets/canlab_help_2c_loading_datasets.html>
- [`atlases_regions_and_patterns.md`](atlases_regions_and_patterns.md) — the parallel registry on the *atlas / region* side, loaded via `load_atlas`.
- [`Object_methods.md`](Object_methods.md) — index of object classes (`fmri_data`, `atlas`, `region`, etc.) and their methods.
- [`fmri_data_methods.md`](fmri_data_methods.md) — functional index of `@fmri_data` and `@image_vector` methods organized by area.
- [`toolbox_folders.md`](toolbox_folders.md) — what lives in each subfolder of `CanlabCore/`.
