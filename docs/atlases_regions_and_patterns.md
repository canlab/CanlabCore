# Atlases, named regions, and signature patterns

CanlabCore consumes a curated collection of brain atlases / parcellations,
named regions of interest, multivariate signature patterns, and
meta-analytic maps. The image files themselves are distributed via the
[canlab/Neuroimaging_Pattern_Masks](https://github.com/canlab/Neuroimaging_Pattern_Masks)
repository, while the keyword registries that resolve a short name to a
specific saved object live inside CanlabCore. Three loader functions are
the main entry points:

- `load_atlas('keyword')` returns an `atlas` object (a labeled
  parcellation with `.dat`, `.probability_maps`, `.labels`,
  `.label_descriptions`, and `.references`). See
  [`atlas_methods.md`](atlas_methods.md) for what you can do with it.
- `load_image_set('keyword')` returns an `fmri_data` object holding
  one or more multivariate signature patterns or meta-analytic maps,
  ready to apply with `apply_mask`, `canlab_pattern_similarity`, or
  `image_similarity_plot`. The full keyword table for `load_image_set`
  (datasets and signatures) lives in
  [`sample_datasets.md`](sample_datasets.md).
- `canlab_load_ROI('region_name')` returns a `region` object for a
  single named region (e.g. `'pag'`, `'lc'`, `'nacc'`), drawn from
  hand-picked sources across atlases and individual papers. Optionally
  also returns a 1 mm `atlas` object covering that region.

For curated descriptions, figures, and usage notes, see the companion
[CANlab Brain Patterns](https://sites.google.com/dartmouth.edu/canlab-brainpatterns/home)
site. The main toolbox documentation is at
[canlab.github.io](https://canlab.github.io).

The `atlas` keywords typically accept suffixes that select the reference
space (`_fmriprep20` for MNI152NLin2009cAsym, `_fsl6` for MNI152NLin6Asym)
and sampling resolution (`_1mm`, `_2mm`); some accept granularity
(`_fine`, `_coarse`). Defaults vary by atlas — see the docstring of
`load_atlas` for the full list.

## Atlases and parcellations

### Whole-brain combined atlases

The `canlab2018`, `canlab2023`, `canlab2024`, and `opencanlab2024`
atlases are dynamically assembled whole-brain parcellations that combine
several published atlases into a single, internally consistent
labeling. Cortex is taken from the Glasser HCP-MMP1 parcellation,
subcortex from CIT168 (and Pauli striatum / Iglesias thalamus in newer
builds), thalamus from Iglesias and/or Morel, cerebellum from SUIT
(Diedrichsen), and brainstem nuclei from Bianciardi (canlab2023/2024,
which therefore require local assembly because the Bianciardi atlas
cannot be redistributed). `opencanlab2024` is a fully open-license
variant that omits Bianciardi (synthesizing approximations of some
brainstem nuclei) and ships pre-assembled. `canlab2018` is the legacy
build, deprecated in favor of `canlab2023`. See the README of the
[Neuroimaging_Pattern_Masks](https://github.com/canlab/Neuroimaging_Pattern_Masks)
repository for the per-region provenance table.

### Cortical parcellations

The **Glasser HCP-MMP1** multimodal parcellation defines 180 areas per
hemisphere from cortical myelin, thickness, function, and connectivity
in the HCP cohort. CanlabCore ships volumetric projections built with
registration fusion to the fmriprep (MNI152NLin2009cAsym) and FSL
(MNI152NLin6Asym) reference spaces.
[Glasser, M.F. et al. (2016). A multi-modal parcellation of human cerebral cortex. Nature, 536, 171-178.](https://doi.org/10.1038/nature18933)

The **Desikan-Killiany** atlas labels 34 cortical gyri per hemisphere
based on sulcal anatomy and is the default Freesurfer cortical labeling.
[Desikan, R.S. et al. (2006). An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest. NeuroImage, 31(3), 968-980.](https://doi.org/10.1016/j.neuroimage.2006.01.021)

The **DKT (Desikan-Killiany-Tourville)** atlas refines the
Desikan-Killiany boundaries to be more anatomically consistent across
labs.
[Klein, A. & Tourville, J. (2012). 101 labeled brain images and a consistent human cortical labeling protocol. Frontiers in Neuroscience, 6, 171.](https://doi.org/10.3389/fnins.2012.00171)

The **Destrieux** atlas (also distributed with Freesurfer) divides the
cortex into ~74 gyral and sulcal parcels per hemisphere, defined by
local curvature.
[Destrieux, C., Fischl, B., Dale, A. & Halgren, E. (2010). Automatic parcellation of human cortical gyri and sulci using standard anatomical nomenclature. NeuroImage, 53(1), 1-15.](https://doi.org/10.1016/j.neuroimage.2010.06.010)

The **Schaefer 2018 / Yeo 17 networks** parcellation provides 100-1000
cortical parcels (CanlabCore ships the 400-parcel and 17-network
versions) optimized for local gradient and global similarity in
resting-state fMRI.
[Schaefer, A. et al. (2018). Local-global parcellation of the human cerebral cortex from intrinsic functional connectivity MRI. Cerebral Cortex, 28(9), 3095-3114.](https://doi.org/10.1093/cercor/bhx179)
The 17-network labeling is from
[Yeo, B.T.T. et al. (2011). The organization of the human cerebral cortex estimated by intrinsic functional connectivity. Journal of Neurophysiology, 106(3), 1125-1165.](https://doi.org/10.1152/jn.00338.2011)

The **Shen 268** parcellation is a functional connectivity-based atlas
covering cortex, subcortex, and cerebellum (CanlabCore ships fmriprep
and FSL projections in addition to the original Colin27 space).
[Shen, X., Tokoglu, F., Papademetris, X. & Constable, R.T. (2013). Groupwise whole-brain parcellation from resting-state fMRI data for network node identification. NeuroImage, 82, 403-415.](https://doi.org/10.1016/j.neuroimage.2013.05.081)

### Subcortical atlases

The **CIT168** subcortical atlas is a probabilistic parcellation of 16
subcortical structures (basal ganglia, amygdala, thalamic and
diencephalic nuclei) derived from high-resolution multi-contrast 7T
templates. CanlabCore distributes v1.1.0 in fmriprep and FSL spaces, plus
a separate amygdala-only sub-atlas (v1.0.3).
[Pauli, W.M., Nili, A.N. & Tyszka, J.M. (2018). A high-resolution probabilistic in vivo atlas of human subcortical brain nuclei. Scientific Data, 5, 180063.](https://doi.org/10.1038/sdata.2018.63)
The CIT168 amygdala parcellation is described in
[Tyszka, J.M. & Pauli, W.M. (2016). In vivo delineation of subdivisions of the human amygdaloid complex in a high-resolution group template. Human Brain Mapping, 37(11), 3979-3998.](https://doi.org/10.1002/hbm.23289)

The **Pauli striatum** atlas decomposes the striatum into five
probabilistic functional zones using diffusion-derived connectivity.
[Pauli, W.M., O'Reilly, R.C., Yarkoni, T. & Wager, T.D. (2016). Regional specialization within the human striatum for diverse psychological functions. PNAS, 113(7), 1907-1912.](https://doi.org/10.1073/pnas.1507610113)

The **Brainnetome** atlas parcellates cortex and subcortex into 246
regions defined by anatomical connectivity from diffusion MRI.
[Fan, L. et al. (2016). The Human Brainnetome Atlas: A new brain atlas based on connectional architecture. Cerebral Cortex, 26(8), 3508-3526.](https://doi.org/10.1093/cercor/bhw157)

### Thalamus

The **Morel / Krauth** thalamic atlas is a histology-based stereotaxic
parcellation of human thalamic nuclei.
[Krauth, A. et al. (2010). A mean three-dimensional atlas of the human thalamus: generation from multiple histological data. NeuroImage, 49(3), 2053-2062.](https://doi.org/10.1016/j.neuroimage.2009.10.042)
The original histological work is
[Morel, A., Magnin, M. & Jeanmonod, D. (1997). Multiarchitectonic and stereotactic atlas of the human thalamus. Journal of Comparative Neurology, 387(4), 588-630.](https://doi.org/10.1002/(SICI)1096-9861(19971103)387:4%3C588::AID-CNE8%3E3.0.CO;2-Z)

The **Iglesias thalamus** atlas is a probabilistic Bayesian segmentation
of 26 thalamic nuclei using ex vivo histology and in vivo MRI; the
CanlabCore build is projected to fmriprep and FSL spaces.
[Iglesias, J.E. et al. (2018). A probabilistic atlas of the human thalamic nuclei combining ex vivo MRI and histology. NeuroImage, 183, 314-326.](https://doi.org/10.1016/j.neuroimage.2018.08.012)

The **Tian 3T subcortex** atlas provides four levels of granularity
(S1-S4, 16-54 regions) of subcortical parcellation from 3T resting-state
fMRI gradients, and is available in both fmriprep and FSL spaces.
[Tian, Y., Margulies, D.S., Breakspear, M. & Zalesky, A. (2020). Topographic organization of the human subcortex unveiled with functional connectivity gradients. Nature Neuroscience, 23(11), 1421-1432.](https://doi.org/10.1038/s41593-020-00711-6)

### Hypothalamus

The **Iglesias / Billot hypothalamus** segmentation uses a CNN to
delineate five hypothalamic subregions from T1 MRI.
[Billot, B. et al. (2020). Automated segmentation of the hypothalamus and associated subunits in brain MRI. NeuroImage, 223, 117287.](https://doi.org/10.1016/j.neuroimage.2020.117287)

### Brainstem and cerebellum

The **SUIT cerebellar atlas (Diedrichsen)** is a high-resolution lobular
parcellation of the cerebellum aligned to a dedicated cerebellar
template.
[Diedrichsen, J., Balsters, J.H., Flavell, J., Cussans, E. & Ramnani, N. (2009). A probabilistic MR atlas of the human cerebellum. NeuroImage, 46(1), 39-46.](https://doi.org/10.1016/j.neuroimage.2009.01.045)
The functional/lobular parcellation used in newer SUIT releases is
[Diedrichsen, J. & Zotow, E. (2015). Surface-based display of volume-averaged cerebellar imaging data. PLoS ONE, 10(7), e0133402.](https://doi.org/10.1371/journal.pone.0133402)

The **Bianciardi mesopontine atlas** provides probabilistic templates for
~30 brainstem nuclei built from 7T multimodal in vivo data. Distribution
is restricted, so CanlabCore assembles a derivative atlas locally on
first use.
[Bianciardi, M. et al. (2018). A probabilistic template of human mesopontine tegmental nuclei from in vivo 7T MRI. NeuroImage, 170, 222-230.](https://doi.org/10.1016/j.neuroimage.2017.04.070)
See also the project page at the
[Brainstem Imaging Lab](https://www.nmr.mgh.harvard.edu/resources/brainstemimagingatlas).

The **Harvard Ascending Arousal Network (AAN) v2.0** atlas labels
brainstem and forebrain nuclei central to arousal and consciousness,
based on histology, immunohistochemistry, and diffusion tractography.
[Edlow, B.L. et al. (2012). Neuroanatomic connectivity of the human ascending arousal system critical to consciousness and its disorders. Journal of Neuropathology and Experimental Neurology, 71(6), 531-546.](https://doi.org/10.1097/NEN.0b013e3182588293)
v2.0 is distributed via the
[Center for Neurotechnology and Neurorecovery AAN page](https://www.martinos.org/resources/aan-atlas).

The **Levinson-Bari Limbic Brainstem Atlas** is a probabilistic, openly
licensed atlas of brainstem limbic structures including VTA, dorsal
raphe, locus coeruleus, NTS, and PAG.
[Levinson, S., Miller, M., Iftekhar, A., Justo, M., Arriola, D., Wei, W., Hazany, S., Avecillas-Chasin, J.M., Kuhn, T.P., Horn, A. et al. (2023). A structural connectivity atlas of limbic brainstem nuclei. Frontiers in Neuroimaging, 1, 1009399.](https://doi.org/10.3389/fnimg.2022.1009399)

The **Kragel 2019 PAG** atlas (`kragel2019pag` keyword) is a 7T-derived
sub-parcellation of the periaqueductal gray into dorsal, lateral, and
ventrolateral columns via unsupervised k-means clustering of voxel
positions.
[Kragel, P.A., Bianciardi, M., Hartley, L., Matthewson, G., Choi, J.-K., Quigley, K.S., Wald, L.L., Wager, T.D., Barrett, L.F., & Satpute, A.B. (2019). Functional involvement of human periaqueductal gray and other midbrain nuclei in cognitive control. Journal of Neuroscience, 39(31), 6180-6189.](https://doi.org/10.1523/JNEUROSCI.2043-18.2019)

### 7T high-resolution

The **Keuken 7T atlas** is a probabilistic atlas of subcortical
structures (STN, SN, RN, GPe, GPi, thalamus) built from 7T quantitative
MRI in 30 healthy adults.
[Keuken, M.C. et al. (2014). Quantifying inter-individual anatomical variability in the subcortex using 7T structural MRI. NeuroImage, 94, 40-46.](https://doi.org/10.1016/j.neuroimage.2014.03.032)

### Functional networks

The **Buckner / Yeo 2011** networks atlas provides 7- and 17-network
cortical parcellations from resting-state fMRI in 1000 subjects.
CanlabCore ships an atlas object that combines the cortical networks
with associated cerebellar (Buckner et al. 2011) and striatal (Choi
et al. 2012) functional networks.
[Yeo, B.T.T. et al. (2011). The organization of the human cerebral cortex estimated by intrinsic functional connectivity. Journal of Neurophysiology, 106(3), 1125-1165.](https://doi.org/10.1152/jn.00338.2011)
[Buckner, R.L. et al. (2011). The organization of the human cerebellum estimated by intrinsic functional connectivity. Journal of Neurophysiology, 106(5), 2322-2345.](https://doi.org/10.1152/jn.00339.2011)
[Choi, E.Y., Yeo, B.T.T. & Buckner, R.L. (2012). The organization of the human striatum estimated by intrinsic functional connectivity. Journal of Neurophysiology, 108(8), 2242-2263.](https://doi.org/10.1152/jn.00270.2012)

### Pain pathways

The **CANlab pain pathways** atlas is a curated combination of
structures involved in nociceptive processing — spinothalamic and
trigeminothalamic targets in thalamus, somatosensory cortex, posterior
and middle insula, dorsal posterior cingulate, amygdala, hypothalamus,
PAG, parabrachial complex, and rostral ventromedial medulla. It is
documented and figured on the
[CANlab Brain Patterns: pain pathways](https://sites.google.com/dartmouth.edu/canlab-brainpatterns/home/pain-pathways)
page. The companion paper for the underlying signature work is
[Wager, T.D. et al. (2013). An fMRI-based neurologic signature of physical pain. New England Journal of Medicine, 368, 1388-1397.](https://doi.org/10.1056/NEJMoa1204471)
A 2024 refresh (`painpathways2024`) re-bases the pathway on canlab2024
parcels and adds finer brainstem coverage; see the
[Neuroimaging_Pattern_Masks](https://github.com/canlab/Neuroimaging_Pattern_Masks)
README.

### Other atlases

The **Faillenot insular atlas** divides the human insula into three
cytoarchitectonic gyri.
[Faillenot, I., Heckemann, R.A., Frot, M. & Hammers, A. (2017). Macroanatomy and 3D probabilistic atlas of the human insula. NeuroImage, 150, 88-98.](https://doi.org/10.1016/j.neuroimage.2017.01.073)

The **Cartmell NAc core/shell atlas** provides probabilistic
delineations of nucleus accumbens core and shell.
[Cartmell, S.C.D. et al. (2019). Multimodal characterization of the human nucleus accumbens. NeuroImage, 198, 137-149.](https://doi.org/10.1016/j.neuroimage.2019.05.019)

The **de la Vega Neurosynth atlas** is a meta-analytically derived
co-activation parcellation of medial frontal cortex (and a related
whole-cortex variant) from automated meta-analysis on Neurosynth.
[de la Vega, A., Chang, L.J., Banich, M.T., Wager, T.D. & Yarkoni, T. (2016). Large-scale meta-analysis of human medial frontal cortex reveals tripartite functional organization. Journal of Neuroscience, 36(24), 6553-6562.](https://doi.org/10.1523/JNEUROSCI.4402-15.2016)
The companion Neurosynth platform is documented at
[neurosynth.org](https://neurosynth.org).

The **Julich-Brain** cytoarchitectonic atlas is a probabilistic
histology-based parcellation covering most of cortex and many
subcortical structures.
[Amunts, K., Mohlberg, H., Bludau, S. & Zilles, K. (2020). Julich-Brain: A 3D probabilistic atlas of the human brain's cytoarchitecture. Science, 369(6506), 988-992.](https://doi.org/10.1126/science.abb4588)
Live updates and download links are at the
[Julich-Brain page on EBRAINS](https://search.kg.ebrains.eu/instances/Project/e39a0407-a98a-480e-9c63-4a2225ddfbe4).

## Named regions (`canlab_load_ROI`)

Named regions are hand-picked single-region masks selected from atlases
or individual papers, returned as `region` objects suitable for display
(`addbrain`, `montage`) or extraction (`extract_roi_averages`). Regions
are grouped here by anatomy.

### Cortex

| Keyword | Description | Source |
|---|---|---|
| `vmpfc` | Ventromedial prefrontal + posterior cingulate, midline | Hand-drawn (Tor Wager) |

### Forebrain (non-basal ganglia)

| Keyword | Description | Source |
|---|---|---|
| `nacc`, `nac` | Nucleus accumbens | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) (CIT168) |
| `amygdala` | Amygdala | [Tyszka & Pauli 2016](https://doi.org/10.1002/hbm.23289) (CIT168) |
| `hipp` | Hippocampus | [Eickhoff et al. 2005](https://doi.org/10.1016/j.neuroimage.2004.12.034) (SPM Anatomy Toolbox via canlab2018) |
| `BST`, `BNST`, `SLEA` | Bed nucleus of stria terminalis / sublenticular extended amygdala | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) |

### Basal ganglia

| Keyword | Description | Source |
|---|---|---|
| `caudate`, `cau` | Caudate nucleus | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) (CIT168) |
| `put`, `putamen` | Putamen | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) (CIT168) |
| `gp` | Globus pallidus (GPe + GPi) | [Keuken et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.03.032) |
| `GPe` | Globus pallidus externa | [Keuken et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.03.032) |
| `GPi` | Globus pallidus interna | [Keuken et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.03.032) |
| `VeP`, `vep`, `vpall` | Ventral pallidum | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) |

### Thalamus, diencephalon, epithalamus

| Keyword | Description | Source |
|---|---|---|
| `thalamus`, `thal` | Thalamus main body | [Krauth et al. 2010](https://doi.org/10.1016/j.neuroimage.2009.10.042) (Morel) |
| `md` | Mediodorsal thalamus | [Krauth et al. 2010](https://doi.org/10.1016/j.neuroimage.2009.10.042) |
| `cm` | Centromedian thalamus | [Krauth et al. 2010](https://doi.org/10.1016/j.neuroimage.2009.10.042) |
| `lgn` | Lateral geniculate nucleus | [Krauth et al. 2010](https://doi.org/10.1016/j.neuroimage.2009.10.042) |
| `mgn` | Medial geniculate nucleus | [Krauth et al. 2010](https://doi.org/10.1016/j.neuroimage.2009.10.042) |
| `VPthal`, `VPL` | Ventral posterior thalamus | [Krauth et al. 2010](https://doi.org/10.1016/j.neuroimage.2009.10.042) |
| `intralaminar_thal` | Intralaminar / midline thalamus | [Krauth et al. 2010](https://doi.org/10.1016/j.neuroimage.2009.10.042) |
| `stn` | Subthalamic nucleus | [Keuken et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.03.032) |
| `habenula`, `HN` | Habenula | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) |
| `mammillary`, `mamm` | Mammillary bodies | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) |
| `hypothalamus`, `hy`, `hythal` | Hypothalamus | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) (CIT168) |

### General brainstem

| Keyword | Description | Source |
|---|---|---|
| `brainstem` | Whole brainstem (cleaned SPM8 tissue probability map) | Hand-cleaned (Tor Wager), SPM8 TPM |
| `midbrain` | Whole midbrain | [Carmack et al. 2004](https://doi.org/10.1118/1.1646232) <!-- TODO: confirm exact Carmack 2004 reference for midbrain mask --> |

### Midbrain

| Keyword | Description | Source |
|---|---|---|
| `pag` | Periaqueductal gray (hand-drawn, aqueduct masked) | Tor Wager 2018; aqueduct from [Keuken et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.03.032) |
| `sc` | Superior colliculus (hand-drawn) | Tor Wager 2018; aqueduct from [Keuken et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.03.032) |
| `ic` | Inferior colliculus (hand-drawn) | Tor Wager 2018; aqueduct from [Keuken et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.03.032) |
| `drn` | Dorsal raphe nucleus | [Beliveau et al. 2015](https://doi.org/10.1016/j.neuroimage.2015.04.065) |
| `mrn` | Median raphe nucleus | [Beliveau et al. 2015](https://doi.org/10.1016/j.neuroimage.2015.04.065) |
| `ncf` | Nucleus cuneiformis | [Zambreanu et al. 2005](https://doi.org/10.1016/j.pain.2005.01.005) |
| `PBP` | Parabrachial pigmented nucleus | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) |
| `sn` | Substantia nigra | [Keuken et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.03.032) |
| `SNc`, `snc` | Substantia nigra compacta | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) |
| `SNr`, `snr` | Substantia nigra reticularis | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) |
| `VTA`, `vta` | Ventral tegmental area | [Pauli et al. 2018](https://doi.org/10.1038/sdata.2018.63) |
| `rn` | Red nucleus | [Keuken et al. 2014](https://doi.org/10.1016/j.neuroimage.2014.03.032) |

### Pons

| Keyword | Description | Source |
|---|---|---|
| `pbn` | Parabrachial complex (current, from canlab2024) | [Bianciardi et al. 2018](https://doi.org/10.1016/j.neuroimage.2017.04.070) via canlab2024 |
| `pbn_old` | Parabrachial complex (legacy) | [Fairhurst et al. 2007](https://doi.org/10.1016/j.pain.2006.09.001) |
| `lc`, `locus_coeruleus` | Locus coeruleus (2-SD probability mask) | [Keren et al. 2009](https://doi.org/10.1016/j.neuroimage.2009.05.071) |
| `nrp_B5` | Nucleus raphe pontis (B5) | [Bär et al. 2016](https://doi.org/10.1016/j.neuroimage.2016.03.071) |

### Medulla

| Keyword | Description | Source |
|---|---|---|
| `rvm`, `rvm_brooks` | Rostral ventromedial medulla | [Brooks et al. 2017](https://doi.org/10.1097/j.pain.0000000000000789) <!-- TODO: confirm year/DOI of Brooks RVM mask --> |
| `rvm_old` | Rostral ventromedial medulla (legacy hand-drawn) | Tor Wager (legacy) |
| `nrm`, `raphe magnus` | Nucleus raphe magnus | [Bär et al. 2016](https://doi.org/10.1016/j.neuroimage.2016.03.071) |
| `medullary_raphe` | Medullary raphe | [Nash et al. 2009](https://doi.org/10.1002/hbm.20805) |
| `ncs_B6_B8` | Nucleus centralis superior (B6/B8) | [Bär et al. 2016](https://doi.org/10.1016/j.neuroimage.2016.03.071) |
| `nuc_ambiguus` | Nucleus ambiguus | [Sclocco et al. 2016](https://doi.org/10.1098/rsta.2015.0189) |
| `dmnx_nts`, `nts` | Dorsal motor nucleus of vagus / nucleus of the solitary tract | [Sclocco et al. 2016](https://doi.org/10.1098/rsta.2015.0189) |
| `spinal_trigeminal` | Spinal trigeminal nucleus | [Nash et al. 2009](https://doi.org/10.1002/hbm.20805) |
| `olive` | Inferior olive (binary mask) | Hand-drawn / legacy <!-- TODO: locate primary source for ROI_inf_olive.img --> |

## Multivariate signature patterns

CanlabCore distributes a curated set of multivariate brain patterns —
spatially-distributed weight maps that, when applied to a new image
(typically with `apply_mask` plus a dot product, or with `canlab_pattern_similarity`
/ `apply_nps` / `apply_all_signatures`), yield a scalar prediction for a
psychological or clinical state. The full keyword table for
`load_image_set` (signatures plus sample fMRI datasets) is in
[`sample_datasets.md`](sample_datasets.md). A few of the most prominent
signatures and their original publications:

- **NPS — Neurologic Pain Signature.** Pattern predictive of acute
  experimental thermal pain, with strong sensitivity and specificity to
  noxious heat.
  [Wager, T.D. et al. (2013). An fMRI-based neurologic signature of physical pain. New England Journal of Medicine, 368, 1388-1397.](https://doi.org/10.1056/NEJMoa1204471)
- **SIIPS — Stimulus Intensity Independent Pain Signature.** A pattern
  for pain that is independent of nociceptive input, capturing
  endogenous / cognitive contributions to pain.
  [Woo, C.-W. et al. (2017). Quantifying cerebral contributions to pain beyond nociception. Nature Communications, 8, 14211.](https://doi.org/10.1038/ncomms14211)
- **PINES — Picture-Induced Negative Emotion Signature.** Predictive
  pattern for the intensity of negative emotion elicited by aversive
  images.
  [Chang, L.J., Gianaros, P.J., Manuck, S.B., Krishnan, A. & Wager, T.D. (2015). A sensitive and specific neural signature for picture-induced negative affect. PLoS Biology, 13(6), e1002180.](https://doi.org/10.1371/journal.pbio.1002180)
- **VPS — Vicarious Pain Signature.** Pattern predictive of vicarious
  (observed-other) pain, dissociable from first-person NPS.
  [Krishnan, A. et al. (2016). Somatic and vicarious pain are represented by dissociable multivariate brain patterns. eLife, 5, e15166.](https://doi.org/10.7554/eLife.15166)
- **Romantic Rejection signature.** Pattern predictive of social
  rejection-related affect.
  [Woo, C.-W. et al. (2014). Separate neural representations for physical pain and social rejection. Nature Communications, 5, 5380.](https://doi.org/10.1038/ncomms6380)
- **Kragel emotion category signatures.** Seven category-specific
  patterns (amusement, anger, awe, contentment, fear, sadness, sexual
  desire, surprise) trained on naturalistic emotional film/music.
  [Kragel, P.A., Reddan, M.C., LaBar, K.S. & Wager, T.D. (2019). Emotion schemas are embedded in the human visual system. Science Advances, 5(7), eaaw4358.](https://doi.org/10.1126/sciadv.aaw4358)
- **NCS — Neurobiological Craving Signature.** Pattern predictive of
  drug and food craving across modalities.
  [Koban, L., Wager, T.D. & Kober, H. (2023). A neuromarker for drug and food craving distinguishes drug users from non-users. Nature Neuroscience, 26(2), 316-325.](https://doi.org/10.1038/s41593-022-01228-w)
- **MPA2 — Multi-aversive Pattern Analysis.** Patterns generalizing
  across pain, negative affect, and other aversive states.
  [Čeko, M. et al. (2022). Common and stimulus-type-specific brain representations of negative affect. Nature Neuroscience, 25(6), 760-770.](https://doi.org/10.1038/s41593-022-01082-w)
- **Fibromyalgia signatures (FM-pain, FM-multisensory).** Patterns
  diagnostic of fibromyalgia from pain and multisensory stimulation
  fMRI.
  [López-Solà, M. et al. (2017). Towards a neurophysiological signature for fibromyalgia. Pain, 158(1), 34-47.](https://doi.org/10.1097/j.pain.0000000000000707)

Many other signatures (e.g. for autonomic arousal, reward, empathy,
threat, working memory, response conflict) are registered in
`load_image_set`. See [`sample_datasets.md`](sample_datasets.md) for the
exhaustive table.

## References to other sources

- [canlab/Neuroimaging_Pattern_Masks](https://github.com/canlab/Neuroimaging_Pattern_Masks) —
  the GitHub repo distributing the atlas, signature, and meta-analytic
  image files.
- [CANlab Brain Patterns](https://sites.google.com/dartmouth.edu/canlab-brainpatterns/home) —
  curated companion site with descriptions, figures, and usage notes
  for each atlas and signature.
- [canlab.github.io](https://canlab.github.io) — main toolbox docs
  including walkthroughs, tutorials, and per-class object pages.
- [`atlas_methods.md`](atlas_methods.md) — methods for working with an
  `atlas` object.
- [`region_methods.md`](region_methods.md) — methods for working with a
  `region` object.
- [`sample_datasets.md`](sample_datasets.md) — full `load_image_set`
  keyword registry, including signatures and shipped sample datasets.
