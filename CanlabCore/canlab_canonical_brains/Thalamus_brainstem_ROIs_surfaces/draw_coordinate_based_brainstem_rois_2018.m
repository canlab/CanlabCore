
% See bottom for notes!!!


% Nash, Paul G., Vaughan G. Macefield, Iven J. Klineberg, Greg M. Murray, and Luke A. Henderson. 2009. ?Differential Activation of the Human Trigeminal Nuclear Complex by Noxious and Non-Noxious Orofacial Stimulation.? Human Brain Mapping 30 (11):3772?82.

medullary_raphe = [0 -36 -50];

spinal_trigeminal = [4 -38 -52; -4 -38 -52];

% Beissner, Florian, Ralf Deichmann, and Simon Baudrexel. 2011. ?fMRI of the Brainstem Using Dual-Echo EPI.? NeuroImage 55 (4). Elsevier:1593?99.
% can't use, not MNI
% trigeminal_motor = [6 4 -6; -6 4 -6];  
% facial_nuc = [7 3 -1; -7 -3 -1];

% Sclocco, Roberta, Florian Beissner, Gaelle Desbordes, Jonathan R. Polimeni, Lawrence L. Wald, Norman W. Kettner, Jieun Kim, et al. 2016. ?Neuroimaging Brainstem Circuitry Supporting Cardiovagal Response to Pain: A Combined Heart Rate Variability/ultrahigh-Field (7 T) Functional Magnetic Resonance Imaging Study.? Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences 374 (2067). rsta.royalsocietypublishing.org. https://doi.org/10.1098/rsta.2015.0189.

nuc_ambiguus = [2 -40 -62; 2 -40 -62];
dmnx_nts = [5 -42 -61; -5 -42 -61];

% Bär, Karl-Jürgen, Feliberto de la Cruz, Andy Schumann, Stefanie Koehler, Heinrich Sauer, Hugo Critchley, and Gerd Wagner. 2016. ?Functional Connectivity and Network Analysis of Midbrain and Brainstem Nuclei.? NeuroImage 134 (July):53?63.

ncs_B6_B8 = [0, -32, -24];
nrp_B5 = [0 -34 -30];
nrm = [0 -34 -40];

% Fairhurst, Merle, Katja Wiech, Paul Dunckley, and Irene Tracey. 2007. ?Anticipatory Brainstem Activity Predicts Neural Processing of Pain in Humans.? Pain 128 (1-2):101?10.
% averaged across two coords reported
pbn = [7 -36 -24; -7 -36 -24];

% Zambreanu, L., R. G. Wise, J. C. W. Brooks, G. D. Iannetti, and I. Tracey. 2005. ?A Role for the Brainstem in Central Sensitisation in Humans. Evidence from Functional Magnetic Resonance Imaging.? Pain 114 (3):397?407.
% nucleus cuneiformis (NCF)
ncf = [10 -28 -18; -10 -28 -18];

xyz_mm = {pbn ncf medullary_raphe spinal_trigeminal nuc_ambiguus dmnx_nts ncs_B6_B8 nrp_B5 nrm};
names = {'pbn' 'ncf' 'medullary_raphe' 'spinal_trigeminal' 'nuc_ambiguus' 'dmnx_nts' 'ncs_B6_B8' 'nrp_B5' 'nrm'};

%% Prep

bstem = fmri_data(which('brainstem_mask_tight_2018.img'));

clear *regions
clear *obj
%% Make regions

myradius = 2.5;

cluster_orthviews();

for i = 1:length(names)
    
cl = sphere_roi_tool_2008(img, myradius, xyz_mm{i}, 'useexisting');

r = cluster2region(cl);
obj = region2imagevec(r);

% Mask with brainstem
obj = apply_mask(obj, bstem);
r = region(obj);
orthviews(r)

if length(r) == 2
    r(1).shorttitle = ['R_' names{i}];
    r(2).shorttitle = ['L_' names{i}];
else
    for j = 1:length(r)  
    r(j).shorttitle = names{i};
    end
end

eval([names{i} '_regions = r;'])
eval([names{i} '_obj = obj;'])

end

%% Save

cd('/Users/tor/Documents/Code_Repositories/CanlabCore/CanlabCore/canlab_canonical_brains/Thalamus_brainstem_ROIs_surfaces')
savename = 'coordinate_brainstem_rois_2018_tor';

save(savename, '*_regions', '*_obj');
%%

% Son, Y. D., Z. H. Cho, E. J. Choi, J. H. Kim, H. K. Kim, S. Y. Lee, J. G. Chi, C. W. Park,
% J. H. Kim, and Y. B. Kim. 2014. 'Individually differentiated serotonergic raphe
% nuclei measured with brain PET/MR imaging', Radiology, 272: 541-8.
% Son, Y. D., Z. H. Cho, H. K. Kim, E. J. Choi, S. Y. Lee, J. G. Chi, C. W. Park, and Y. B.
% Kim. 2012. 'Glucose metabolism of the midline nuclei raphe in the brainstem
% observed by PET-MRI fusion imaging', NeuroImage, 59: 1094-7

% Bär, Karl-Jürgen, Feliberto de la Cruz, Andy Schumann, Stefanie Koehler, Heinrich Sauer, Hugo Critchley, and Gerd Wagner. 2016. ?Functional Connectivity and Network Analysis of Midbrain and Brainstem Nuclei.? NeuroImage 134 (July):53?63.
% alternate:
% nucleus raphes dorsalis (DRN, B7, MNI-coordinates, x = 2, y = -26, z = -18), 
% the nucleus centralis superior (B6 + B8,MNI-coordinates, x = 0, y = ?32, z = ?24) 
% and the nucleus raphes pontis (B5, MNI-coordinates, x = 0, y = ?34, z = ?30). 
% The nucleus raphes magnus [would overlap with RVM] was drawn as a box (B3, 6 Å~ 6 Å~ 12 mm centered at MNI-coordinates, x = 0, y = ?34, z = ?40

% Sclocco, 2017 review:
% Sclocco, Roberta, Florian Beissner, Marta Bianciardi, Jonathan R. Polimeni, and Vitaly Napadow. 2017. ?Challenges and Opportunities for Brainstem Neuroimaging with Ultrahigh Field MRI.? NeuroImage, February. Elsevier. https://doi.org/10.1016/j.neuroimage.2017.02.052.

% Current probabilistic in vivo brain atlases (Keuken et al. 2014; Destrieux et al. 2010;
% Desikan et al. 2006; Tzourio-Mazoyer et al. 2002) released with common neuroimaging
% software (e.g FSL, FreeSurfer, SPM) include several cortical and subcortical regions,
% but not most brainstem nuclei nor their subdivisions. Atlases including a limited number
% of brainstem nuclei, such as the substantia nigra (SN), red nucleus (RN), and
% subthalamic nucleus (STh) (Keuken et al. 2014; Chowdhury et al. 2013; Kwon et al.
% 2012; Menke et al. 2010; Mori et al. 2009) do exist, though their subdivisions are not
% reported. Mori et al. developed and publicly released an ex vivo diffusion tensor imaging
% (DTI) atlas of brainstem white matter tracts (Mori et al. 2008; Mori et al. 2009) and in
% recent work, this group also showed the feasibility of ex vivo DTI imaging of several
% nuclei important for motor and cranial nerve functional systems in a single postmortem
% brainstem specimen (Aggarwal et al. 2013).

%  Son and colleagues (Son et al. 2012) successfully exploited
% fluorodeoxyglucose (FDG) PET imaging to identify four distinct clusters with high
% metabolic activity located in the brainstem midline region, consistent with putative raphe
% nuclei. UHF-MRI images were then fused with the metabolic data to produce a single
% image with major brainstem anatomical landmarks, which allowed for identification of
% clusters putatively consistent with the Ncl. raphe dorsalis, pontis, magnus, and
% reticularis centralis superior. Importantly, a follow-up study from the same group
% reported significant correlation between the standard uptake of FDG and
% nondisplaceable binding potential of DASP, a radioligand designed to study the
% serotonin transporter (Son et al. 2014). Both measures identified a group of clusters
% consistent with the previous results, but also including the smaller medullary raphe
% obscurus and pallidus.

% Interestingly, a
% recent study exploited the presence of neuromelanin, a pigment produced in
% noradrenergic neurons that exhibits ferrous properties (Enochs et al. 1997), in order to
% visualize the LC using a T1-weighted Turbo Spin Echo sequence (Sasaki et al. 2006).
% This approach was subsequently extended to a larger cohort (N = 44), in order to obtain
% a probabilistic map of LC in standard MNI space (Keren et al. 2009). Using this regionof-interest
% to guide localization, a recent 7T fMRI study found increased temporal
% correlation between LC fMRI signal and the high-frequency component of heart rate 
% variability, an index of parasympathetic activity, during sustained evoked-pain stimulus
% (Sclocco et al. 2016)

%  neurons containing acetylcholine (ACh) abound
% particularly in the central portions of the reticular core. In the midbrain, a clinically
% important cholinergic center is the pedunculopontine nucleus (PPN), corresponding to
% the CH5 group, as defined by Mesulam and co-authors (Mesulam et al. 1984). The PPN
% is part of the mesencephalic locomotor region, which also includes the cuneiform and
% subcuneiform nuclei. It mainly gives rise to ascending projections, particularly to the
% pars compacta of the substantia nigra and the subthalamic nucleus. stereotactic atlases usually identify only the rostral component of PPN (Zrinzo
% et al. 2008). However, new studies suggest the possibility of an increased effectiveness
% of DBS in the caudal portion (Thevathasan et al. 2012).  the PPN has been shown to be involved in motor planning and activating
% during imagined gait, unlike the subthalamic nucleus (Lau et al. 2015).

% A recent 3T fMRI study using restingstate
% fMRI to investigate functional connectivity in PAG revealed connections between
% the vlPAG and brain regions associated with descending pain modulation (anterior
% cingulate cortex, upper pons/medulla), whereas lPAG and dlPAG were connected with
% brain regions implicated in executive functions (prefrontal cortex, striatum,
% hippocampus) (Coulombe et al. 2016).
%  Faull et al. explored the role of
% PAG in respiratory control using UHF fMRI with 1 mm isotropic voxel size (Faull et al.
% 2015). Their results show deactivation in lPAG and dmPAG columns during short
% (around 10 s) breath hold blocks, suggesting the involvement of these PAG subdivisions
% for conscious respiratory control.  lPAG involvement in resistive respiratory loading, whereas vlPAG was
% associated to anticipation of breathlessness (Faull et al. 2016). 
    


