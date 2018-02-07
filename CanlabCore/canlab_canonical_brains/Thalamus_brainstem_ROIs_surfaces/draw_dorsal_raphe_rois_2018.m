
% From Beliveau 2015 Neuroimage
% Beliveau, Vincent, Claus Svarer, Vibe G. Frokjaer, Gitte M. Knudsen, Douglas N. Greve, and Patrick M. Fisher. 2015. ?Functional Connectivity of the Dorsal and Median Raphe Nuclei at Rest.? NeuroImage 116 (August). Elsevier:187?95.
%
% dorsal raphe (DR) and median raphe (MR)
% Histological studies performed by Baker et al. (1991a, 1991b, 1990)
% have provided in-depth knowledge of the morphology and location of
% the DR and the MR in the ex vivo human brain. However, to perform
% seed-based FC, accurate in vivo segmentation of DR and MR are needed
% (Kalbitzer and Svarer, 2009). This presents a challenge (Kranz and
% Hahn, 2012), as the raphe nuclei are composed of sparse neurons
% surrounded by white matter and they have no well-defined boundaries
% visible in MRI (Baker et al., 1996, 1991a, 1990).
% We have adopted a method similar to Schain et al. (2013) in which
% liberal search volumes were defined on the structural MRI and then
% refined using the PET image. The DR lies on the midline of the brainstem
% and extends from the oculomotor nucleus to the middle of the pons
% (Baker et al., 1990). It can be subdivided at the level of the isthmus
% into two groups, a midbrain (B7) group and a pontine (B6) group
% (Dahlström and Fuxe, 1964) which meet near the inferior opening of
% the cerebral aqueduct (CA). The B7 group is adjacent to the CA. The B6
% group is only about 0.5 mm in radius, well below current scanner resolution
% for fMRI. For this reason, we focused on the B7 group as the seed
% region for our analysis. 

% The average centroid of the DR was (0, ? 31 ? 9) in MNI305 space; 
%that of the MR was (0, ? 35, ? 21). Although the volume of the seeds 
% defined on the PET images was constant, the reintroduction of GD slightly 
% affected the final volume of the seeds from subject to subject; the volume 
% (mean ± std) was 118 ± 11 mm3 for DR and 65 ± 8 mm3 for MR. For the same reason, 
% the number of functional voxels was 18 ± 3 for DR and 11 ± 2 for MR.

myradius = 2;

cluster_orthviews();
xyz = [0 -31 -9];
cl = sphere_roi_tool_2008(img, myradius, xyz, 'useexisting');

r = cluster2region(cl);
obj = region2imagevec(r);

% Mask with brainstem
bstem = fmri_data(which('brainstem_mask_tight_2018.img'));
obj = apply_mask(obj, bstem);
r = region(obj);
orthviews(r)

r(1).shorttitle = 'Dorsal_raphe_DR';

drn_regions = r;
drn_obj = obj;

%%
xyz = [0 -35 -21];
cl = sphere_roi_tool_2008(img, myradius, xyz, 'useexisting');

r = cluster2region(cl);
obj = region2imagevec(r);

% Mask with brainstem
%bstem = fmri_data(which('brainstem_mask_tight_2018.img'));
obj = apply_mask(obj, bstem);
r = region(obj);
orthviews(r)

r(1).shorttitle = 'Median_raphe_MR';

mrn_regions = r;
mrn_obj = obj;

%% Save

cd('/Users/tor/Documents/Code_Repositories/CanlabCore/CanlabCore/canlab_canonical_brains/Thalamus_brainstem_ROIs_surfaces')
savename = 'dorsal_median_raphe_roi_2018_tor';

pag_regions = r;
pag_obj = obj;

save(savename, 'drn_regions', 'drn_obj', 'mrn_regions', 'mrn_obj');
