
cluster_orthviews();

% zoom in to 40 mm

%% Define sphere centers

img = which('keuken_2014_enhanced_for_underlay.img');
% 
% cl = sphere_roi_tool_2008(img, 4, [], 'useexisting');
% 
% xyz = cat(1, cl.mm_center);



%% PAG

xyz = [    0.3892  -29.1969   -6.7374
    0.3946  -30.0161   -7.6839
    0.3955  -31.3315   -9.1658
    0.4000  -32.2671  -10.6220
    0.4000  -33.2795  -11.9429];

cluster_orthviews();

cl = sphere_roi_tool_2008(img, 4, xyz, 'useexisting');

r = cluster2region(cl);
obj = region2imagevec(r);

% Mask with brainstem
bstem = fmri_data(which('brainstem_mask_tight_2018.img'));

obj = apply_mask(obj, bstem);
r = region(obj);

orthviews(r)
%% Save

cd('/Users/tor/Documents/Code_Repositories/CanlabCore/CanlabCore/canlab_canonical_brains/Thalamus_brainstem_ROIs_surfaces')
savename = 'pag_roi_2018_tor';

r(1).shorttitle = 'PAG';

pag_regions = r;
pag_obj = obj;

save(savename, 'pag_regions', 'pag_obj');
