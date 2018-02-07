img = which('keuken_2014_enhanced_for_underlay.img');

% sup colliculus
xyz = [5.1767  -32.3069   -4.2592
-5.1767  -32.3069  -4.2592];

myradius = 3.5;

cluster_orthviews();

cl = sphere_roi_tool_2008(img, myradius, xyz, 'useexisting');

r = cluster2region(cl);
obj = region2imagevec(r);

% Mask with brainstem
bstem = fmri_data(which('brainstem_mask_tight_2018.img'));
obj = apply_mask(obj, bstem);
r = region(obj);
orthviews(r)

r(1).shorttitle = 'R_SC';
r(2).shorttitle = 'L_SC';

sc_regions = r;
sc_obj = obj;


%% inf coll

xyz = [4.8186  -37.2304   -9.9917
      -4.8186  -37.2304   -9.9917];

  myradius = 3;

%cluster_orthviews();

cl = sphere_roi_tool_2008(img, myradius, xyz, 'useexisting');

r = cluster2region(cl);
obj = region2imagevec(r);

% Mask with brainstem
bstem = fmri_data(which('brainstem_mask_tight_2018.img'));
obj = apply_mask(obj, bstem);
r = region(obj);
orthviews(r)

r(1).shorttitle = 'R_IC';
r(2).shorttitle = 'L_IC';

ic_regions = r;
ic_obj = obj;


%% Save

cd('/Users/tor/Documents/Code_Repositories/CanlabCore/CanlabCore/canlab_canonical_brains/Thalamus_brainstem_ROIs_surfaces')
savename = 'sup_inf_colliculus_roi_2018_tor';

pag_regions = r;
pag_obj = obj;

save(savename, 'sc_regions', 'sc_obj', 'ic_regions', 'ic_obj');

