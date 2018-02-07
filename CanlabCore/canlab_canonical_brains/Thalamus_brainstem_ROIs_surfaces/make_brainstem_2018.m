mask = fmri_data(which('canlab_atlas_brainstem.nii'));
img = fmri_data(which('CIT168toMNI152_T1w_head_1mm.nii'));
mask = resample_space(mask, img);
img = apply_mask(img, mask);
img = threshold(img, [1.7 Inf], 'raw-between');
orthviews(img)

%% remove small clusters

r = region(img);
r = r(1);

img2 = region2fmri_data(r, img);

% omit above
%     5.1885
%   -22.5394
%    -2.4894
   
%% 

figure; 
isosurface(img)

%% Remove a bunch of edge stuff

imgname = which('CIT168toMNI152_T1w_head_1mm.nii');

%[cl, all_data] = sphere_roi_tool_2008(imgname, 4, [], 'useexisting');

xyz = [11.6562  -37.0938  -21.0312
   11.6216  -38.6216  -28.3243
   19.6216  -25.6757  -18.6216
   19.6389  -24.2778  -15.3611
   0.8582  -42.6978  -21.8060
   14.9179  -34.0224  -21.1716
   14.8951  -35.6517  -23.2060
   18.5000  -35.6370  -26.7407
  -13.6493  -37.5560  -26.7239
   13.9251  -30.7453   -3.7228
   11.0000  -30.8052   -0.1461
   20.0830  -30.7585  -24.4868
   23.3820  -30.7453  -29.3820
   12.9361  -42.7481  -29.3571
   12.9363  -39.5243  -24.4906
   10.6741  -41.1852  -24.5000
    9.3941  -41.2082  -19.6059
   12.6231  -37.8806  -19.6119
   12.6222  -35.6370   -7.8963
   12.5404  -37.5404  -12.4877
    6.0947   -4.1098   -5.9697
   11.2873   -6.1119   -5.9701
   16.5346   -9.0154   -5.9615
   21.3698  -12.2830   -5.9472
   24.6679  -17.4403   -5.9366
   18.0712  -40.2547  -28.7491
   23.7361  -40.2528  -33.9033
   24.5132  -36.2415  -33.9170
  -15.7500  -39.4552  -29.8619
   18.0749  -44.2884  -45.1948
   20.5000  -37.8134  -45.1866
  -17.3407  -45.1556  -43.5778
   12.4772   -2.4526  -14.5228
   12.5000   -4.0000   -9.7328
   -5.2943   -4.0000  -10.5321
    1.1536   -4.0000  -11.2846
    5.2424   -4.0000  -10.5568
    8.4091   -4.0000  -11.3182
   26.9331  -15.2528   -8.1338
   19.6900   -9.6384   -8.1439
   23.6815  -15.3000  -10.5296
   23.6891  -15.3109   -4.9101
  -14.9660  -32.1358   -4.9321
   18.0418  -29.8745   -4.9582
   26.1567  -22.5560   -9.7313
   25.3630  -22.5516   -6.5231
    3.6007  -22.5746   -0.1343
    5.1992  -28.1992   -0.1090
    4.4050   -9.5950   -2.4875
    4.4067  -13.6343    0.7201
    4.3993  -18.5092    0.7289
   23.7068  -23.3271   -8.9323];

cl = sphere_roi_tool_2008(imgname, 4, xyz, 'useexisting');
   
r_remove = cluster2region(cl);
img_remove = region2fmri_data(r_remove, img2);
img2.dat(img_remove.dat > 0) = 0;

orthviews(region(img2), 'color', {[0 0 1]})

%% Smooth, threshold and save

img3 = preprocess(img2, 'smooth', 3);

img3 = threshold(img3, [.2 Inf], 'raw-between');

myname = fullfile('/Users/tor/Documents/Code_Repositories/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces', 'canlab_brainstem.img');
img3.fullpath = myname;
% write(img3)

%%

figure; 
p = isosurface(img3, 'sd', 2, 'thresh', .1);
set(p, 'Facecolor', [0 0 1]);
addbrain('hires left')

drawnow

%% Make tight mask for masking things
% Mask with brainstem
% probabilistic mask is still a bit off/asymmetrical.
% use Keuken data to make it match.
bstem = atlas(which('canlab_brainstem.img'));
bstemimg = fmri_data(which('keuken_2014_enhanced_for_underlay.img'));
bstemimg = apply_mask(bstemimg, bstem); % values from Keuken, existing mask
bstemimg = threshold(bstemimg, [120 Inf], 'raw-between');

bstemimg.fullpath = fullfile(pwd, 'brainstem_mask_tight_2018.img');
%write(bstemimg);

