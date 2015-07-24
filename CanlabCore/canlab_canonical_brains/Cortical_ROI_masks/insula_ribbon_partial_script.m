% procedure: use draw_anatomical_ROI_2008 to draw saggital slices
% then mask with cortical ribbon from segmentation
% then clean up using draw...

cl = mask2clusters('insula_ribbon_R.img');
cl = mask2clusters('SPM8_insula_ribbon_L.img');
mask_union([],'SPM8_insula_ribbon_LR.img','SPM8_insula_ribbon_L.img','insula_ribbon_R.img');
dat = fmri_data('SPM8_insula_ribbon_LR.img'); orthviews(dat)

dat = apply_mask(dat, fmri_data('SPM8_colin27T1_cortical_ribbon.img'));
orthviews(dat)

dat.fullpath = fullfile(pwd, 'insula_ribbon1.img');
write(dat);

%% clean up

draw_anatomical_roi_2008('load', 'insula_ribbon1.img')
draw_anatomical_roi_2008('moveslice')
draw_anatomical_roi_2008('write')
SPM8_insula_ribbon_LR.img

% smooth and threshold
dat = fmri_data('SPM8_insula_ribbon_LR.img')
dat = preprocess(dat, 'smooth', 1);
dat = threshold(dat, [.02 Inf], 'raw-between');
orthviews(dat)


%%
% viz
cl = mask2clusters('SPM8_insula_ribbon_LR.img');

cluster_orthviews(cl, {[1 .7 0]});
cluster_orthviews(cl, {[1 .7 0]}, 'add');

spm_orthviews('xhairs', 'off');
delete(findobj(gcf, 'Type', 'text'));


%% 
dat = fmri_data('SPM8_insula_ribbon_LR.img')
orthviews(dat)
surface(dat)