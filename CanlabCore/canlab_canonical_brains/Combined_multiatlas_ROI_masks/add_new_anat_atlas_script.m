name = '/Users/tor/Downloads/Cerebellum-MNIsegment-MRICroN/Cerebellum-MNIsegment.nii'
labelimg = which('anat_lbpa_thal.img'); %% 'lpba40.spm5.avg152T1.label.nii');
[volInfo, labels] = iimg_read_img(labelimg, 2);
figure; plot(labels)
newlabels = scn_map_image(name, labelimg);
newlabels = round(newlabels(:));

%% show us what we're replacing

disp('Replacing some or all of the voxels originally labeled with these numbers:')
u = unique(labels(newlabels > 0))
spm_image('init', labelimg)

clear cl
for i = 1:length(u)
    
    replaced_regions = zeros(size(labels));
replaced_regions(labels == u(i)) = u(i);

cl{i} = iimg_indx2clusters(replaced_regions, volInfo);

cluster_orthviews(cl{i}, {rand(1, 3)}, 'add');

end

clnew = iimg_indx2clusters(newlabels, volInfo);

cluster_orthviews(clnew, {[1 0 0]}, 'add', 'solid');

%% do the replacement
% add an offset in integer numbers so that we don't have mutliple regions
% with the same labels

offset_num = max(max(u) + 1, 200);  

labels(newlabels > 0) = offset_num + newlabels(newlabels > 0);

%% write output

[dd, ff, ee] = fileparts(labelimg);

outname = [dd filesep 'atlas_labels_combined.img']

 iimg_reconstruct_vols(labels, volInfo, 'outname', outname);
 
%% check it

names = char(which('T1.nii'), labelimg, name, outname);
spm_check_registration(names)

%% custom for CEREBELLAR diedrichsen atlas

wh_cblm = labels == 181; sum(wh_cblm)  
% these are in the original cerebellum in LBPA40, but outside the Diedrichesen 2009 atlas in SPM5 segment space

cl = iimg_indx2clusters(wh_cblm, volInfo);
% white matter or outside new cerebellar atlas

