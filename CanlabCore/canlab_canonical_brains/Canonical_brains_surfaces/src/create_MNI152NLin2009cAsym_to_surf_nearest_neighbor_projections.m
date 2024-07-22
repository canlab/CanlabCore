addpath(genpath('/home/bogdan/.matlab/canlab/CanlabCore')) % for barplot_columns
addpath('/home/bogdan/.matlab/spm/spm12'); % for giftis

addpath(genpath('/home/bogdan/.matlab/canlab/Neuroimaging_Pattern_Masks')) % for barplot_columns
close all; clear all;

% we need to load up some data for testing these projections. the
% destrieux, dkt or desikan killiany atlases are all good choices because
% they're designed to label distinct anatomical features (gyri and sulci),
% and should should be easily evaluated in terms of accuracy on the
% cortical surface projections we obtain.
atlas_obj = load_atlas('destrieux_fsl6');

%% format surface meshes for canlab
% we got the surface meshes from either the stnadard freesurfer's subject 
% folder (fsaverage meshes), running recon-all on a template (MNI surfaces)
% or from HCP. Surfaces were converted from freesurfer format to gifti
% files using freesurfer's mris_convert where needed.

% template white
create_canlab_surf_file(fullfile('templates','MNI152NLin2009cAsym_white_lh.surf.gii'), 'MNI152NLin2009cAsym_white_lh.mat');
create_canlab_surf_file(fullfile('templates','MNI152NLin2009cAsym_white_rh.surf.gii'), 'MNI152NLin2009cAsym_white_rh.mat');

% template sphere
create_canlab_surf_file(fullfile('templates','MNI152NLin2009cAsym_sphere.reg_lh.surf.gii'), 'MNI152NLin2009cAsym_sphere.reg_lh.mat');
create_canlab_surf_file(fullfile('templates','MNI152NLin2009cAsym_sphere.reg_rh.surf.gii'), 'MNI152NLin2009cAsym_sphere.reg_rh.mat');

% HCP ref sphere (rotated relative to fsaverage)
% the sphere is aligned to fsaverage space
% source: https://github.com/rmldj/hcp-utils
create_canlab_surf_file(fullfile('/home/bogdan/Downloads/hcp_utils/hcp_utils/data/','S1200.L.sphere.32k_fs_LR.surf.gii'), 'S1200.L.sphere.32k_fs_LR.mat');
create_canlab_surf_file(fullfile('/home/bogdan/Downloads/hcp_utils/hcp_utils/data/','S1200.R.sphere.32k_fs_LR.surf.gii'), 'S1200.R.sphere.32k_fs_LR.mat');

% HCP ref sphere in fsaverage space (not rotated)
% the sphere is aligned to fsaverage space and not rotated
% source:  https://github.com/Washington-University/HCPpipelines/tree/master/global/templates/standard_mesh_atlases
create_canlab_surf_file(fullfile('templates','fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii'), 'fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.mat');
create_canlab_surf_file(fullfile('templates','fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.surf.gii'), 'fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.mat');

% freesurfer ref sphere
% the sphere is aligned to fsaverage space
create_canlab_surf_file(fullfile('templates','fsavg_sphere_lh.surf.gii'), 'fsavg_sphere_lh.mat');
create_canlab_surf_file(fullfile('templates','fsavg_sphere_rh.surf.gii'), 'fsavg_sphere_rh.mat');


%% interpolate volume to white surface
%src_img = mpa2.get_wh_image(3);
src_img = atlas_obj;
tic;
surf_white_lh = load(which('MNI152NLin2009cAsym_white_lh.mat'),'faces','vertices');
[~,~,mesh_struct] = reconstruct_image(src_img);
c_white_lh = interp3(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, mesh_struct.voldata, ...
            surf_white_lh.vertices(:,1), surf_white_lh.vertices(:,2), surf_white_lh.vertices(:,3),'nearest');
surf_white_rh = load(which('MNI152NLin2009cAsym_white_rh.mat'));
[~,~,mesh_struct] = reconstruct_image(src_img);
c_white_rh = interp3(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, mesh_struct.voldata, ...
            surf_white_rh.vertices(:,1), surf_white_rh.vertices(:,2), surf_white_rh.vertices(:,3),'nearest');
toc
%% freesurfer resampling
tic
sphere_reg_lh = load('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/MNI152NLin2009cAsym_sphere.reg_lh.mat');
sphere_reg_lh.vertices = double(sphere_reg_lh.vertices);
freesurfer_sphere_lh = load('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/fsavg_sphere_lh.mat');
freesurfer_sphere_lh.vertices = double(freesurfer_sphere_lh.vertices);
[c_fs_sphere_lh, BcC_MNI6_to_fsavg_lh, BcV_MNI6_to_fsavg_lh] = spherical_icosahedral_interpolation(sphere_reg_lh.vertices - mean(sphere_reg_lh.vertices), sphere_reg_lh.faces, c_white_lh, freesurfer_sphere_lh.vertices,'nearest');

sphere_reg_rh = load('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/MNI152NLin2009cAsym_sphere.reg_rh.mat');
sphere_reg_rh.vertices = double(sphere_reg_rh.vertices);
freesurfer_sphere_rh = load('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/fsavg_sphere_rh.mat');
freesurfer_sphere_rh.vertices = double(freesurfer_sphere_rh.vertices);
[c_fs_sphere_rh, BcC_MNI6_to_fsavg_rh, BcV_MNI6_to_fsavg_rh] = spherical_icosahedral_interpolation(sphere_reg_rh.vertices, sphere_reg_rh.faces, c_white_rh, freesurfer_sphere_rh.vertices,'nearest');

%% HCP resampling
% note that this is not using the standard HCP sphere. For whatever reason
% the standard HCP sphere is rotated relative to fsaverage. This sphere
% here was obtained from
% https://github.com/Washington-University/HCPpipelines/tree/master/global/templates/standard_mesh_atlases
% and is rotated relative to the standard HCP sphere for alignment with
% fsaverage. Think of the standard HCP sphere as analagous to a subject
% specific *.sphere.* file from freesurfer recon-all outputs and the sphere
% below as equivalent to *.sphere.reg.* outputs from recon-all.

fsLR_sphere_lh = load('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/S1200.L.sphere.32k_fs_LR.mat');
fsLR_sphere_lh.vertices = double(fsLR_sphere_lh.vertices);
fsLR_to_fs_sphere_lh = load('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.mat');
fsLR_to_fs_sphere_lh.vertices = double(fsLR_to_fs_sphere_lh.vertices);
[c_hcp_sphere_lh, BcC_fsavg_to_fslr32_lh, BcV_fsavg_to_fslr32_lh]  = spherical_icosahedral_interpolation(fsLR_to_fs_sphere_lh.vertices, fsLR_to_fs_sphere_lh.faces, c_fs_sphere_lh, fsLR_sphere_lh.vertices,'nearest');

fsLR_sphere_rh = load('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/S1200.R.sphere.32k_fs_LR.mat');
fsLR_sphere_rh.vertices = double(fsLR_sphere_rh.vertices);
fsLR_to_fs_sphere_rh = load('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.mat');
fsLR_to_fs_sphere_rh.vertices = double(fsLR_to_fs_sphere_rh.vertices);
[c_hcp_sphere_rh, BcC_fsavg_to_fslr32_rh, BcV_fsavg_to_fslr32_rh] = spherical_icosahedral_interpolation(fsLR_to_fs_sphere_rh.vertices, fsLR_to_fs_sphere_rh.faces, c_fs_sphere_rh, fsLR_sphere_rh.vertices,'nearest');
toc

% evaluate projection

o2 = src_img.montage('freesurfer inflated') % thermal
f1 = gcf;
f2 = figure;
copyobj(f1.Children,f2)

[c_lh_colored, c_rh_colored] = map_to_colors(c_fs_sphere_lh, c_fs_sphere_rh, o2, 'doindexmap');

c_orig = o2.surface{1}.object_handle.FaceVertexCData;
o2.surface{1}.object_handle.FaceVertexCData = c_rh_colored;
o2.surface{2}.object_handle.FaceVertexCData = c_lh_colored;
o2.surface{3}.object_handle.FaceVertexCData = c_rh_colored;
o2.surface{4}.object_handle.FaceVertexCData = c_lh_colored;

% evaluate projection

o2 = src_img.montage('hcp inflated') % thermal
f1 = gcf;
f2 = figure;
copyobj(f1.Children,f2)

[c_lh_colored, c_rh_colored] = map_to_colors(c_hcp_sphere_lh, c_hcp_sphere_rh, o2, 'doindexmap');

c_orig = o2.surface{1}.object_handle.FaceVertexCData;
o2.surface{1}.object_handle.FaceVertexCData = c_rh_colored;
o2.surface{2}.object_handle.FaceVertexCData = c_lh_colored;
o2.surface{3}.object_handle.FaceVertexCData = c_rh_colored;
o2.surface{4}.object_handle.FaceVertexCData = c_lh_colored;

%% save projection

weights_lh = zeros(size(BcC_fsavg_to_fslr32_lh,1),1);
vertices_lh = zeros(size(BcC_fsavg_to_fslr32_lh,1),1);
val_lh = zeros(size(BcC_fsavg_to_fslr32_lh,1),1);
for target_vx = 1:size(BcC_fsavg_to_fslr32_lh,1)
    MNI_vx = BcV_MNI6_to_fsavg_lh(BcV_fsavg_to_fslr32_lh(target_vx,:),:);
    MNI_weights = BcC_MNI6_to_fsavg_lh(BcV_fsavg_to_fslr32_lh(target_vx,:),:).*BcC_fsavg_to_fslr32_lh(target_vx,:)';
    
    weights_lh(target_vx,:) = MNI_weights(:);
    vertices_lh(target_vx,:) = MNI_vx(:);

    val_lh(target_vx) = MNI_weights(:)'*c_white_lh(MNI_vx(:));
end

weights_rh = zeros(size(BcC_fsavg_to_fslr32_rh,1),1);
vertices_rh = zeros(size(BcC_fsavg_to_fslr32_rh,1),1);
val_rh = zeros(size(BcC_fsavg_to_fslr32_rh,1),1);
for target_vx = 1:size(BcC_fsavg_to_fslr32_rh,1)
    MNI_vx = BcV_MNI6_to_fsavg_rh(BcV_fsavg_to_fslr32_rh(target_vx,:),:);
    MNI_weights = BcC_MNI6_to_fsavg_rh(BcV_fsavg_to_fslr32_rh(target_vx,:),:).*BcC_fsavg_to_fslr32_rh(target_vx,:)';
    
    weights_rh(target_vx,:) = MNI_weights(:);
    vertices_rh(target_vx,:) = MNI_vx(:);

    val_rh(target_vx) = MNI_weights(:)'*c_white_rh(MNI_vx(:));
end
    

% QC, should be the same as the last plot
o2 = src_img.montage('hcp inflated'); % thermal
f1 = gcf;
f2 = figure;
copyobj(f1.Children,f2)

val_lh = diag(weights_lh*c_white_lh(vertices_lh)');
val_rh = diag(weights_rh*c_white_rh(vertices_rh)');
[c_lh_colored, c_rh_colored] = map_to_colors(val_lh, val_rh, o2,'doindexmap');

c_orig = o2.surface{1}.object_handle.FaceVertexCData;
o2.surface{1}.object_handle.FaceVertexCData = c_rh_colored;
o2.surface{2}.object_handle.FaceVertexCData = c_lh_colored;
o2.surface{3}.object_handle.FaceVertexCData = c_rh_colored;
o2.surface{4}.object_handle.FaceVertexCData = c_lh_colored;

descrip = 'this struct provides nearest neighbor vertices indexed from CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/MNI152NLin2009cAsym_sphere.reg_lh.mat (a freesurfer generated surface reconstruction of MNI152NLin2009cAsym) to remap vertices into HCP''s fs_LR 32k template space.';
resample_from_MNI152NLin2009cAsym_to_fsavg_164k = struct('weights_lh', BcC_MNI6_to_fsavg_lh, 'vertices_lh', BcV_MNI6_to_fsavg_lh, 'weights_rh', BcC_MNI6_to_fsavg_rh, 'vertices_rh', BcV_MNI6_to_fsavg_rh, 'description', descrip);

descrip = 'this struct provides nearest neighbor vertices indexed into CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/MNI152NLin2009cAsym_sphere.reg_lh.mat (a freesurfer generated surface reconstruction of MNI152NLin2009cAsym) to remap vertices into HCP''s fs_LR 32k template space.';
resample_from_MNI152NLin2009cAsym_to_fsLR_32k = struct('weights_lh', weights_lh, 'vertices_lh', vertices_lh, 'weights_rh', weights_rh, 'vertices_rh', vertices_rh, 'description', descrip);

save('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/resample_from_MNI152NLin2009cAsym_to_fsLR_32k_nearestneighbor.mat', 'resample_from_MNI152NLin2009cAsym_to_fsLR_32k');
save('/home/bogdan/.matlab/canlab/CanlabCore/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/resample_from_MNI152NLin2009cAsym_to_fsavg_164k_nearestneighbor.mat', 'resample_from_MNI152NLin2009cAsym_to_fsavg_164k');