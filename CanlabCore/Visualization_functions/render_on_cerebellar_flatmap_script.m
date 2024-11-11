% function render_on_cerebellar_flatmap(target_obj)
% This script contains steps to render an image stored in a nii file, or in
% an fmri_data object, onto a standard cerebellar flat map.
%
% - requires SPM, the SUIT toolbox, and some standard templates on your Matlab path.
% - saves a Gifti .gii image of your map, called <image_name>.cblm_surf.gii
%
% Diedrichsen, J. & Zotow, E. (2015). Surface-based display of volume-averaged cerebellar data. PLoS One, 7, e0133402.
% https://www.diedrichsenlab.org/imaging/suit_flatmap.htm


% g = gifti('FLAT.surf.gii')
% figure; han = patch('Faces',g.faces,'Vertices',g.vertices);
% set(han, 'EdgeColor', 'none');

% 2021 Spisak Placebo meta N = 603 maps
fname = which('full_pain_g_pperm_FWE05.nii.gz');
pain603 = fmri_data(fname);

fname = which('full_pla_g_pperm_tfce_FWE05.nii.gz');
placebo603 = fmri_data(fname);

fname = which('full_pla_rrating_pperm_tfce_FWE05.nii.gz');
placebocorr603 = fmri_data(fname);

gunzip(pain603.fullpath);
pain603fname = strrep(pain603.fullpath, '.gz', '');

gunzip(placebo603.fullpath);
placebo603fname = strrep(placebo603.fullpath, '.gz', '');

gunzip(placebocorr603.fullpath);
placebocorr603fname = strrep(placebocorr603.fullpath, '.gz', '');



targetniifile = '/home/bogdan/Downloads/full_pla_g_unthresh.nii';

% targetniifile = fullfile(pwd, 'nps.nii');

targetniifile = pain603fname;
targetniifile = placebo603fname;
targetniifile = placebocorr603fname;

target_obj = load_image_set('nps');

%%

% Turn target object into a temporary .nii file
% ------------------------------
if size(target_obj.dat, 2) > 1
    error('Target object should contain only a single image')
end

target_obj.fullpath = fullfile(pwd, 'tmp_target_nii.nii');
write(target_obj);

targetniifile = fullfile(pwd, 'tmp_target_nii.nii');
fprintf('created %s\n', targetniifile)

% Check for SUIT toolbox
% ------------------------------
suitname = which('suit_reslice_dartel.m');
if isempty(suitname)
    disp('Can''t find suit_reslice_dartel.m')
    error('Add SUIT folder to path')
end

% Check for SPM toolbox
% ------------------------------
if isempty(which('spm_dartel_norm.m')), error('Start SPM to find spm_dartel_norm.m'), end

% Find SUIT in correct place and add to path
% Note: SUIT must be in spm12/toolbox/suit/flatmap' for suit_map2surf to work
% ------------------------------
spmtoolboxdir = fileparts(fileparts(which('spm_dartel_norm.m')));
suitdir = fullfile(spmtoolboxdir, 'suit');
if ~isdir(suitdir), error('Add suit toolbox to spm12/toolbox/suit'), end
g = genpath(suitdir); addpath(g);

% Check for target .nii file
% ------------------------------
targetniifile = which(targetniifile);       % get full path
if isempty(targetniifile)
    disp('Can''t find target .nii file on path')
    error('Check path and file names.')
end

[dd, ff, ee] = fileparts(targetniifile);

% Set output names
% ------------------------------

wdtargetniifile = fullfile(dd, ['wd' ff ee]);

giftiname = fullfile(dd, [ff '.cblm_surf.gii']);

% Find templates for flat map transformation
% ------------------------------
affinename = which('Affine_MNI152NLin6Asym_T1_1mm_seg1.mat');
if isempty(affinename)
    disp('Can''t find Affine_MNI152NLin6Asym_T1_1mm_seg1.mat')
    error('Add Canlab Neuroimaging_Pattern_Masks github subfolders to path')
end

flowname = which('u_a_MNI152NLin6Asym_T1_1mm_seg1.nii');
if isempty(flowname)
    disp('Can''t find u_a_MNI152NLin6Asym_T1_1mm_seg1.nii')
    error('Add Canlab Neuroimaging_Pattern_Masks github subfolders to path')
end

cerebmaskfile = which('c_MNI152NLin6Asym_T1_1mm_pcereb.nii');
if isempty(cerebmaskfile)
    disp('Can''t find c_MNI152NLin6Asym_T1_1mm_pcereb.nii')
    error('Add Canlab Neuroimaging_Pattern_Masks github subfolders to path')
end

% Some code to create templates above
% ------------------------------
% Normalize template to SUIT cerebellar template - define mapping to SUIT space
% suit_normalize_dartel(struct('subjND', struct('gray', {{'MNI152NLin6Asym_T1_1mm_seg1.nii'}}, 'white', {{'MNI152NLin6Asym_T1_1mm_seg2.nii'}}, 'isolation', {{'c_MNI152NLin6Asym_T1_1mm_pcereb.nii'}})))

% Apply normalization : Reslice the target image, already in MNI space, in SUIT space

% See also:
% apply_spm_warp(mvg_img0, fxd_img0, pre_affine_mat, warp_img, post_affine_mat, out_img, interp)

% Reslice the target image
% ------------------------------

suit_input_struct = struct('subj', struct('affineTr', {{affinename}}, ...
    'flowfield', {{flowname}}, ...
    'resample', {{targetniifile}}, ...
    'mask', {{cerebmaskfile}}));

% creates wdtargetniifile
suit_reslice_dartel(suit_input_struct)
fprintf('created %s\n', wdtargetniifile)

% Map to cerebellar surface
% ------------------------------

C = [];
C.cdata = suit_map2surf(wdtargetniifile, 'space','SUIT', 'stats', @nanmean);

% Save GIFTI gii file, if requested
% ------------------------------
if ~isempty(giftiname)

    save(gifti(C), giftiname)
    fprintf('created %s\n', giftiname)
    
end


% Plot the map
% ------------------------------

create_figure('Cerebellar flatmap')

cmap = colormap_tor([.5 0 1], [1 1 0]);

mymax = prctile(abs(C.cdata), 95);
suit_plotflatmap(C.cdata, 'cmap', cmap, 'cscale',[-mymax mymax]);

% f = gifti(giiname)
% mymax = prctile(abs(f.cdata), 95);
% suit_plotflatmap(f.cdata, 'cmap', cmap, 'cscale',[-mymax mymax]);

% Clean up
% ------------------------------

delete(targetniifile)
delete(wdtargetniifile)



