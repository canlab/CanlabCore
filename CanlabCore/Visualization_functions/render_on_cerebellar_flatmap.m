
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

targetniifile = fullfile(pwd, 'nps.nii');

targetniifile = pain603fname;
targetniifile = placebo603fname;
targetniifile = placebocorr603fname;

create_figure('Cerebellar flatmap')

%%

targetniifile = which(targetniifile);       % get full path
if isempty(targetniifile)
    disp('Can''t find target .nii file on path')
    error('Check path and file names.')
end

[dd, ff, ee] = fileparts(targetniifile);

wdtargetniifile = fullfile(dd, ['wd' ff ee]);

giiname = fullfile(dd, [ff '.func.gii']);

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

if isempty(which('spm_dartel_norm.m')), error('Start SPM to find spm_dartel_norm.m'), end

% Normalize template to SUIT cerebellar template - define mapping to SUIT space
% suit_normalize_dartel(struct('subjND', struct('gray', {{'MNI152NLin6Asym_T1_1mm_seg1.nii'}}, 'white', {{'MNI152NLin6Asym_T1_1mm_seg2.nii'}}, 'isolation', {{'c_MNI152NLin6Asym_T1_1mm_pcereb.nii'}})))

% Apply normalization : Reslice the target image, already in MNI space, in SUIT space

% See also:
% apply_spm_warp(mvg_img0, fxd_img0, pre_affine_mat, warp_img, post_affine_mat, out_img, interp)

suit_input_struct = struct('subj', struct('affineTr', {{affinename}}, ...
    'flowfield', {{flowname}}, ...
    'resample', {{targetniifile}}, ...
    'mask', {{cerebmaskfile}}));


suit_reslice_dartel(suit_input_struct)

C = [];
C.cdata = suit_map2surf(wdtargetniifile, 'space','SUIT', 'stats', @nanmean);

save(gifti(C), giiname)


cmap = colormap_tor([.5 0 1], [1 1 0]);


mymax = prctile(abs(C.cdata), 95);
suit_plotflatmap(C.cdata, 'cmap', cmap, 'cscale',[-mymax mymax]);

% f = gifti(giiname)
% mymax = prctile(abs(f.cdata), 95);
% suit_plotflatmap(f.cdata, 'cmap', cmap, 'cscale',[-mymax mymax]);


