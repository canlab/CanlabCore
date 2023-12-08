
% Diedrichsen, J. & Zotow, E. (2015). Surface-based display of volume-averaged cerebellar data. PLoS One, 7, e0133402.
% https://www.diedrichsenlab.org/imaging/suit_flatmap.htm

% g = gifti('FLAT.surf.gii')
% figure; han = patch('Faces',g.faces,'Vertices',g.vertices);
% set(han, 'EdgeColor', 'none');

targetniifile = '/home/bogdan/Downloads/full_pla_g_unthresh.nii';

targetniifile = fullfile(pwd, 'nps.nii');

%%

targetniifile = which(targetniifile);       % get full path
[dd, ff, ee] = fileparts(targetniifile);

wdtargetniifile = fullfile(dd, ['wd' ff ee]);

giiname = fullfile(dd, [ff '.func.gii']);

affinename = which('Affine_MNI152NLin6Asym_T1_1mm_seg1.mat');

flowname = which('u_a_MNI152NLin6Asym_T1_1mm_seg1.nii');

cerebmaskfile = which('c_MNI152NLin6Asym_T1_1mm_pcereb.nii');

if isempty(which('spm_dartel_norm.m')), error('Start SPM to find spm_dartel_norm.m'), end

% Normalize template to SUIT cerebellar template - define mapping to SUIT space
% suit_normalize_dartel(struct('subjND', struct('gray', {{'MNI152NLin6Asym_T1_1mm_seg1.nii'}}, 'white', {{'MNI152NLin6Asym_T1_1mm_seg2.nii'}}, 'isolation', {{'c_MNI152NLin6Asym_T1_1mm_pcereb.nii'}})))

% Apply normalization : Reslice the target image, already in MNI space, in SUIT space

% See also:
% apply_spm_warp(mvg_img0, fxd_img0, pre_affine_mat, warp_img, post_affine_mat, out_img, interp)

suit_reslice_dartel(struct('subj', struct('affineTr', {{affinename}}, ...
    'flowfield', {{flowname}}, ...
    'resample', {{targetniifile}}, ...
    'mask', {{cerebmaskfile}})))

C = [];
C.cdata = suit_map2surf(wdtargetniifile, 'space','SUIT', 'stats', @nanmean);

save(gifti(C), giiname)

% f = gifti(giiname)

cmap = colormap_tor([.5 0 1], [1 1 0]);

mymax = prctile(abs(C.cdata), 95);
suit_plotflatmap(C.cdata, 'cmap', cmap, 'cscale',[-mymax mymax]);


mymax = prctile(abs(f.cdata), 95);
suit_plotflatmap(f.cdata, 'cmap', cmap, 'cscale',[-mymax mymax]);


