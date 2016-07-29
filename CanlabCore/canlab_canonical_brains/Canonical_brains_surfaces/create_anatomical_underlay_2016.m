% paths

p = '/Users/tor/Google_Drive/CanlabDataRepository/ROI_masks_and_parcellations/Anatomical_Rendering';
addpath(p)

p = '/Users/tor/Google_Drive/CanlabDataRepository/ROI_masks_and_parcellations/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo';
addpath(genpath(p))

p = '/Users/tor/Google_Drive/CanlabDataRepository/ROI_masks_and_parcellations/Subcortex_Keuken2014';
addpath(genpath(p))


single_subj = 'SPM8_colin27T1_seg.img';
icbm = 'icbm152_2009_symmetric_for_underlay.img';

dartel1 = 'dartel_spm12_mni152_gray_matter_mask.img';
dartel2 = 'dartel_spm12_mni152_brainmask.img';
dartel3 = which('dartel_spm12_icbm152.nii');

glasser1 = 'clean_Q1-Q6_RelatedParcellation210_AverageT1w_restore.nii';
glasser2 = 'Q1-Q6_RelatedParcellation210_AverageT1w_restore.nii';

keuken = which('MNI152_T1_04mm_brain.nii');

%% COMPARE

imgs = { icbm dartel2 dartel3 glasser1};


n = length(imgs);

shiftvals = [0:.17:n]; % more than we need, but works

create_figure('BRAINS!!!'); 
axis off

for i = 1:n
    
    o2 = fmridisplay('overlay', which(imgs{i}));
    
    
    % saggital
    axh = axes('Position', [-0.02 0.15+shiftvals(i) .17 .17]);
    o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'onerow', 'noverbose', 'existing_axes', axh);
    drawnow
    
    % axial
    axh = axes('Position', [.015 0.15+shiftvals(i) .17 .17]);
    o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
    drawnow
    
end


%% CLEAN UP MASK
% take icbm09, contrast enhance based on dartel3, mask with dartel2

dat = fmri_data(which(icbm));
orthviews(dat)

mask = resample_space(fmri_data(which(dartel2)), dat);
con = resample_space(fmri_data(which(dartel3)), dat);

mask = get_wh_image(mask, 1);

dat_orig = dat;

%%
dat = dat_orig;

gray = con.dat(:, 1) > .5;
white = con.dat(:, 2) > .5;
edge = con.dat(:, 1) < .1 & con.dat(:, 2) < .1 ;

csf = dat.dat < 40 & con.dat(:, 1) < .2 & con.dat(:, 2) < .2 ;
% to-do midline cortex: exclude 

% brainstem mask separate: add to gray
% striatal gray mask: add to gray (keuken 2014)

grayprop = con.dat(:, 1);

%dat.dat(gray) = dat.dat(gray) .* 1.2;
dat.dat(gray) = dat.dat(gray) .* ((1 - grayprop(gray)) + .6);

dat.dat(white) = dat.dat(white) .* 1;
dat.dat(edge) = dat.dat(edge) .* 0;
dat.dat(csf) = dat.dat(csf) .* 0;

dat.dat(mask.dat < .99) = 0;

orthviews(dat);

d = fileparts(dat.fullpath);
dat.fullpath = fullfile(d, 'icbm152_2009_symm_enhanced_for_underlay.img');
write(dat);

spm_image('init', dat.fullpath);

o2 = canlab_results_fmridisplay([], 'compact2', 'overlay', dat.fullpath);


