% The two brains that this script creates are:

% ----------------------------------------------------------------------
% icbm152_2009_symm_enhanced_for_underlay.img
% keuken_2014_enhanced_for_underlay.img

% Another good underlay may be:
% 'Q1-Q6_RelatedParcellation210_AverageT2w_restore.nii'
% ----------------------------------------------------------------------

% paths

p = '/Users/tor/Google_Drive/CanlabDataRepository/ROI_masks_and_parcellations/Anatomical_Rendering';
addpath(p)

p = '/Users/tor/Google_Drive/CanlabDataRepository/ROI_masks_and_parcellations/Glasser_et_al_2016_HCP_MMP1.0_RVVG/HCP_PhaseTwo';
addpath(genpath(p))

p = '/Users/tor/Google_Drive/CanlabDataRepository/ROI_masks_and_parcellations/Subcortex_Keuken2014';
addpath(genpath(p))


single_subj = 'SPM8_colin27T1_seg.img';

% the icbm is a good one, symmetric
icbm = 'icbm152_2009_symmetric_for_underlay.img';

dartel1 = 'dartel_spm12_mni152_gray_matter_mask.img';
dartel2 = 'dartel_spm12_mni152_brainmask.img';
dartel3 = which('dartel_spm12_icbm152.nii');

glasser1 = 'clean_Q1-Q6_RelatedParcellation210_AverageT1w_restore.nii';
glasser2 = 'Q1-Q6_RelatedParcellation210_AverageT1w_restore.nii';

keuken = which('MNI152_T1_04mm_brain.nii');

ribbon = 'SPM8_colin27T1_cortical_ribbon.img';

% this is a good one
t2 = which('Q1-Q6_RelatedParcellation210_AverageT2w_restore.nii')

% This is the underlay we are creating
underlay = 'icbm152_2009_symm_enhanced_for_underlay.img';

% load regions for known gray areas
load Keuken_2014_7T_regions

% Morel?
thalold = load('thal_cl.mat');

rvm = fmri_data(which('RVMmask_symm_2mm.nii'));

thal = load('morel_thalamus_atlas_whole.mat');

%% COMPARE

imgs = { icbm underlay };


n = length(imgs);

shiftvals = [0:.17:n]; % more than we need, but works

create_figure('BRAINS!!!'); 
set(gcf, 'Position', [105         444        1188         274])
enlarge_axes(gcf, 1.4)
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

%% SHOW REGIONS TO ENHANCE
% add them to a "Known Gray" vector

% Thalamus, PAG, LC

% Keuken 2014 : Striatum at 7T
o2 = removeblobs(o2);
o2 = addblobs(o2, striatum, 'color', [1 0 0]);
o2 = addblobs(o2, GPe, 'color', [0 1 0]);
o2 = addblobs(o2, GPi, 'color', [0 1 1]);
o2 = addblobs(o2, SN, 'color', [1 .5 0]); % too small?
o2 = addblobs(o2, STN, 'color', [1 0 1]); % too small?
o2 = addblobs(o2, RN, 'color', [1 0 0]);

dat = fmri_data(which(icbm));

known_gray = region2imagevec(striatum, dat);

tmp = region2imagevec(GPe, dat);
known_gray.dat = known_gray.dat + tmp.dat;  % add this to the known gray regions

tmp = region2imagevec(GPi, dat);
known_gray.dat = known_gray.dat + tmp.dat;  % add this to the known gray regions

tmp = region2imagevec(RN, dat);
known_gray.dat = known_gray.dat + tmp.dat;  % add this to the known gray regions



%% CLEAN UP MASK
% take icbm09, contrast enhance based on dartel3, mask with dartel2

dat = fmri_data(which(icbm));
orthviews(dat)

mask = resample_space(fmri_data(which(dartel2)), dat);
con = resample_space(fmri_data(which(dartel3)), dat);

mask = get_wh_image(mask, 1);

dat_orig = dat;
con_orig = con;

%%
dat = dat_orig;
con = con_orig;

gray = con.dat(:, 1) > .5;
white = con.dat(:, 2) > .5;
edge = con.dat(:, 1) < .1 & con.dat(:, 2) < .1 ;

csf = dat.dat < 40 & con.dat(:, 1) < .2 & con.dat(:, 2) < .2 ;

known_gray_vox = known_gray.dat > 0;

% to-do midline cortex: exclude 

% brainstem mask separate: add to gray

% striatal gray mask: add to gray (keuken 2014)
gray(known_gray_vox) = true;
con.dat(known_gray_vox, 1) = prctile(con.dat(gray, 1), 80);% con.dat(known_gray_vox, 1) + .2; % prctile(con.dat(gray, 1), 80);

grayprop = con.dat(:, 1);

%dat.dat(gray) = dat.dat(gray) .* 1.2;
dat.dat(gray) = dat.dat(gray) .* ((1 - grayprop(gray)) + .6);
dat.dat(white) = dat.dat(white) .* 1;
dat.dat(edge) = dat.dat(edge) .* 0;
dat.dat(csf) = dat.dat(csf) .* 0;

%dat.dat(known_gray_vox) = mean(dat.dat(gray));

dat.dat(mask.dat < .99) = 0;

orthviews(dat);

d = fileparts(dat.fullpath);
dat.fullpath = fullfile(d, 'icbm152_2009_symm_enhanced_for_underlay.img');
write(dat);

spm_image('init', dat.fullpath);

%o2 = canlab_results_fmridisplay([], 'compact2', 'overlay', dat.fullpath);
o2 = fmridisplay('overlay', which(underlay));
o2 = montage(o2, 'slice_range', [-10 10]);


% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%% KEUKEN BRAIN CLEAN UP MASK
% Keuken brain - clean up/enhance contrast
% contrast enhance based on dartel3, mask with dartel2

% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

t2dat = fmri_data(t2);

dat = fmri_data(which(keuken));

dat = resample_space(dat, t2dat);  % downsample to 0.7 mm space
orthviews(dat)

mask = resample_space(fmri_data(which(dartel2)), dat);
con = resample_space(fmri_data(which(dartel3)), dat);

mask = get_wh_image(mask, 1);

dat_orig = dat;
con_orig = con;

%%

dat = dat_orig;
con = con_orig;

% gray = con.dat(:, 1) > .5;
% white = con.dat(:, 2) > .5;
edge = con.dat(:, 1) < .1 & con.dat(:, 2) < .1 ;

csf = dat.dat < 40 & con.dat(:, 1) < .2 & con.dat(:, 2) < .2 ;

%known_gray_vox = known_gray.dat > 0;

% striatal gray mask: add to gray (keuken 2014)
gray(known_gray_vox) = true;
%con.dat(known_gray_vox, 1) = prctile(con.dat(gray, 1), 80);% con.dat(known_gray_vox, 1) + .2; % prctile(con.dat(gray, 1), 80);

grayprop = con.dat(:, 1);

%dat.dat(gray) = dat.dat(gray) .* 1.2;
%     dat.dat(gray) = dat.dat(gray) .* ((1 - grayprop(gray)) + .6);
%     dat.dat(white) = dat.dat(white) .* 1;
dat.dat(edge) = dat.dat(edge) .* 0;
dat.dat(csf) = dat.dat(csf) .* 0;

%dat.dat(known_gray_vox) = mean(dat.dat(gray));

dat.dat(mask.dat < .99) = 0;

orthviews(dat);

d = fileparts(dat.fullpath);
dat.fullpath = fullfile(d, 'keuken_2014_enhanced_for_underlay.img');
write(dat);

spm_image('init', dat.fullpath);

%o2 = canlab_results_fmridisplay([], 'compact2', 'overlay', dat.fullpath);
o2 = fmridisplay('overlay', which(dat.fullpath));
o2 = montage(o2, 'slice_range', [-10 10]);

CMAP = contrast(dat.dat)
colormap(CMAP)
