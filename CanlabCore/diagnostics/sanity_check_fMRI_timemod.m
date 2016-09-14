function fig = sanity_check_fMRI_timemod(imgs, SPM, roi)
% Sanity check your time modeling is accurate
% Plots average signal from a specified sensory ROI (e.g., visual cortex)
%  with your design file over it. Visually inspect that increased in signal
%  are time locked to stimulus presentation.
%  For example, if you are presenting images, you'd expect avg visual
%  cortex signal to increase after each image presentation.
%  User specifies sensory cortex.
%  SPM anatomy toolbox must be on path to
%  use the flags: "visual" "motor" "auditory" or the user can input a path
%  to an ROI of their choosing.
%  This code takes an SPM.mat file. Single subj only.
%  Don't forget to give all images if regressors are concatenated.
%  Requires: CANLabCore, SPM, SPM Anatomy toolbox
%
%  Note: This is for a cursory visual sanity check only!
%
%  by Marianne, 2016
%
% :Usage:
% ::
%     sanity_check_fMRI_timemod(imgs, SPM, roi)
%
% :Examples:
% ::
%    imgs={'/r1IAPS/swdIAPS_pa__32ch_mb8_v01_r01_0005_AQ.nii';'/r2IAPS/swdIAPS_pa_32ch_mb8_v01_r02_0007_AQ.nii'};
%    load('SPM.mat');roi='visual';
%    sanity_check_fMRI_timemod(imgs, SPM, roi)
%

%%
% set up test images
dat = fmri_data(imgs);

%%
% set up ROI Mask

if ischar(roi)
    switch roi
        case 'visual'
            roi = fmri_mask_image(which('Visual_hOc1.nii')); % add warning SPM Anat toolbox not on path
        case 'motor'
            roi = fmri_mask_image(which('Motor_4p.nii'));
        case 'auditory'
            roi = fmri_mask_image(which('Auditory_Te3.nii'));
        otherwise
            roi = fmri_mask_image(roi);
    end
else
    error('Please specify an ROI whose signal you want to check.')
end

%%
% set up DSGN file
% give it an SPM dsgn file
wh_cols = scn_spm_get_events_of_interest(SPM, 'events_only');
reg=(SPM.xX.X(:, wh_cols));
event_count=0;
for c=1:size(reg,2)
    [p l]=findpeaks(reg(:,c));
    reg(l,c)=100;
    event_count=event_count+length(l);
end
reg(find(reg<100))=-3;
reg(find(reg==100))=3;
figure(1);imagesc(reg);title(sprintf('%d Total Events; %d Regressor Onsets',size(reg,2),event_count),'FontSize',15)
disp(sprintf('There are %d regressors and %d total events counted.',size(reg,2),event_count));
saveas(gcf,'RegressorCheck.png')
%%
% Set up plot

pexp = apply_mask(dat,roi,'pattern_expression','correlation','ignore_missing');
pexp=zscore(pexp);
if size(pexp,1) ~= size(reg,1)
    error('Your regressor is not the same size as your image. Did you concatenate in the SPM model but only input one runs image?')
end

fig=figure(2); 
set(fig, 'Position', [10 900 1200 200])
fig=plot(pexp);hold on;
for i = 1:size(reg,2)
    fig=plot(reg(:,i),'Color',rand(1,3));hold on;%[.05*i .05*i .05*i]
end
title('Avg ROI signal time course','FontSize', 22);
saveas(fig,'ROItimecourseCheck.png')

disp('User should zoom in and visually examine if the signal is increased post stimulation.')
end
