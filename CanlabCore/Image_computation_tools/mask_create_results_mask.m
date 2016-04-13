function [cl_out, out_name] = mask_create_results_mask(mask, study_img, kern_size, varargin)
% Takes a mask image name (mask)
% smooths it by smooth_mm (optional; enter 0 for kern_size to avoid this)
% Resamples it to study image dimensions (study_img)
% writes output image
%
% :Usage:
% ::
%
%     [cl_out, out_name] = mask_create_results_mask(mask, study_img, kern_size, [opt: expression to evaluate on mask])
%
% :Optional: Applies an expression to be evaluated to the image, 
% in spm_imcalc format, e.g., 'i1 < .05'  (uses spm_imcalc_ui.m)
%
% :Outputs: Returns mask_clusters (mask_cl) and mask name (out_name)
%
% :Examples:
% ::
%
%    mask = which('spm2_amy.img')
%    study_img = '/Users/tor/Documents/Tor_Documents/PublishedProjects/inpress_2007_Emotion_handbook_2006/an_metaFWE_rad10/exp_vs_percept/Activation_FWE_all.img';
%    kern_size = 3;
%    cl_out = mask_create_results_mask(mask, study_img, kern_size);
%    cluster_orthviews(cl_out, {[0 1 0]}, 'add', 'handle', 1);
%
%    [cl_out, out_name] = mask_create_results_mask('X-Y_total_pvals.img', 'X-Y_total_pvals.img', 0, 'i1 < .05 & i1 > eps');
%    spm_image('init', 'X-Y_total_pvals.img');
%    cluster_orthviews(cl_out, {[1 0 0]}, 'trans', 'add');
%

    [dd, ff, ee] = fileparts(mask);

    out_name = ['mask_' ff ee];
    disp(['Creating: ' out_name]);

    % warning: spm_smooth will introduce non-zero values even with a kernel size
    % of 0. avoid by conditional statement:
    if kern_size
        spm_smooth(mask, out_name, kern_size)
        scn_map_image(out_name, study_img, 'write', out_name);

    else
        scn_map_image(mask, study_img, 'write', out_name);
    end
    

    %spm_conv_vol(maskData,maskData,2,2,2,[0 0 0]);


    % optional: threshold/eval
    if ~isempty(varargin)
        out_name = spm_imcalc_ui(out_name, out_name, varargin{1});
    end


    cl_out = mask2clusters(out_name);

end
