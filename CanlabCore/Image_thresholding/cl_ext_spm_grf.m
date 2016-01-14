function [cl_ext_spm, fwhm] = cl_ext_spm_grf(corrected_p, prim_p, residual_images, mask, varargin) 
% This function is designed to estimate a cluster extent size for 
% the correction for multiple comparisons based on a Gaussian Random Field 
% Theory using SPM toolboxes. 
%
% :Usage:
% ::
%
%     [cl_ext, fwhm] = cl_ext_spm_grf(corrected_p, prim_p, residual_images, mask, varargin) 
%
% :Inputs:
%
%   **corrected_p:**
%        corrected p value
%
%        e.g.) cluster-extent corrected p < .05: corrected_p = .05
%
%   **prim_p:**
%        primary threshold for height (i.e., cluster-defining threshold)
%
%        e.g.) prim_p = [0.01 0.005 0.001];
%
%   **residual_images:**
%        residual image names; if you used
%
%        cl_ext_make_resid.m, this should be 'Res4d.nii'
%
%   **mask:**
%        mask image name
%
% :Optional Inputs:
%
%   **'doplot':**
%
%   **'twotail':**
%        default is one-tail - with this option, primary_p/2 will be used 
%        for all clsuter extent estimations. 
%
% Output:
%
%   **cl_ext_spm:**
%        cl_ext_spm is the cluster size that makes a corrected p value under 
%       corrected_p (e.g., 0.05). 
%
%   **fwhm (x, y, z in voxels):**
%        intrinsic smoothness level estimated by SPM (spm_est_smoothness.m)
%       If you want to convert this into mm, you need to multiply these
%       values by voxel sizes in mm. 
%
% ..
%    Choong-Wan (Wani) Woo, 08/13/2012
%    modified by Wani, 05/18/2013
% ..
 
doplot = false;
isTwoTailed = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'doplot', doplot = true;
            case 'twotail', isTwoTailed = true;
        end
    end
end

con_num = size(expand_4d_filenames(residual_images),1);
ndf = [1 con_num-1];

[fwhm, dummy, r] = spm_est_smoothness(residual_images, mask, [con_num con_num-1]);

V2R = 1/prod(fwhm);

cl_ext = zeros(length(prim_p),2);

if doplot
    scrsz = get(0, 'ScreenSize');
    create_figure('cluster_extent_spm');
    set(gcf, 'Position', [1 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/2]);
end

for i = 1:length(prim_p)
    if isTwoTailed
        u(i) = spm_u(prim_p(i)/2, ndf, 'T');
    else
        u(i) = spm_u(prim_p(i), ndf, 'T');
    end
end

for i = 1:length(prim_p) % cluster-extent threshold (k) based on RF in SPM

    P = [];
        
    for j = 1:100000
        k = j*V2R; % convert the number of voxels into the number of resels
        
        P(end+1,1) = spm_P_RF(1,k,u(i),ndf,'T',r,1);
        P(end,2) = j;

        if P(end,1) <= corrected_p
            cl_ext(i,1) = j;
            cl_ext(i,2) = P(end,1);
            break
        end 
    end

    if doplot
        hh(i) = subplot(1,length(prim_p),i);
        plot(P(:,2), P(:,1), '-b', 'LineWidth', 1.5);
        xlabel('cluster extent size', 'FontSize', 16);
        ylabel('corrected P value', 'FontSize', 16);
        eval(['title(hh(i), ''Primary P:' num2str(prim_p(i)), ''', ''fontsize'', 16);']);
        set(gca, 'FontSize', 15);
        hold on;
        plot(cl_ext(i,1), cl_ext(i,2), 'r+', 'MarkerSize', 15);
        eval(['text(cl_ext(i,1)/2, cl_ext(i,2), ''+ cl size:' num2str(cl_ext(i,1)) ''', ''FontSize'', 14);']);
    end

end

cl_ext_spm = cl_ext(:,1);

return

