function [cl_ext_ClustSim, fwhm_vox] = cl_ext_3dClustSim(corrected_p, prim_p, residual_images, mask, voxelsize_mm, ClustSim_dir, varargin)
% :Usage:
% ::
%
%     [cl_ext_ClustSim, fwhm] = cl_ext_3dClustSim(corrected_p, prim_p, residual_images, mask, voxelsize_mm, ClustSim_dir, varargin)
% 
% :Inputs:
%
%   **corrected_p:**
%        cluster-extent corrected p value
%
%        e.g.) if cluster-extent corrected p < .05: corrected_p = .05
%
%   **prim_p:**
%        primary threshold for height (i.e., cluster-defining threshold)
%
%        e.g.) prim_p = [0.01 0.005 0.001];
%
%   **residual_images:**
%        residual image names; if you used
%
%        cl_ext_make_resid.m, this should be 'Res4d.nii'. 
%
%        e.g.) residual_images = filenames('Res4d.hdr', 'char', 'absolute');
%
 %             residual_images = filenames('Res4d.nii', 'char', 'absolute');   
%
%   **mask:**
%        mask image name (should have header)
%
%        e.g.) mask = filenames('mask.hdr', 'char', 'absolute');
%
%              mask = filenames('mask.nii', 'char', 'absolute');
%
%   **voxelsize_mm:**
%        voxel sizes in milimeter. e.g) voxelsize_mm = [2 2 2];
%
%   **3dClustSim_dir:**
%        directory where alphasim is installed.
%
%        e.g.) 3dClustSim_dir = '/Users/clinpsywoo/abin/macosx_10.6_Intel_64';
%
%        If you don't have 3dClustSim, see http://afni.nimh.nih.gov/pub/dist/HOWTO/howto/ht00_inst/html/index.shtml
%
% :Optional Inputs:
%
%   **'iter':**
%        you can set up the iteration number for Monte Carlo simulation.
%        default is doing 1000 iterations.
%
%   **'twotail':**
%        default is one-tail - with this option, primary_p/2 will be used 
%        for all clsuter extent estimations. 
%
%   **'fwhm':**
%        you can add fwhm manually
%
% :Outputs: 
%
%   **cl_ext_ClustSim:**
%        cl_ext_ClustSim is the cluster size that makes a corrected p value under 
%        corrected_p (e.g., 0.05). 
%
%   **fwhm (x, y, z in voxel):**
%        intrinsic smoothness level estimated by AFNI(3dFWHMx). 
%        If you want to convert this into mm, you need to multiply these
%        values by voxel sizes in mm. 
%
% ..
%    Choong-Wan (Wani) Woo, 01/21/2013
%    modified by Wani, 05/18/2013
% ..


curr_dir = pwd; %% go to the alphasim directory
cd(ClustSim_dir);

%% defaults
iter = 1000;
IsTwoTailed = 0;
manual_fwhm = 0;
usemask = 1;

%% get options
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'iter', iter = varargin{i+1};
            case {'twotail', 'twotails', 'twotailed'}, IsTwoTailed = 1;
            case 'fwhm', fwhm = varargin{i+1}; manual_fwhm = 1;
        end
    end
end

%% estimate smoothness 
if ~manual_fwhm
    if usemask
        eval_function_fwhm = ['unset DYLD_LIBRARY_PATH; ./3dFWHMx -mask ' mask ' -detrend -dset ' residual_images];
    else
        eval_function_fwhm = ['unset DYLD_LIBRARY_PATH; ./3dFWHMx -detrend -dset ' residual_images];
    end
    [status, res_fwhm] = unix(eval_function_fwhm);
    
    k = textscan(res_fwhm, '%s'); k = k{1};
    b = []; for i = 1:length(k), b = cat(1,b,str2num(k{i})); end
    
    fwhm = b(end-2:end);
end

fwhm = fwhm';
fwhm_vox = fwhm./voxelsize_mm;
% 06/27/13 Wani: 3dFWHMx gives fwhm in mm, not vox.  

% if sum(size(voxelsize_mm) == size(fwhm)) > 0
%     fwhm_mm = voxelsize_mm.*fwhm;
% else
%     fwhm = fwhm';
%     fwhm_mm = voxelsize_mm.*fwhm;
% end

%% calculate cluster extent size

if IsTwoTailed
    eval_function_cl_ext = ['unset DYLD_LIBRARY_PATH; ./3dClustSim -mask ' mask ' -dxyz ' num2str(voxelsize_mm) ...
        ' -iter ' num2str(iter) ' -pthr ' num2str(prim_p/2) ' -fwhmxyz ' num2str(fwhm) ' -athr ' num2str(corrected_p)];
else
    eval_function_cl_ext = ['unset DYLD_LIBRARY_PATH; ./3dClustSim -mask ' mask ' -dxyz ' num2str(voxelsize_mm) ...
        ' -iter ' num2str(iter) ' -pthr ' num2str(prim_p) ' -fwhmxyz ' num2str(fwhm) ' -athr ' num2str(corrected_p)];
end
    
[status, res_ext] = unix(eval_function_cl_ext);

clear k b;
k = textscan(res_ext, '%s'); k = k{1};
b = []; for i = 1:length(k), if isnumeric(str2num(k{i})), b = cat(1,b,str2num(k{i})); end, end

for i = 1:length(prim_p)
    if IsTwoTailed
        ii = find(b == prim_p(i)/2);
    else
        ii = find(b == prim_p(i));
    end
    
    if numel(ii) == 2
        prim_p_idx(i) = ii(2);
    elseif numel(ii) == 1
        prim_p_idx(i) = NaN;
    elseif numel(ii) > 3
        if i > 1
            if sum(ii == prim_p_idx(i-1)+2) == 0, prim_p_idx(i) = NaN;
            else prim_p_idx(i) = ii(ii == prim_p_idx(i-1)+2);
            end
        else
            prim_p_idx(i) = NaN;
        end
    end
end

cl_ext_ClustSim = zeros(size(prim_p))';
cl_ext_ClustSim(isnan(prim_p_idx)) = NaN;

cl_ext_ClustSim(~isnan(prim_p_idx)) = round(b(prim_p_idx(~isnan(prim_p_idx))+1));

cd(curr_dir);

return
