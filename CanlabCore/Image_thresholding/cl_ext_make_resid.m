function cl_ext_make_resid(con_files, varargin)
% This function will create residual images (4d) and mask image in 
% a current or assigned directory in order to use them in estimating smoothness 
% (relevant functions: spm_est_smoothness (SPM), 3dFWHMx (AFNI), smoothest (FSL). 
%
% :Usage:
% ::
%
%     function cl_ext_make_resid(conimgs, varargin)
%
% :Inputs:
%
%   **con_files:**
%        contrast image file names; This could be a cell array or
%        strings. This could be 4d images. 
%
%        Best: Input a cell string. e.g., for a string matrix:
%
%        Use cl_ext_make_resid(cellstr(imgs)); % save residual images
%
%   * If you are not providing the absolute paths of the images, you need to
%   be in the directory that has the image files. 
%
% :Outputs: 
%
%   **Res4d.nii:**
%        residual images saved by SPM. 
%
%   **mask.nii:**
%        the mask image that was used.
%
% :Options for varargin:
%
%   **'mask'**
%        This option can be used to estimate a cluster size for the correction for multiple
%        comparisons "within the mask". You can put in a ROI mask or gray matter,
%        whatever. If you don't specify a mask image, brainmask.nii (default) will be
%        used, but the image has to be in your path.
%         e.g.)
%         ::
%
%              mask = fullfile(basedir, 'ROI_image.img'); 
%              mask = which('scalped_avg152T1_graymatter_smoothed.img'); % limited to gray matter
%
%   **'outputdir'**
%        With this option, this will save residual and mask images and in the 
%        outputdir directory. If you don't give outputdir, the current directory
%        will be used (default). 
%
% This function calls cl_ext_spm_spm.m, which is a modified spm_spm not to
% delete residual images. 
%
% ..
%    Choong-Wan (Wani) Woo, 01/22/2013
%    modified by Wani, 05/18/2013
% ..

outputdir = pwd;
mask = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'mask', mask = varargin{i+1};
            case 'outputdir', outputdir = varargin{i+1}; 
        end
    end
end

temp_dir = fullfile(outputdir, 'cl_extent_estimation');
mkdir(temp_dir);

if isempty(mask)
    try
        mask = which('brainmask.nii');
    catch
        warning('We tried to use a default mask (brainmask.nii), but we couldn''t use it. Maybe you need to add "brainmask.nii" into your path. We are going to use the first contrast image as a mask.');
        mask = conimgs{1};
    end
end

matlabbatch{1}.spm.stats.factorial_design.dir{1} = temp_dir; % result dir

% file names should be in cell arrays. 

if ~iscell(con_files)
    if size(con_files, 1) == 1 
        con_files = filenames(con_files, 'absolute', 'char');
        con_files = expand_4d_filenames(con_files);
        for ii = 1:size(con_files,1), conimgs{ii} = deblank(con_files(ii,:)); end
    else
        for ii = 1:size(con_files,1)
            conimgs(ii) = filenames(deblank(con_files(ii,:)), 'absolute'); % wani added this line for efficiency
%             con_files_temp(ii,:) = filenames(con_files(ii,:), 'absolute', 'char'); 
        end
%         con_files = con_files_temp;
    end
%     for i = 1:size(con_files,1)
%         conimgs{i} = deblank(con_files(i,:));
%     end
    conimgs = conimgs';
else
    if length(con_files) == 1 
        con_files = expand_4d_filenames(con_files);
        for i = 1:size(con_files,1)
            conimgs{i} = deblank(con_files(i,:));
        end
        conimgs = conimgs';
    else
        conimgs = con_files;
    end        
end

con_num = length(conimgs);
if con_num == 1
    error('The number of contrast images is one. Please check the image names.');
end

matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = conimgs; % each con image for 2nd level
% matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', [], 'cname',[],'iCFI',[],'iCC',[]);

matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;

matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = mask;

matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

cd(temp_dir);
save spm_2nd matlabbatch;

ndf = [1 con_num-1];

if exist('SPM.mat', 'file'), delete('SPM.mat'); end
spm_jobman('run', 'spm_2nd.mat');

load SPM.mat;

try
    cl_ext_spm_spm(SPM); % you need cl_ext_spm_spm.m in your path
catch
    error('You need "cl_ext_spm_spm.m" in your path.');
end

mask_img = filenames('mask*.img','char');
Res_imgs = filenames('ResI_*.img', 'char');

dat_mask = fmri_data(mask_img); 
dat_mask.fullpath = fullfile(outputdir, 'mask.nii');
write(dat_mask);

dat_res = fmri_data(Res_imgs);
dat_res.fullpath = fullfile(outputdir, 'Res4d.nii');
write(dat_res);

% rmdir(temp_dir, 's');

cd(outputdir);

end
