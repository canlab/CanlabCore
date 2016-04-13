function searchlight_dream(dat, dist_n, mask, varargin)
% This function generates codes for submitting multiple distributed jobs to 
% clusters to run a searchlight analysis on multiple chunks of the brain.
%
% :Usage:
% ::
%
%     searchlight_dream(dat, dist_n, mask, 'algorithm_name', 'cv_svm' (or 'cv_lassopcr'), 'cv_assign', whfolds, [optional input])
%
% :Features:
%   - generates dist_n scripts in modeldir (or current directory)
%   - can run a searchlight analysis on one dataset, or two datasets
%     (cross-classification)
%   - can apply different radius
%   - can obtain cross-validated results with 'cv_assign' option
%   - can use SVM (linear svm is a default) and LASSO-PCR. You need to have
%     a spider toolbox and lasso rocha toolbox in your path
%   - you can save predictive weights for each searchlight or discard them
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2014  Wani Woo
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
%
% :Inputs:
%
%   **dat:**
%        image_vector or fmri_data object with data
%
%   **dat1.Y(:,1):**
%        for svm: true(1) or false(-1) for each observation (image) in Y(:,1)
%        for lassopcr: continuous value for Y(:,1)
%
%   **dat1.Y(:,[2:n]):**
%        Test sets: could be binary: and true(1), false(-1),
%        ignore(0) or continuous values
%
%   **dist_n:**
%        The number of jobs (brain chunks) you want to create
%
%   **mask:**
%        This will be run on voxels within the mask
%        e.g., which('scalped_avg152T1_graymatter.img')
%
%   **'algorithm_name':**
%        should be followed by 'cv_svm' or 'cv_lassopcr'
%
%   **'cv_assign':**
%        a vector of integers for membership in custom holdout set 
%        of each fold
% 
% :Optional Inputs: 
%
%   **'dat2':**
%        cross-classification; should be followed by dat2 and dat2.Y
%        dat2.Y(:,1) - for training/testing, dat2(:,[2:n]) - for testing
%
%   **'r':**
%        searchlight sphere radius (in voxel) (default: r = 3 voxels)
%
%   **'modeldir':**
%        the directory where all the variables and results will be 
%        saved; should be followed by a directory location 
%        (default: the current directory)
%
%   **'scale':**
%        z-scored input data in image_vector or fmri_data object
%        (default = false)
%
%   **'balanced':**
%        use the balanced ridge option - balanced ridge value should
%        be followed. (default = false)
%
%   **'outcome_method':**
%        followed by the following options
%          - 'correlation' - "default" for for continuous measures
%       'twochoice'- "default" for binary outcome
%       'singleinterval'  - for binary outcome
%
%   **'save_weights':**
%        save weights for each searchlight (default = false)
%
%   **'email':**
%        should be followed by an email adress (default = false)
%
% :Outputs:
%
% This function will generate codes that call "searchlight_disti.m", which
% will save the following output variables.
%
%   **out.test_results:**
%        accuracy, p, and se for binary classification, and
%        correlation (pearson), p for prediction of continuous values
%        For the cross-classification, test_results will save
%        four values for each searchlight. The order of the test 
%        results are [dat1-on-dat1, dat1-on-dat2, dat2-on-dat1, 
%        dat2-on-dat2]. All results are cross-validated results.
%
%   **out.stats1 & stats2:**
%        stats1 and stats2 are similar to the outputs of predict function. 
%
% :Examples for lassopcr:
% ::
%
%    % data preparation
%    dat = fmri_data(which('brainmask.nii'));
%    dat.dat = randn(dat.volInfo.n_inmask, 30);
%
%    % setting up training values
%    dat.Y = dat.dat(111111, :)' + .3 * randn(30, 1);
%
%    % setting up testing values
%    dat.Y(:,2) = [ones(10,1); zeros(10,1); -ones(10,1)];
%    dat.Y(:,3) = dat.dat(111111, :)' + .3 * randn(30, 1);
%    mask = which('scalped_avg152T1_graymatter.img');
%    dist_n = 50;
%
%    % data for cross classification
%    dat2 = fmri_data(which('brainmask.nii'));
%    dat2.dat = randn(dat.volInfo.n_inmask, 30);
%    dat2.Y = dat2.dat(111, :)' + .3 * randn(30, 1);
%    dat2.Y(:,2) = dat.Y(:,2);
%    dat2.Y(:,3) = dat2.dat(111111, :)' + .3 * randn(30, 1);
%
%    % setting up other variables
%    r = 6;
%    modeldir = '/dreamio/home/chwo9116/Documents/searchlight_dream_test';
%    holdout_set = ones(6, 1); for i = 2:5, holdout_set = [holdout_set; i*ones(6, 1)]; end
%
%    % generate scripts in modeldir
%    searchlight_dream(dat, dist_n, mask, 'dat2', dat2, 'algorithm_name', ...
%    'cv_lassopcr', 'r', 6, 'modeldir', modeldir, 'cv_assign', holdout_set, ...
%    'save_weights', 'outcome_method', 'singleinterval');
%
% :See Also:
% searchlight_disti.m, xval_cross_classify.m, fmri_data.predict.m, 


modeldir = pwd; % where to save? 
send_email = false;
email_address = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'modeldir'
                modeldir = varargin{i+1};
                varargin{i} = []; varargin{i+1} = [];
            case 'email'
                send_email = true;
                email_address = varargin{i+1};
                varargin{i} = []; varargin{i+1} = [];
        end
    end
end

%% distribution prep
dat = apply_mask(dat, mask);
dat = trim_mask(dat);
vox_num = dat.volInfo.n_inmask;

dist_indx = distribution_prep(dist_n, vox_num);

% save basic variables
if exist(modeldir, 'dir')
    fprintf('modeldir is already existing. deleting and recreating the directory...\n');
    rmdir(modeldir, 's'); mkdir(modeldir); 
else
    mkdir(modeldir); 
end

savename = fullfile(modeldir, 'searchlight_dream_variables.mat');
save(savename, '-v7.3', 'dat', 'mask', 'dist_indx', 'modeldir', 'varargin');

for i = 1:dist_n
    generate_scripts(savename, modeldir, i, send_email, email_address, varargin);
end

fprintf('Scripts have been generated in %s.\n', modeldir);

end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Prep distributing jobs: Get dist_indx

function dist_indx = distribution_prep(dist_n, vox_num)

% preparation of distribution indx

unit_num = ceil(vox_num/dist_n);
unit = true(unit_num,1);

for i = 1:dist_n
    start_point = (unit_num.*(i-1)+1);
    end_point = min(start_point+unit_num-1, vox_num);
    dist_indx(start_point:end_point,i) = unit(1:end_point-start_point+1);
end

dist_indx = logical(dist_indx);

end

% -------------------------------------------------------------------------
% save scripts in modeldir

function generate_scripts(savename, modeldir, i, send_email, email_address, varargin)

% write the codes
postfix = date;
postfix(strfind(postfix, '-')) = [];
distrib_script = fullfile(modeldir, ['searchlight_dream_variables_' num2str(i) '_' postfix '.m']);

% if overwrite
%     prev_script = fullfile(modeldir, ['searchlight_dream_variables_*.m']);
%     delete(prev_script);
% end

FID = fopen(distrib_script, 'w');
    
% description
fprintf(FID, ['%% searchlight analysis: variables are saved in ' savename '\n']);
fprintf(FID, '\n');
   
fprintf(FID, '%% path definition \n');
fprintf(FID, 'reposdir = ''/dreamio3/wagerlab/Repository'';\n');
fprintf(FID, 'addpath(genpath(reposdir));\n');
fprintf(FID, 'spmdir = ''/usr/local/spm/spm8/'';\n');
fprintf(FID, 'addpath(genpath(spmdir));\n');
fprintf(FID, 'lassodir = ''/dreamio/home/chwo9116/codes_ineed/MATLAB/machine_learning/lasso_rocha'';\n');
fprintf(FID, 'addpath(genpath(lassodir));\n');
fprintf(FID, 'spiderdir = ''/dreamio/home/chwo9116/codes_ineed/MATLAB/spider'';\n');
fprintf(FID, 'addpath(genpath(spiderdir));\n');
fprintf(FID, 'mycodedir = ''/dreamio/home/chwo9116/codes_ineed/MATLAB/mycodes'';\n');
fprintf(FID, 'addpath(mycodedir);\n');
fprintf(FID, '\n');
fprintf(FID, 'use_spider;\n');
fprintf(FID, '\n');

% load variables
fprintf(FID, '%% load variables \n');
fprintf(FID, ['load(''' savename ''');\n']);

% run the analysis
fprintf(FID, '%% running searchlight \n');
fprintf(FID, '\n');
fprintf(FID, ['out = searchlight_disti(dat, mask, dist_indx(:,' num2str(i) '), varargin);\n']);

fprintf(FID, 'fout = fields(out);\n');
fprintf(FID, 'save_vars = '''';\n');
fprintf(FID, 'for i = 1:numel(fout)\n');
fprintf(FID, '\teval([fout{i} ''= out.'' fout{i} '';'']);\n');
fprintf(FID, '\tsave_vars = [save_vars '''''''' char(fout{i}) '''''', ''];\n');
fprintf(FID, 'end\n');
fprintf(FID, 'clear out;\n');
fprintf(FID, '\n');

% save data
fprintf(FID, '%% save output\n');
fprintf(FID, ['save_newname = fullfile(''' modeldir '/searchlight_results_' num2str(i) '_' postfix '.mat'');\n']);
fprintf(FID, 'eval([''save(save_newname, '' save_vars ''''''-v7.3'''');'']);\n');
fprintf(FID, '\n');
if send_email
    fprintf(FID, ['content = ''' modeldir '/searchlight_results_' num2str(i) '_' postfix ''';\n']);
    fprintf(FID, ['sendmail_wani(''' email_address ''', ''searchlight_job_done'', content);\n']);
end
fclose(FID);

end
