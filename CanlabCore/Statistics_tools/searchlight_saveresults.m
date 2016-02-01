function searchlight_saveresults(modeldir)
% This function combines and saves searchlight analysis results. 
%
% :Usage:
% ::
%
%     searchlight_saveresults(modeldir)
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
% :Input:
%
%   **modeldir:**
%        the directory where the searchlight result mat files (e.g.,
%        searchlight_results_*_01Aug2014.mat)
% 
% :Outputs:
%
%   This function saves the brain maps where each voxel contains a summary 
%   stat value (e.g., accuracy or outcome correlation) in the modeldir. The
%   naming convention of the result maps are as follow:
%
%   results_searchlight_(a)_(b)_dat(c).nii
%      1. will be a number - the number of test results
%      2. acc: accuracy, r: outcome correlation, p: p-values, se: standard error, 
%         and thr: threshold for the single-interval test
%      3. If the test was on one dataset, (c) will be empty, but if the
%         test was cross-prediction, (c) will be 11, 12, 21, or 22. 
%           - 11: trained on the first dataset, and tested on the first dataset
%           - 12: trained on the first dataset, and tested on the second dataset
%           - 21: trained on the second dataset, and tested on the first dataset
%           - 22: trained on the second dataset, and tested on the second dataset
% 
% :Output Examples:
% ::
%
%    results_searchlight_1_r_dat11.nii (outcome correlation)
%    results_searchlight_1_p_dat11.nii (p value for the correlation values)
%
%    results_searchlight_2_acc_dat11.nii (accuracy)
%    results_searchlight_2_p_dat11.nii (p value for the accuracy)
%    results_searchlight_2_se_dat11.nii (standard error for the accuracy)
%    results_searchlight_2_thr_dat11.nii (threshold for the accuracy test)


load(fullfile(modeldir, 'searchlight_dream_variables.mat'));
data = load(fullfile(modeldir, 'searchlight_dream_variables.mat'), 'varargin');

do_cross = false;

% parsing varargin to know it's cross_classify or not
for i = 1:numel(data.varargin)
    if ischar(data.varargin{i})
        switch data.varargin{i}
            case 'dat2'
                do_cross = true;
        end
    end
end

res_files = filenames(fullfile(modeldir, 'searchlight_results_*mat'));

% load one test_results to set up results variables
load(res_files{1}, 'test_results');

% image frame
if do_cross
    dat.dat = zeros(size(dat.dat,1),4);
else
    dat.dat = zeros(size(dat.dat,1),1);
end

% results fields
for m = 1:numel(test_results)
    res_fields{m} = fields(test_results{m});
    for j = 1:numel(res_fields{m}), eval(['dat_' num2str(m) '_'  res_fields{m}{j} ' = dat;']); end
end

% combining results
for ii = 1:numel(res_files)
    load(res_files{ii}, 'test_results');
    
    for m = 1:numel(test_results)
        for j = 1:numel(res_fields{m})
            eval(['test_results{' num2str(m) '}.' res_fields{m}{j} '(isnan(test_results{' num2str(m) '}.' res_fields{m}{j} ')) = 0;']);
            eval(['dat_' num2str(m) '_' res_fields{m}{j} '.dat = dat_' num2str(m) '_' res_fields{m}{j} '.dat + test_results{' num2str(m) '}.' res_fields{m}{j} ';']);
        end
    end
end

% image names for each column
if do_cross
    col_names = {'dat11', 'dat12', 'dat21', 'dat22'};
else
    col_names = {'dat'};
end

% write images
for m = 1:numel(test_results)
    for jj = 1:numel(res_fields{m})
        eval(['for ii = 1:size(dat_' num2str(m) '_' res_fields{m}{jj} '.dat,2), dat.dat = dat_' ...
            num2str(m) '_' res_fields{m}{jj} '.dat(:,ii); dat.fullpath = fullfile(modeldir, [''results_searchlight_' ...
            num2str(m) '_' res_fields{m}{jj} '_'' col_names{ii} ''.nii'']); write(dat); end']);
    end
end

end
