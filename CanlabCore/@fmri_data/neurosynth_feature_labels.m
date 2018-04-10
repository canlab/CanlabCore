function [image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(test_dat, varargin)
% Calculates top term match with Neurosynth 2013 feature set (or topic maps)
%
% :Usage:
% ::
%
%     [image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(test_dat,[optional inputs]);
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2017 Tor Wager and Mustafa Salman and Eswar Damaraju
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
% :Inputs:
%
%   **test_dat**
%        fMRI_data object with one or more images.
%        These can be replicates (i.e., same task, diff subjects), or
%        independent images (i.e., different component maps).
%
% :Optional Inputs:
%   **'display_output'**
%        Followed by true or false; default is true
%
%   **'images_are_replicates'**
%        Followed by true or false; default is true
%        Flag for whether images are replicates (i.e., same task, diff subjects)
%        If so, does t-test across them and reports t, p-values in one
%        table
%
%   **topk:**
%        Followed by number of most negative and most positive associations
%        to save in summary table(s).  e.g., topk = 5 saves 5 most negative
%        and 5 most positive associations.
%
%   **noverbose:**      
%       Suppress printing of all loaded image names. Default is to print
%       all image names.
%
% :Outputs:
%
%   **image_by_feature_correlations:**
%        Images x Neurosynth features, Pearson's correlation values
%
%   **top_feature_tables:**
%        Cell vector of Matlab table objects with top hits for each test
%        image.  If images_are_replicates, returns one table with summary
%        statistics.
%
% :Examples:
% ::
%
% % This example loads a  set of 75 ICA component maps from Allen et al. 2012
% % and identifies the most similar Neurosynth feature maps for each ICA component:
% test_image_name = which('rest_hcp_agg__component_ica_.nii');
% test_dat = fmri_data(test_image_name);
% [image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(test_dat);
%
% ----------------------------------------------------------------------------
% % This example loads a standard set of 30 subjects' data on an emotion
% % regulation task, and identifies the most similar Neurosynth feature maps:
% test_dat = load_image_set('emotionreg');
% [image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(test_dat, 'images_are_replicates', true);


% cd('/Users/tor/Google_Drive/CanlabDataRepository/Neuroimaging_Autolabeler')
% p = genpath(pwd)
% addpath(p)


% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%
%   Tor, Mustafa Salman, and Eswar Damaraju: created 2017
%
%   2017/09/07 Stephan
%       - added (no)verbose option
%       - adapted to allow for binary mask input. results were all NaN before.
% ..



%%
% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

display_output = true;
images_are_replicates = false;               % images in test_dat are replicates of same underlying task for, e.g., diff subjects
topk = 10;                                   % take top and bottom k for each test image
verbosestr = 'verbose';                      % print loaded neurosynth image names (long list)


% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'display_output', display_output = varargin{i+1}; varargin{i+1} = [];
                
            case 'images_are_replicates', images_are_replicates = varargin{i+1}; varargin{i+1} = [];
                
            case 'topk', topk = varargin{i+1}; varargin{i+1} = [];
                
            case 'noverbose', verbosestr = 'noverbose'; 
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


%% Load Neurosynth 2013 Term Feature Set
% -------------------------------------------------------------------------
% feature_file = which('Yarkoni_2013_Neurosynth_featureset1.mat');
% if ~exist(feature_file, 'file'), error('Cannot find Yarkoni_2013_Neurosynth_featureset1.mat on Matlab path'); end
% feature_dat = load(feature_file);
% feature_dat = feature_dat.dat;
% words = feature_dat.Y_names;

[feature_dat, words] = load_image_set('neurosynth',verbosestr);
words = words';


% Reslice
test_dat = resample_space(test_dat, feature_dat);

%% Get matrix of correlations between images and features
% -------------------------------------------------------------------------
% r = images x features matrix of correlations
if numel(unique(test_dat.dat(:)))==2 && sum(unique(test_dat.dat(:))-[0; 1]) == 0 % binary masks
    
    image_by_feature_correlations = canlab_pattern_similarity(feature_dat.dat, test_dat.dat, 'correlation');
    image_by_feature_correlations = image_by_feature_correlations';
    
else % other images
    image_by_feature_correlations = canlab_pattern_similarity(test_dat.dat, feature_dat.dat, 'correlation', 'ignore_missing');
end

%% If replicates, do t-test across images

if images_are_replicates
    
    z = fisherz(image_by_feature_correlations); % Fisher's r to z
    
    [h, p, ci, stats] = ttest(z);
    
    t = stats.tstat;
    
    image_by_feature_correlations = t; % will sort on t-values later   %mean(image_by_feature_correlations);
    
end


%% Extract top features (neg and pos loadings) for each test image
% -------------------------------------------------------------------------


ntest = size(image_by_feature_correlations, 1);

top_feature_tables = cell(1, ntest);

for i = 1:ntest                             % note: ntest always == 1 if images_are_replicates
    
    testr = image_by_feature_correlations(i, :);
    [testr_sorted, indx] = sort(testr, 'ascend');
    
    len = length(testr_sorted);
    lowindx = 1:min(topk, len);
    highindx = len:-1:len-min(topk, length(testr_sorted))+1;
    
    testr_low = testr_sorted(lowindx)';
    testr_high = testr_sorted(highindx)';
    
    test_words = words(indx);               % re-sort from low to high loadings
    
    words_low = test_words(lowindx);
    words_high = test_words(highindx);
    
    top_feature_tables{i} = table(testr_low, words_low, testr_high, words_high);
    
end

% Clean up names and add P-values if appropriate
if images_are_replicates
    
    t_low = testr_low;
    t_high = testr_high;
    
    p = p(indx);    % sort.
    p_low = p(lowindx)';
    p_high = p(highindx)';
    
    top_feature_tables{1} = table(t_low, p_low, words_low, t_high, p_high, words_high);
    
end
%% Display output if requested
% -------------------------------------------------------------------------

ustr = '_____________________________________________________________________';

if display_output
    
    if isempty(test_dat.fullpath), test_dat.fullpath = 'fullpath_was_empty'; end
    
    for i = 1:ntest
        
        if length(test_dat.fullpath(i, :)) < 70
            imgPrintName = test_dat.fullpath(i, :);
        else
            imgPrintName = spm_file(test_dat.fullpath(i, :),'short70');
        end
        fprintf('Input image %i\n%s\n%s\n', i, imgPrintName, ustr);
        
        disp(top_feature_tables{i})
        
    end
    
end

end % main function
