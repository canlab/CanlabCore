function [image_obj, networknames, imagenames] = load_image_set(image_names_or_keyword, varargin)
% Locate a series of images on the path and load them into an fmri_data
% object.  Useful for loading sets of canonical masks or patterns.
%
% - Checks whether images exist on path
% - Returns full image names with path names
% - Returns formatted networknames for plot labels
%
% Usage:
% ::
%
%    [imgs, names] = load_image_set(image_names_or_keyword)
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2016 Tor Wager
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
%   **image_names_or_keyword:**
%        A string matrix with images to load, or a keyword:
%        'bucknerlab': 7 network parcellation from Yeo et al. 
%        'kragelemotion': 7 emotion-predictive models from Kragel & LaBar 2015 
%        'allengenetics': Five maps from the Allen Brain Project human gene expression maps
%                         from Luke Chang (unpublished)
%        'npsplus': Wager lab multivariate patterns:
%                   NPS, PINES, Romantic Rejection, VPS
%
% :Optional inputs:
%   None yet.
%
% :Outputs:
%
%   **image_obj:**
%        fmri_data object with the maps loaded
%
%   **networknames:**
%        cell array of names based on the image names or custom titles
%
%   **imagenames:**
%        cell array of names of images loaded
%
% :Examples:
% ::
%
% imagenames = {'weights_NSF_grouppred_cvpcr.img' ...  % NPS
%     'Rating_Weights_LOSO_2.nii'  ...  % PINES
%     'dpsp_rejection_vs_others_weights_final.nii' ... % rejection
%     'bmrk4_VPS_unthresholded.nii'};
% 
% [obj, netnames, imgnames] = load_image_set(imagenames);
% 
% The above loads the same images as:
%
% [obj, netnames, imgnames] = load_image_set('npsplus');
%
% :See also:
%
% image_similarity_plot, fmri_data

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%
%   Tor: created, July 2016
%
% ..

% ..
% DEFAULTS AND INPUTS
% ..

docustom = 0;

% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isa(image_names_or_keyword, 'fmri_data')
    % We already have images loaded - just get the names
    image_obj = image_names_or_keyword;
    imagenames = image_obj.image_names;
    networknames = format_strings_for_legend(imagenames);
    return

elseif iscell(image_names_or_keyword) || (ischar(image_names_or_keyword) && size(image_names_or_keyword, 1) > 1)
    % We have custom image input

    docustom = 1;
    
else
    % we have a standard named map set
    
    switch image_names_or_keyword
        
        case 'bucknerlab'
            [image_obj, networknames, imagenames] = load_bucknerlab_maps;
            networknames=networknames';
            
        case 'npsplus'
            [image_obj, networknames, imagenames] = load_npsplus;
            
        case 'kragelemotion'
            [image_obj, networknames, imagenames] = load_kragelemotion;
            
        case 'allengenetics'
            [image_obj, networknames, imagenames] = load_allengenetics;
            
            
        otherwise
            error('Unknown mapset keyword.');
            
    end % switch
    
end % custom or not

if docustom

    [image_obj, networknames, imagenames] = load_custom(image_names_or_keyword);
    
end

disp('Loaded images:');
fprintf('%s\n', imagenames{:});




end % function



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



function [image_obj, networknames, imagenames] = load_custom(imagenames)

% Load images, whatever they are
% ------------------------------------------------------------------------
if ~iscell(imagenames), imagenames = cellstr(imagenames); end

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns

networknames = format_strings_for_legend(image_obj.image_names);

end  % function




function [image_obj, networknames, imagenames] = load_bucknerlab_maps

% Load Bucker Lab 1,000FC masks
% ------------------------------------------------------------------------

names = load('Bucknerlab_7clusters_SPMAnat_Other_combined_regionnames.mat');
img = which('rBucknerlab_7clusters_SPMAnat_Other_combined.img');

image_obj = fmri_data(img, [], 'noverbose');  % loads image with integer coding of networks

networknames = names.rnames(1:7);
k = length(networknames);

newmaskdat = zeros(size(image_obj.dat, 1), k);

for i = 1:k  % breaks up into one map per image/network
    
    wh = image_obj.dat == i;
    
    nvox(1, i) = sum(wh);
    
    newmaskdat(:, i) = double(wh);
    
    
end

image_obj.dat = newmaskdat;

imagenames = {img};
end  % function




function [image_obj, networknames, imagenames] = load_npsplus

% Load NPS, PINES, Rejection, VPS,
% ------------------------------------------------------------------------

networknames = {'NPS' 'PINES' 'Rejection' 'VPS'};

imagenames = {'weights_NSF_grouppred_cvpcr.img' ...  % NPS
    'Rating_Weights_LOSO_2.nii'  ...  % PINES
    'dpsp_rejection_vs_others_weights_final.nii' ... % rejection
    'bmrk4_VPS_unthresholded.nii'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns

end  % function







function [image_obj, networknames, imagenames] = load_kragelemotion

% Load NPS, PINES, Rejection, VPS,
% ------------------------------------------------------------------------

networknames = {'Amused' 'Angry' 'Content' 'Fearful' 'Neutral' 'Sad' 'Surprised'};

imagenames = { ...
    'mean_3comp_amused_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_angry_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_content_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_fearful_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_neutral_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_sad_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_surprised_group_emotion_PLS_beta_BSz_10000it.img'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns

end % function


function [image_obj, networknames, imagenames] = load_allengenetics

% Load Allen Brain Atlas project human genetic maps (from Luke Chang)
% ------------------------------------------------------------------------

networknames = {'5HT' 'Opioid' 'Dopamine' 'NEalpha' 'NEbeta'};

imagenames = { ...
    'Serotonin.nii' ...
    'Opioid.nii' ...
    'Dopamine.nii' ...
    'AdrenoAlpha.nii' ...
    'AdrenoBeta.nii' ...
};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns

end % function




function imagenames = check_image_names_get_full_path(imagenames)

for i = 1:length(imagenames)
    
    if exist(imagenames{i}, 'file')
        % do nothing. Sometimes which returns empty even though file
        % exists. Do not use which if returns empty. Otherwise, do.
        
        if ~isempty(which(imagenames{i}))
            imagenames{i} = which(imagenames{i});
        end
        
    else 
        fprintf('CANNOT FIND %s \n', imagenames{i})
        error('Exiting.');

    end
    
end

end % function


