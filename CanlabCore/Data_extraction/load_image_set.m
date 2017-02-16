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
%        'bucknerlab': 7 network parcellation from Yeo et al., cortex only
%        'bucknerlab_wholebrain': 7 networks in cortex, BG, cerebellum
%        'bucknerlab_wholebrain_plus': 7 networks in cortex, BG, cerebellum
%        + SPM Anatomy Toolbox regions + brainstem
%        'kragelemotion': 7 emotion-predictive models from Kragel & LaBar 2015
%        'allengenetics': Five maps from the Allen Brain Project human gene expression maps
%                         from Luke Chang (unpublished)
%        'npsplus': Wager lab multivariate patterns:
%                   NPS, PINES, Romantic Rejection, VPS
%        'emotionreg' : N = 30 emotion regulation sample dataset from Wager
%        et al. 2008. Each image is a contrast image for the contrast [reappraise negative vs. look negative]
%        'bgloops', 'pauli' : 5-basal ganglia parcels and 5 associated cortical
%        networks from Pauli et al. 2016
%        'bgloops17', 'pauli17' : 17-parcel striatal regions only from Pauli et al. 2016
%        'bgloops_cortex' : Cortical regions most closely associated with
%        the Pauli 5-region striatal clusters
%        'fibromyalgia':  patterns used to predict FM from Lopez Sola et al.:
%                   NPSp, FM-pain, FM-multisensory
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
    
elseif isa(image_names_or_keyword, 'image_vector')
    error('Load image set only tested for input fmri_data objects now.');
    
elseif iscell(image_names_or_keyword) || (ischar(image_names_or_keyword) && size(image_names_or_keyword, 1) > 1)
    % We have custom image input
    
    docustom = 1;
    
else
    % we have a standard named map set
    
    switch image_names_or_keyword
        
        case 'bucknerlab'
            [image_obj, networknames, imagenames] = load_bucknerlab_maps;
            networknames=networknames';
            
        case 'bucknerlab_wholebrain'
            
            [image_obj, networknames, imagenames] = load_bucknerlab_maps_wholebrain;
            networknames=networknames';
            
        case 'bucknerlab_wholebrain_plus'
            
            [image_obj, networknames, imagenames] = load_bucknerlab_wholebrain_plus_subctx;
            networknames=networknames';
          
        case 'npsplus'
            [image_obj, networknames, imagenames] = load_npsplus;
            
        case 'kragelemotion'
            [image_obj, networknames, imagenames] = load_kragelemotion;
            
        case 'allengenetics'
            [image_obj, networknames, imagenames] = load_allengenetics;
            
        case {'emotionreg' 'emotionregulation'}
            [image_obj, networknames, imagenames] = load_emotion_reg_sample;
            
        case {'bgloops17', 'pauli17'}
            [image_obj, networknames, imagenames] = load_pauli_bg17;
            
        case {'bgloops', 'pauli'}
            [image_obj, networknames, imagenames] = load_pauli_bg;

        case {'bgloops_cortex', 'pauli_cortex'}
            [image_obj, networknames, imagenames] = load_pauli_bg_cortex;
            
        case {'fibromyalgia','fibro','fm'}    
            [image_obj, networknames, imagenames] = load_fibromyalgia;
       
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


function imagenames = check_image_names_get_full_path(imagenames)

if ~iscell(imagenames), imagenames = cellstr(imagenames); end

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



function image_obj = integer_coded_image_to_separate_images(image_obj)

u = unique(image_obj.dat);
u(u == 0) = [];
k = length(u);

newmaskdat = zeros(size(image_obj.dat, 1), k);

for i = 1:k  % breaks up into one map per image/network
    
    wh = image_obj.dat == i;
    
    newmaskdat(:, i) = double(wh);
    
    
end

image_obj.dat = newmaskdat;

end % function

function [image_obj, networknames, imagenames] = load_custom(imagenames)

% Load images, whatever they are
% ------------------------------------------------------------------------
if ~iscell(imagenames), imagenames = cellstr(imagenames); end

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns

networknames = format_strings_for_legend(image_obj.image_names);

end  % function



% -------------------------------------------------------------------------
%
% Specific, named image sets
%
% -------------------------------------------------------------------------

function [image_obj, networknames, imagenames] = load_bucknerlab_maps

% Load Bucker Lab 1,000FC masks
% ------------------------------------------------------------------------

names = load('Bucknerlab_7clusters_SPMAnat_Other_combined_regionnames.mat');
img = which('rBucknerlab_7clusters_SPMAnat_Other_combined.img');

image_obj = fmri_data(img, [], 'noverbose');  % loads image with integer coding of networks

networknames = names.rnames(1:7); % cortex only
k = length(networknames);

newmaskdat = zeros(size(image_obj.dat, 1), k);

for i = 1:k  % breaks up into one map per image/network
    
    wh = image_obj.dat == i;
    
    %nvox(1, i) = sum(wh);
    
    newmaskdat(:, i) = double(wh);
    
    
end

image_obj.dat = newmaskdat;

imagenames = {img};
end  % function



function [image_obj, networknames, imagenames] = load_bucknerlab_maps_wholebrain

% Load Bucker Lab 1,000FC masks
% ------------------------------------------------------------------------

names = load('Bucknerlab_7clusters_SPMAnat_Other_combined_regionnames.mat');
img = which('rBucknerlab_7clusters_SPMAnat_Other_combined.img');

image_obj = fmri_data(img, [], 'noverbose');  % loads image with integer coding of networks

networknames = names.rnames(1:7); 

k = length(networknames);       % cortex, striatum, cerebellum, same names

newmaskdat = zeros(size(image_obj.dat, 1), k);

% Cortex, BG, CBLM
for i = 1:7  % breaks up into one map per image/network
    
    wh = image_obj.dat == i;
    newmaskdat(:, i) = double(wh);
  
    wh = image_obj.dat == i + 7;
    newmaskdat(:, i) = newmaskdat(:, i) + double(wh);
    
    wh = image_obj.dat == i + 14;
    newmaskdat(:, i) = newmaskdat(:, i) + double(wh);
end

% ADD OTHER REGIONS
% *****


image_obj.dat = newmaskdat;

imagenames = {img};
end  % function




function [image_obj, networknames, imagenames] = load_bucknerlab_wholebrain_plus_subctx

% Load Bucker Lab 1,000FC masks
% ------------------------------------------------------------------------

names = load('Bucknerlab_7clusters_SPMAnat_Other_combined_regionnames.mat');
img = which('rBucknerlab_7clusters_SPMAnat_Other_combined.img');

image_obj = fmri_data(img, [], 'noverbose');  % loads image with integer coding of networks

networknames = names.rnames(1:7); 

m = 5;                              % number of other regions
                                    % Amy, Thal, Hy, Brainstem, Hippocampus
k = length(networknames) + m;       % cortex, striatum, cerebellum

newmaskdat = zeros(size(image_obj.dat, 1), k);

% Cortex, BG, CBLM
for i = 1:7  % breaks up into one map per image/network
    
    wh = image_obj.dat == i;
    newmaskdat(:, i) = double(wh);
  
    wh = image_obj.dat == i + 7;
    newmaskdat(:, i) = newmaskdat(:, i) + double(wh);
    
    wh = image_obj.dat == i + 14;
    newmaskdat(:, i) = newmaskdat(:, i) + double(wh);
end

% ADD OTHER REGIONS
% *****
% find hipp:
ishipp = ~cellfun(@isempty, strfind(names.rnames, 'Hipp'));
newmaskdat(:, end + 1) = double(any(image_obj.dat(:, ishipp), 2));


image_obj.dat = newmaskdat;

imagenames = {img};  % ***add to names
end  % function





function [image_obj, networknames, imagenames] = load_npsplus

% Load NPS, PINES, Rejection, VPS,
% ------------------------------------------------------------------------

networknames = {'NPS' 'NPSpos' 'NPSneg' 'SIIPS' 'PINES' 'Rejection' 'VPS' 'GSR' 'Heart' 'FM-Multisens' 'FM-pain'};

imagenames = {'weights_NSF_grouppred_cvpcr.img' ...     % NPS   - somatic pain
    'NPSp_Lopez-Sola_2017_PAIN.img' ...                 % 2017 Lopez-Sola positive NPS regions only
    'NPSn_Lopez-Sola_2017_PAIN.img' ...                 % 2017 Lopez-Sola negative NPS regions only, excluding visual
    'nonnoc_v11_4_137subjmap_weighted_mean.nii' ...     % SIIPS - stim-indep pain
    'Rating_Weights_LOSO_2.nii'  ...                    % PINES - neg emo
    'dpsp_rejection_vs_others_weights_final.nii' ...    % romantic rejection
    'bmrk4_VPS_unthresholded.nii' ...                   % Vicarious pain
    'ANS_Eisenbarth_JN_2016_GSR_pattern.img' ...        % autonomic - GSR
    'ANS_Eisenbarth_JN_2016_HR_pattern.img' ...         % autonomic - heart rate (HR)
    'FM_Multisensory_wholebrain.nii' ...                % 2017 Lopez-Sola fibromyalgia 
    'FM_pain_wholebrain.nii'};                          % 2017 Lopez-Sola fibromyalgia 

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose', 'sample2mask');  % loads images with spatial basis patterns

end  % function







function [image_obj, networknames, imagenames] = load_kragelemotion

% Load Kragel 2015 emotion maps
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




function [image_obj, networknames, imagenames] = load_emotion_reg_sample

% Load Wager et al. 2008 Emotion Regulation sample dataset
% ------------------------------------------------------------------------

myfile = which('con_00810001.img');
mydir = fileparts(myfile);

if isempty(mydir)
    disp('Uh-oh! I can''t find the data.')
else
    disp('Data found.')
end

imagenames = filenames(fullfile(mydir, 'con_008100*img'));

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');  % loads images

networknames = format_strings_for_legend(image_obj.image_names);

end % function






function [image_obj, networknames, imagenames] = load_pauli_bg

% Load Pauli et al. 2016 basal ganglia 5-cluster solution
% ------------------------------------------------------------------------

networknames = {'Post. Caudate (Cp)' 'Ant. Putamen (Pa)' 'Ant. Caudate (Ca)' 'Ventral striatum (VS)' 'Post. Putamen (PP)'}; 

imagenames = {'Pauli_bg_cluster_mask_5.nii'};
imagenames = check_image_names_get_full_path(imagenames);

% load image with integer coding of networks
image_obj = fmri_data(imagenames, [], 'noverbose');

% Break up into one image per region
% -------------------------------------------------
image_obj = integer_coded_image_to_separate_images(image_obj);


end % function



function [image_obj, networknames, imagenames] = load_pauli_bg17

% Load Pauli et al. 2016 basal ganglia 17-cluster solution (no labels
% given in the paper)
% ------------------------------------------------------------------------
networknames = {'cluster 1','cluster 2','cluster 3','cluster 4','cluster 5','cluster 6',...
    'cluster 7','cluster 8','cluster 9','cluster 10','cluster 11','cluster 12','cluster 13',...
    'cluster 14','cluster 15','cluster 16','cluster 17',};

imagenames = {'Pauli_bg_cluster_mask_17.nii'};
imagenames = check_image_names_get_full_path(imagenames);

% load image with integer coding of networks
image_obj = fmri_data(imagenames, [], 'noverbose');

% Break up into one image per region
% -------------------------------------------------
image_obj = integer_coded_image_to_separate_images(image_obj);


end % function



function [image_obj, networknames, imagenames] = load_pauli_bg_cortex

% Load Pauli et al. 2016 basal ganglia 5-cluster solution
% ------------------------------------------------------------------------
imagenames = {'Pauli_bg_nb_param_rank_fst_Cp.nii' ...
                'Pauli_bg_nb_param_rank_fst_Pa.nii' ...
                'Pauli_bg_nb_param_rank_fst_Ca.nii' ...
                'Pauli_bg_nb_param_rank_fst_VS.nii' ...
                'Pauli_bg_nb_param_rank_fst_Pp.nii' ...
};

networknames = {'Post. Caudate (Cp)' 'Ant. Putamen (Pa)' 'Ant. Caudate (Ca)' 'Ventral striatum (VS)' 'Post. Putamen (PP)'}; 

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function



function [image_obj, networknames, imagenames] = load_fibromyalgia
  

% Load Lopez Sola et al. 2017 neural classifier maps
% ------------------------------------------------------------------------
imagenames = {'FM_pain_wholebrain.nii' ...
                'FM_Multisensory_wholebrain.nii' ...
                'rNPS_fdr_pospeaks_smoothed.img' };

networknames = {'FM-pain' 'FM-multisensory' 'NPSp'}; 

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function