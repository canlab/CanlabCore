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
%        A string matrix with images to load, or a keyword.
%        keywords load pre-defined image sets, as indicated below.
%        NOTE: you will need to have these images on your Matlab path!
%        Some are in the CANlab Neuroimaging_Pattern_Masks repository, 
%        some in Masks_Private repository 
% 
%        'bucknerlab': 7 network parcellation from Yeo et al., cortex only
% 
%        'bucknerlab_wholebrain': 7 networks in cortex, BG, cerebellum
% 
%        'bucknerlab_wholebrain_plus': 7 networks in cortex, BG, cerebellum + SPM Anatomy Toolbox regions + brainstem
% 
%        'kragelemotion': 7 emotion-predictive models from Kragel & LaBar 2015
% 
%        'allengenetics': Five maps from the Allen Brain Project human gene expression maps
%                         from Luke Chang (unpublished)
% 
%        'npsplus': Wager lab published multivariate patterns:
%                   NPS (incl NPSpos & NPSpos), SIIPS, PINES, Romantic Rejection, VPS, more
% 
%        'painsig': NPS (incl NPSpos & NPSpos) and SIIPS only
% 
%        'emotionreg' : N = 30 emotion regulation sample dataset from Wager et al. 2008. 
%                       Each image is a contrast image for the contrast [reappraise negative vs. look negative]
% 
%        'bgloops', 'pauli' : 5-basal ganglia parcels and 5 associated cortical
%                             networks from Pauli et al. 2016
% 
%        'bgloops17', 'pauli17' : 17-parcel striatal regions only from Pauli et al. 2016
% 
%        'bgloops_cortex' : Cortical regions most closely associated with
%                           the Pauli 5-region striatal clusters
% 
%        'fibromyalgia':  patterns used to predict FM from Lopez Sola et al.:
%                         NPSp, FM-pain, FM-multisensory
% 
%        'neurosynth', 'neurosynth_featureset1': 525 "Reverse inference" z-score maps from Tal Yarkoni's
%                                                Neurosynth, unthresholded, 2013 
% 
%        'pain_cog_emo', 'kragel18': Partial least squares maps for generalizable representations of pain, cog control, emotion 
% 
%        'pain_pdm', 'pdm': High-dimensional mediators of pain. 10 individual PDM maps and a joint
%                           PDM, which is a weighted combination of the 10. From Geuter et al. (in prep)
%           
% 
%
% :Optional inputs:
%
%   **noverbose:**
%       Suppress printing of all loaded image names. Default is to print
%       all image names.
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
%  % Example 1: Load NPS (private) and several other signatures
% % ------------------------------------------------------------------------- 
% imagenames = {'weights_NSF_grouppred_cvpcr.img' ...  % NPS
%     'Rating_Weights_LOSO_2.nii'  ...  % PINES
%     'dpsp_rejection_vs_others_weights_final.nii' ... % rejection
%     'bmrk4_VPS_unthresholded.nii'};
%
% [obj, netnames, imgnames] = load_image_set(imagenames);
%
% The above loads a subset of the same images as:
%
% [obj, netnames, imgnames] = load_image_set('npsplus');
%
% % Example 2: Apply the PLS signatures from Kragel et al. 2018 to the emotion regulation dataset
% % -------------------------------------------------------------------------
% % Load PLS signatures from Kragel et al. 2018
% [obj, names] = load_image_set('pain_cog_emo');
% bpls_wholebrain = get_wh_image(obj, [8 16 24]);
% names_wholebrain = names([8 16 24]);
% bpls_subregions = get_wh_image(obj, [1:6 9:14 17:22]);
% names_subregions = names([1:6 9:14 17:22]);
% 
% % Load test data: Emotion regulation from Wager et al. 2008
% test_data_obj = load_image_set('emotionreg');
% 
% %  Make plots
% % Yellow: positive associations. Blue: Negative associations.  Plot shows mean +- std. error for each pattern of interest
%   
% create_figure('Kragel Pain-Cog-Emo maps', 1, 2);
% stats = image_similarity_plot(test_data_obj, 'average', 'mapset', bpls_wholebrain, 'networknames', names_wholebrain, 'nofigure');
% subplot(1, 2, 2)
% stats = image_similarity_plot(test_data_obj, 'average', 'mapset', bpls_subregions, 'networknames', names_subregions, 'nofigure');
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
%   2017/09/07 Stephan
%       - added (no)verbose option
% ..

% ..
% DEFAULTS AND INPUTS
% ..

docustom = 0;
verbose = 1;

% optional inputs with default values
% -----------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'noverbose', verbose = 0;
            
            case 'verbose', verbose = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isa(image_names_or_keyword, 'fmri_data') || isa(image_names_or_keyword, 'atlas')
    % We already have images loaded - just get the names
    image_obj = image_names_or_keyword;
    imagenames = image_obj.image_names;
    networknames = format_strings_for_legend(imagenames);
    if iscolumn(networknames), networknames = networknames'; end

    return
    
elseif isa(image_names_or_keyword, 'image_vector')
    error('Load image set only tested for input fmri_data objects now.');
    
elseif iscell(image_names_or_keyword) || (ischar(image_names_or_keyword) && size(image_names_or_keyword, 1) > 1)
    % We have custom image input
    
    docustom = 1;
    
else
    % we have a standard named map set
    
    switch lower(image_names_or_keyword)
        
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
            
        case 'painsig'
            [image_obj, networknames, imagenames] = load_painsig;
            
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
          
           case {'pauli_subcortical'}
            [image_obj, networknames, imagenames] = load_pauli_subcortical;
                   
        case {'fibromyalgia','fibro','fm'}    
            [image_obj, networknames, imagenames] = load_fibromyalgia;
       
        case {'neurosynth', 'neurosynth_featureset1'}
            [image_obj, networknames, imagenames] = load_neurosynth_featureset1;
            
        case {'pain_cog_emo', 'kragel18'}
             [image_obj, networknames, imagenames] = load_kragel18;
             
        case {'pdm','pain_pdm'}
             [image_obj, networknames, imagenames] = load_pain_pdm;

        otherwise
            error('Unknown mapset keyword.');
            
    end % switch
    
end % custom or not

if docustom
    
    [image_obj, networknames, imagenames] = load_custom(image_names_or_keyword);
    
end

if verbose
    disp('Loaded images:');
    fprintf('%s\n', imagenames{:});
end

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

% remove trailing commas indexing volume number from any 4-d image names
pat = ',\w*';
imagenames = regexprep(imagenames, pat, '');

for i = 1:length(imagenames)
    
    if exist(imagenames{i}, 'file')
        
        % do nothing. Sometimes which returns empty even though file
        % exists. Do not use which if returns empty. Otherwise, do.
        
        if ~isempty(which(imagenames{i}))
            imagenames{i} = which(imagenames{i});
        end
        
    else
        
        % check for .gz version; we entered .img/.nii but have .gz
        myimg = which(deblank([imagenames{i} '.gz']));
        
        if ~isempty(myimg)
            imagenames{i} = myimg;
        end
            
    end
    
    if ~exist(imagenames{i}, 'file')
        
        fprintf('CANNOT FIND IMAGES %s \n', imagenames{i})
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

if iscolumn(networknames), networknames = networknames'; end

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

networknames = {'NPS' 'NPSpos' 'NPSneg' 'SIIPS' 'PINES' 'Rejection' 'VPS' 'VPS_nooccip' 'GSR' 'Heart' 'FM-Multisens' 'FM-pain' 'Empathic_Dist' 'Empathic_Care'};

imagenames = {'weights_NSF_grouppred_cvpcr.img' ...     % Wager et al. 2013 NPS   - somatic pain
    'NPSp_Lopez-Sola_2017_PAIN.img' ...                 % 2017 Lopez-Sola positive NPS regions only
    'NPSn_Lopez-Sola_2017_PAIN.img' ...                 % 2017 Lopez-Sola negative NPS regions only, excluding visual
    'nonnoc_v11_4_137subjmap_weighted_mean.nii' ...     % Woo 2017 SIIPS - stim-indep pain
    'Rating_Weights_LOSO_2.nii'  ...                    % Chang 2015 PINES - neg emo
    'dpsp_rejection_vs_others_weights_final.nii' ...    % Woo 2014 romantic rejection
    'bmrk4_VPS_unthresholded.nii' ...                   % Krishnan 2016 Vicarious pain VPS
    'Krishnan_2016_VPS_bmrk4_Without_Occipital_Lobe.nii' ... % Krishnan 2016 no occipital
    'ANS_Eisenbarth_JN_2016_GSR_pattern.img' ...        % Eisenbarth 2016 autonomic - GSR
    'ANS_Eisenbarth_JN_2016_HR_pattern.img' ...         % Eisenbarth 2016 autonomic - heart rate (HR)
    'FM_Multisensory_wholebrain.nii' ...                % 2017 Lopez-Sola fibromyalgia 
    'FM_pain_wholebrain.nii' ...                        % 2017 Lopez-Sola fibromyalgia 
    'Ashar_2017_empathic_care_marker.nii' ...           % 2017 Ashar et al. Empathic care and distress
    'Ashar_2017_empathic_distress_marker.nii'};         

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose', 'sample2mask');  % loads images with spatial basis patterns

end  % function


function [image_obj, networknames, imagenames] = load_painsig

% Load pain signatures 
% ------------------------------------------------------------------------

networknames = {'NPS' 'NPSpos' 'NPSneg' 'SIIPS'};

imagenames = {'weights_NSF_grouppred_cvpcr.img' ...     % Wager et al. 2013 NPS   - somatic pain
    'NPSp_Lopez-Sola_2017_PAIN.img' ...                 % 2017 Lopez-Sola positive NPS regions only
    'NPSn_Lopez-Sola_2017_PAIN.img' ...                 % 2017 Lopez-Sola negative NPS regions only, excluding visual
    'nonnoc_v11_4_137subjmap_weighted_mean.nii'};    % Woo 2017 SIIPS - stim-indep pain         

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose', 'sample2mask');  % loads images with spatial basis patterns

end  % function


function [image_obj, networknames, imagenames] = load_kragel18
% Load NPS, PINES, Rejection, VPS,
% ------------------------------------------------------------------------

domains = {'Pain' 'Cognitive Control' 'Negative Emotion'};
rois = {'pMCC' 'aMCC' 'pACC' 'sgACC' 'vmPFC' 'dMFC' 'MFC' 'Wholebrain'};

networknames = {'Pain pMCC' 'Pain aMCC' 'Pain pACC' 'Pain sgACC' 'Pain vmPFC' 'Pain dMFC' 'Pain MFC' 'Pain Wholebrain' ,...
    'Cog pMCC' 'Cog aMCC' 'Cog pACC' 'Cog sgACC' 'Cog vmPFC' 'Cog dMFC' 'Cog MFC' 'Cog Wholebrain',...
    'Emo pMCC' 'Emo aMCC' 'Emo pACC' 'Emo sgACC' 'Emo vmPFC' 'Emo dMFC' 'Emo MFC' 'Emo Wholebrain'};
imagenames=cell(21,1);
cs=0;
for d=1:length(domains)
    for r=1:length(rois)
        cs=cs+1;
        imagenames{cs} = ['bPLS_' rois{r} '_' strrep(domains{d}, ' ', '_') '.nii'];
    end
end

% for i=1:length(imagenames) %preserve order..
% imagenames{i} = which(imagenames{i});
% end

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

% myfile = which('con_00810001.img');
% mydir = fileparts(myfile);

myfile = which('Wager_2008_emo_reg_vs_look_neg_contrast_images.nii.gz');

if isempty(myfile)
    myfile = which('Wager_2008_emo_reg_vs_look_neg_contrast_images.nii');
end

if isempty(myfile)
    disp('Uh-oh! I can''t find the data.')
else
    % found data ok
end

% imagenames = filenames(fullfile(mydir, 'con_008100*img'));
% 
% imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(myfile, [], 'noverbose');  % loads images

imagenames = image_obj.fullpath;
imagenames = check_image_names_get_full_path(imagenames);

networknames = format_strings_for_legend(image_obj.image_names);

if length(networknames) == 1 % 4-d, expand
   
    networknames = repmat(networknames, size(image_obj.dat, 2), 1);
    
end
    
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


function [image_obj, networknames, imagenames] = load_pauli_subcortical

% Load Pauli et al. 2016 basal ganglia 17-cluster solution (no labels
% given in the paper)
% ------------------------------------------------------------------------
networknames = {'Pu' 'Ca' 'NAC' 'EXA' 'GPe' 'GPi' 'VeP' 'VTA' 'SNc' 'SNr' 'PBP' 'RN' 'HTH' 'STH' 'HN' 'MN'};

imagenames = {'pauli_subcortical.nii'};
imagenames = check_image_names_get_full_path(imagenames);

% load image with integer coding of networks
image_obj = fmri_data(imagenames, [], 'noverbose');

% Break up into one image per region
% -------------------------------------------------


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




function [image_obj, networknames, imagenames] = load_pain_pdm
  

% Load Geuter et al. 2018 high-dimensional pain mediator map (PDM)
% ------------------------------------------------------------------------
imagenames = {'JointPDM_unthresholded.nii' 'All_PDM10_unthresholded.nii'};

networknames = {'GeuterJointPDM','PDM1','PDM2','PDM3','PDM4','PDM5','PDM6','PDM7','PDM8','PDM9','PDM10'}; 

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function



function [image_obj, networknames, imagenames] = load_neurosynth_featureset1
  
% Load Yarkoni_2013_Neurosynth_featureset1
% ------------------------------------------------------------------------

datfilename = 'Yarkoni_2013_Neurosynth_featureset1.mat';
fullfilename = which('Yarkoni_2013_Neurosynth_featureset1.mat');

if isempty(fullfilename)
    disp('Cannot find required data image file.')
    fprintf('Find and add the file %s to your Matlab path.', datfilename);
    error('Exiting');
end

% Generic load fmri_data object with error checking
image_obj = [];
tmpstruct = load(fullfilename);

N = fieldnames(tmpstruct);

for i = 1:length(N)
   if isa(tmpstruct.(N{i}), 'fmri_data')
       image_obj = tmpstruct.(N{i});
   end
end
if isempty(image_obj)
        fprintf('File %s does not contain any fmri_data objects.', datfilename);
    error('Exiting');
end

imagenames = cellstr(image_obj.image_names);

networknames = imagenames';
networknames = cellfun(@(x) strrep(x, '_pFgA_z.nii', ''), networknames, 'UniformOutput', false);

end % function
