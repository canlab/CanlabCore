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
%        some in Masks_Private repository, other (unlisted) datasets can be
%        loaded if you have a load_<dataset>.m file in your path. This is
%        for simplified extensions to this method by other libraries.
%
% Sample test datasets - one image per subject
% ------------------------------------------------------------------------
%       'list' : List of signatures (enter any keyword from list as input)
%                Note: returns table as 1st output instead of image_obj
%
% Sample test datasets - one image per subject
% ------------------------------------------------------------------------
%        'emotionreg' : N = 30 emotion regulation sample dataset from Wager et al. 2008.
%                       Each image is a contrast image for the contrast [reappraise negative vs. look negative]
%
%        'bmrk3', 'pain' : 33 participants, with brain responses to six levels of heat (non-painful and painful).
%                  NOTE: requires access to bmrk3_6levels_pain_dataset.mat,
%                  on figshare (see canlab.github.io/walkthroughs)
%
%        'kragel18_alldata' : 270 subject maps from Kragel 2018;
%                             These are saved in kragel_2018_nat_neurosci_270_subjects_test_images.mat
%                             if not found, will attempt to download from Neurovault using
%                             retrieve_neurovault_collection().
%
% Sample test datasets - one image per trial (single trial datasets)
% ------------------------------------------------------------------------
%     A set of single-trial datasets for pain studies have been compiled by Bogdan Petre and stored here:
%     https://github.com/canlab/canlab_single_trials
%
%     Each dataset has a name (e.g., 'nsf', 'exp', 'bmrk3pain'), and you can enter
%     any of these names as keywords, or 'all_single_trials' to load all of
%     them. The canlab_single_trials repo must be on your matlab path.
%     Each study file loads as an fmri_data object, with a metadata_table field
%     that stores trial info in a Matlab table-class object.
%
%     Single-trial datasets include:
%     % 'nsf' 'bmrk3pain' 'bmrk3warm' 'bmrk4' 'exp' 'ie' 'ie2' 'ilcp'
%     'romantic' 'scebl' 'stephan'
%
%     Each contains a single-trial object from a differnt study.
%
% Parcellations and large-scale networks/patterns
% ------------------------------------------------------------------------
%        'bucknerlab': 7 network parcellation from Yeo et al., cortex only
%
%        'bucknerlab_wholebrain': 7 networks in cortex, BG, cerebellum
%
%        'bucknerlab_wholebrain_plus': 7 networks in cortex, BG, cerebellum + SPM Anatomy Toolbox regions + brainstem
%
%        'allengenetics': Five maps from the Allen Brain Project human gene expression maps
%                         from Luke Chang (unpublished)
%
%        'bgloops', 'pauli' : 5-basal ganglia parcels and 5 associated cortical
%                             networks from Pauli et al. 2016
%
%        'bgloops17', 'pauli17' : 17-parcel striatal regions only from Pauli et al. 2016
%
%        'bgloops_cortex' : Cortical regions most closely associated with
%                           the Pauli 5-region striatal clusters
%
%        'pet_nr_map', 'hansen22' 'pet' 'receptorbinding' :     2022_Hansen_PET_tracer_maps, 36 maps with
%        different combinations of tracers and neurotransmitter receptors
%
%        'emometa' 'emotionmeta' '2015emotionmeta' : 2015 Wager/Kang et al.
%        Meta-analysis maps for 5 basic categories.
%        'Anger' 'Disgust' 'Fear' 'Happy' 'Sad'
%
%        'marg' or 'transmodal' or 'principalgradient': 2016 PNAS Margulies
%                    MNI152NLin2009cAsym_margulies_grad1.nii.gz
%        'margfsl' : MNI152NLin6Asym_margulies_grad1.nii.gz
%
%        'transcriptomic_gradients':
%               Principal transcriptomic gradients from
%               Hawrylycz et al. An anatomically comprehensive atlas of the adult human brain transcriptome. (2012) Nature
%               Vogel et al. "Deciphering the functional specialization of whole-brain spatiomolecular gradients in the adult brain" (2024) PNAS
%               [transcriptomic_grads transcriptomic_names] = load_image_set('transcriptomic_gradients');
%
% 'Signature' patterns and predictive models
% ------------------------------------------------------------------------
%        'list'     : Return list of signatures in a table
%        'all'      : Load all signatures in table registry (see 'list')
%        'nps'      : Wager et al. 2013 Neurologic Pain Signature
%        'vps'      : Krishnan et et al. 2016 Vicarious Pain Signature
%        'rejection': Woo et et al. 2014 Romantic Rejection
%        'siips'    : Woo et et al. 2017 Stimulus intensity-independent pain Signature
%        'pines'    : Chang et et al. 2015 Picture-induced negative emotion Signature
%        'gsr'      : Eisenbarth et al. 2016 Stress-induced skin conductance
%        'hr'       : Eisenbarth et al. 2016 Stress-induced heart rate
%        'multisensory' : Lopez-sola et al. 2016 Fibromyalgia multisensory pattern
%        'fmpain'   : Lopez-sola et al. 2016 Fibromyalgia pain-period pattern
%        'plspain'  : Kragel el al. 2018 PLS pain-related
%        'cpdm'     : Geuter et al. 2020 multivariate mediation pain-related
%
%        'npsplus'  : Wager lab published multivariate patterns:
%                   NPS (incl NPSpos & NPSpos), SIIPS, PINES, Romantic Rejection, VPS, more
%
%        'painsig'  : NPS (incl NPSpos & NPSpos) and SIIPS only
%
%        'fibromyalgia':  patterns used to predict FM from Lopez Sola et al.:
%                         NPSp, FM-pain, FM-multisensory
%
%        'guilt'    : a multivariate fMRI pattern related to guilt behavior
%                   Yu, Koban et al. 2019, Cerebral Cortex
%                   Yu_guilt_SVM_sxpo_sxpx_EmotionForwardmask.nii.gz
%
%        'neurosynth', 'neurosynth_featureset1': 525 "Reverse inference" z-score maps from Tal Yarkoni's
%                                                Neurosynth, unthresholded, 2013
%
%        'neurosynth_topics_forwardinference' 'neurosynth_topics_reverseinference'
%                   54 topic maps from Yarkoni & Poldrack 2014 topic modeling analysis
%                   Selected from 100 topics for psychological relevance
%                   and given ChatGPT-based summary topic labels by Ke et al. 2024, Nat Neurosci
%
%        'pain_cog_emo', 'kragel18': Partial least squares maps for generalizable
%                   representations of pain, cog control, emotion. From
%                   Kragel et al. 2018, Nature Neuroscience
%
%        'pain_pdm', 'pdm': High-dimensional mediators of pain. 10
%        individual PDM maps and a combined PDM, which is a weighted
%        combination of the 10. From Geuter et al. (2020) Cerebral Cortex
%        (see cpdm above)
%
%        'kragelemotion': 7 emotion-predictive models from Kragel & LaBar 2015
%
%        'kragelschemas': 20 visual emotion-schemas from Kragel et al. 2019
%
%        {'reddanCSplus' 'threat'}: Reddan et al. 2018 Neuron CS+ vs. CS- classifier map
%        'zhouvps':  Zhou et al. 2020 eLife generalized vicarious pain signature  
%
%        'multiaversive', 'mpa2': Ceko et al. multiple predictive patterns
%        for aversive experience: General, Mechanical pain,
%        Aversive Sounds, Thermal pain, Visual aversive images
%
%        'stroop': Silvestrini et al. 2020 Stroop-demand SVM. stroop_pattern_wani_121416.nii
%
%        'ncs':     drug and food craving signature(s)
%                   Koban et al, Nature Neuroscience 2022
%                   craving_wmapN99_boot10K_02-May-2022.img
%                   wmap_onlyDRUGS_l2nGM_N99_20220428.img
%                   wmap_onlyFOOD_l2nGM_N99_20220428.img
%
%
% :Optional inputs:
%
%   **noverbose:**
%       Suppress printing of all loaded image names. Default is to print
%       all image names.
%
%  **md5check:**
%       Perform md5 hash check if supported for dataset. If verbosity is
%       enabled md5 check results will be returned to stdout.
%
%  **forcedl:**
%       Force download without prompting for permission if dataset is
%       missing.
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
%
%   2020/10/28 Tor
%       - Corrected error in Kragel270 all data thermal vs. visceral labels
%       switched.<Note: switched back as phil updated Neurovault order,
%       fixed another bug in study order saving in metadata table.
%
%   2020/11/10 Tor
%       - Corrected kragel270 image load to retrieve from Google drive.
%       Added more descriptive study name codes.
%
%   2023/01/23 Lukas
%       -   Added NCS craving signatures (Koban et al, Nature Neuroscience 2020)
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
                
            case 'md5check', continue; % only supported by custom functions in other extension repositories
                
            case 'forcedl', continue; % only supported by custom functions in other extension repositories.
                
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
            
        case {'nps', 'pines' 'vps', 'rejection', 'siips' 'hr' 'gsr' 'cpdm' 'fmpain' 'multisensory' 'plspain'}
            
            [image_obj, networknames, imagenames] = load_signature(image_names_or_keyword);
            
        case 'npsplus'
            [image_obj, networknames, imagenames] = load_npsplus;
            
        case 'painsig'
            [image_obj, networknames, imagenames] = load_painsig;
            
        case 'kragelemotion'
            [image_obj, networknames, imagenames] = load_kragelemotion;
            
        case 'kragelschemas'
            [image_obj, networknames, imagenames] = load_kragelschemas;
            
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
            
        case {'neurosynth_topics_forwardinference' 'neurosynth_topics_fi'}

            [image_obj, networknames, imagenames] = load_neurosynth_topics_fi;

        case {'neurosynth_topics_reverseinference' 'neurosynth_topics_ri'}

            [image_obj, networknames, imagenames] = load_neurosynth_topics_ri;

        case {'pain_cog_emo', 'kragel18'}
            [image_obj, networknames, imagenames] = load_kragel18;
            
        case {'pdm','pain_pdm'}
            [image_obj, networknames, imagenames] = load_pain_pdm;
            
        case {'bmrk3', 'pain'}
            [image_obj, networknames, imagenames] = load_bmrk3;
            
        case {'kragel270' 'kragel2018_alldata' 'kragel18_alldata' 'kragel18_testdata'}
            
            [image_obj, networknames, imagenames] = load_kragel18_alldata;
            
        case {'guilt', 'guilt_behavior'}
            
            [image_obj, networknames, imagenames] = load_guilt;
            
        case {'multiaversive', 'mpa2'}
            
            [image_obj, networknames, imagenames] = load_mpa2;
            
        case {'emometa' 'emotionmeta' '2015emotionmeta'}
            
            [image_obj, networknames, imagenames] = load_emotionmeta;

        case {'pet_nr_map', 'hansen22' 'pet' 'receptorbinding'}
            
            [image_obj, networknames, imagenames] = load_hansen22;

        case 'stroop'
            
            [image_obj, networknames, imagenames] = load_stroop;
            
        case 'ncs'
            
            [image_obj, networknames, imagenames] = load_ncs;
            
        case {'marg', 'transmodal', 'principalgradient'}
            
            [image_obj, networknames, imagenames] = load_marg;

        case {'transcriptomic_gradients'}
            [image_obj, networknames, imagenames] = load_transcriptomic_gradients;

        case {'margfsl'}
            
            [image_obj, networknames, imagenames] = load_margfsl;

        case 'list'
            
            [networknames, imagenames] = deal({});

            table_list = list_signatures;
            disp(table_list);
            
            image_obj = table_list; % return as 1st output
            return
            
        otherwise
            % Try to load if we have a function name
            % This will be true for the single trial data (Bogdan Petre
            % repository) if you have it on your Matlab path
            
            if which(['load_', lower(image_names_or_keyword)])
                
                [image_obj, networknames, imagenames] = feval(['load_', lower(image_names_or_keyword)],...
                    'verbose', verbose, varargin{:});
                
            else
                % Try to load if it matches any in list
                
                table_list = list_signatures;
                
                wh = find(strcmp(image_names_or_keyword, table_list.keyword)); %#ok<*EFIND>
                
                if ~isempty(wh)
                    [image_obj, networknames, imagenames] = load_signature(image_names_or_keyword);
                    
                else
                    error('Unknown mapset keyword. If entering image names, use a cell array. Try ''list'' for a list of signatures');
                end
            end
            
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



function [image_obj, networknames, imagenames] = load_signature(keyword)

% Load one of a set of signatures in table (see function at bottom for
% table registry)

% ------------------------------------------------------------------------

table_list = list_signatures;

networknames = table_list.keyword';
imagenames = table_list.imagenames';

switch lower(keyword)
    case 'nps'
        wh = 1;
    case 'siips'
        wh = 4;
    case 'pines'
        wh = 5;
    case 'rejection'
        wh = 6;
    case 'vps'
        wh = 7;
    case 'gsr'
        wh = 9;
    case 'hr'
        wh = 10;
    case 'multisensory'
        wh = 11;
    case 'fmpain'
        wh = 12;
    case 'plspain'
        wh = 15;
    case 'cpdm'
        wh = 27;
        
    otherwise
        wh = find(strcmp(keyword, table_list.keyword));
        
end

if isempty(wh)
    
    error('Only some signature keywords are implemented in load_image_set, and the one you entered is invalid. Try ''npsplus''');
    
end

imagenames = check_image_names_get_full_path(imagenames(wh));

image_obj = fmri_data(imagenames, [], 'noverbose', 'sample2mask');  % loads images with spatial basis patterns

networknames = networknames(wh);

end  % function


function [image_obj, networknames, imagenames] = load_all(varargin)

% Load all signatures
% ------------------------------------------------------------------------

table_list = list_signatures;

networknames = table_list.keyword';
imagenames = table_list.imagenames';

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose', 'sample2mask');  % loads images with spatial basis patterns

end  % function


function [image_obj, networknames, imagenames] = load_npsplus

% Load NPS, PINES, Rejection, VPS, more
% ------------------------------------------------------------------------

networknames = {'NPS' 'NPSpos' 'NPSneg' 'SIIPS' 'PINES' 'Rejection' 'VPS' 'VPS_nooccip' 'GSR' 'Heart' ...
    'FM-Multisens' 'FM-pain' 'Empathic_Care' 'Empathic_Dist'};

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


function [image_obj, networknames, imagenames] = load_kragel18(varargin)
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



function [image_obj, networknames, imagenames] = load_kragelemotion(varargin)

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


function [image_obj, networknames, imagenames] = load_kragelschemas(varargin)

% Load Kragel 2019 emotion maps
% ------------------------------------------------------------------------

networknames = {'Adoration'	'Aesthetic Appreciation' 'Amusement' 'Anxiety'	'Awe'	'Boredom'	'Confusion'  'Craving'		'Disgust'	'Empathic Pain'	'Entrancement'		'Excitement'	'Fear'	'Horror'	'Interest'	'Joy'	'Romance'	'Sadness'	  'Sexual Desire'	'Surprise'};
s=dir(which('PLS_betas_Adoration.nii.gz'));
all_imgs=dir([s.folder filesep '*.gz']);
imagenames={all_imgs(:).name}';

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns

end % function




function [image_obj, networknames, imagenames] = load_guilt(varargin)

% Load Yu 2019 CerCtx map
% ------------------------------------------------------------------------

networknames = {'Guilt_behavior'};

imagenames = {'Yu_guilt_SVM_sxpo_sxpx_EmotionForwardmask.nii'};

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






function [image_obj, networknames, imagenames] = load_pauli_bg(varargin)

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



function [image_obj, networknames, imagenames] = load_pauli_bg17(varargin)

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


function [image_obj, networknames, imagenames] = load_pauli_subcortical(varargin)

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



function [image_obj, networknames, imagenames] = load_pauli_bg_cortex(varargin)

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



function [image_obj, networknames, imagenames] = load_fibromyalgia(varargin)

% Load Lopez Sola et al. 2017 neural classifier maps
% ------------------------------------------------------------------------
imagenames = {'FM_pain_wholebrain.nii' ...
    'FM_Multisensory_wholebrain.nii' ...
    'rNPS_fdr_pospeaks_smoothed.img' };

networknames = {'FM-pain' 'FM-multisensory' 'NPSp'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function



function [image_obj, networknames, imagenames] = load_csplus(varargin)

% Load Reddan et al. 2018 Neuron CS+ vs. CS- classifier map
% ------------------------------------------------------------------------
imagenames = {'IE_ImEx_Acq_Threat_SVM_nothresh.nii.gz'};

networknames = {'Reddan18CSplus_vs_CSminus'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function



function [image_obj, networknames, imagenames] = load_zhouvps(varargin)

% Load Zhou et al. 2020 eLife generalized vicarious pain signature
% ------------------------------------------------------------------------
imagenames = {'General_vicarious_pain_pattern_unthresholded.nii'};

networknames = {'ZhouVPS'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function



function [image_obj, networknames, imagenames] = load_pain_pdm


% Load Geuter et al. 2018 high-dimensional pain mediator map (PDM)
% ------------------------------------------------------------------------
imagenames = {'Geuter_2020_cPDM_combined_pain_map.nii' 'All_PDM10_unthresholded.nii'};

networknames = {'GeuterPaincPDM','PDM1','PDM2','PDM3','PDM4','PDM5','PDM6','PDM7','PDM8','PDM9','PDM10'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function


% ------------------------------------------------------------------------
% CEKO MPA2 PREDICTIVE MODELS
% ------------------------------------------------------------------------

function [image_obj, networknames, imagenames] = load_mpa2


% Load MPA2 Ceko patterns - multiaversive
% ------------------------------------------------------------------------
imagenames = {'General_bplsF_unthr.nii'
    'Mechanical_bplsF_unthr.nii'
    'Thermal_bplsF_unthr.nii'
    'Sound_bplsF_unthr.nii'
    'Visual_bplsF_unthr.nii'};

networknames = {'General aversive' 'Mech pain' 'Thermal pain' 'Aversive Sound' 'Aversive Visual'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function



% ------------------------------------------------------------------------
% PET TRACERS
% ------------------------------------------------------------------------
function [image_obj, networknames, imagenames] = load_hansen22
datfilename = 'Hansen_2022_PET_tracer_maps_masked.mat';
fullfilename = which(datfilename);    

load(fullfilename)
image_obj=obj_masked_zscored;
imagenames=image_obj.image_names;
imagenames=cellstr(imagenames);
networknames=image_obj.metadata_table(:,1);

networknames = networknames.target'; % format to cell for uniformity with other map sets

end % function


% ------------------------------------------------------------------------
% Stroop
% ------------------------------------------------------------------------
function [image_obj, networknames, imagenames] = load_stroop

imagenames = {'stroop_pattern_wani_121416.nii'}    ;

networknames = {'Stroop'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');
    
end % function


function [image_obj, networknames, imagenames] = load_ncs

% Load NCS craving signatures
% ------------------------------------------------------------------------

networknames = {'NCS' 'NCSdrugs' 'NCSfood'};

imagenames = {'craving_wmapN99_boot10K_02-May-2022.img' ...
    'wmap_onlyDRUGS_l2nGM_N99_20220428.img' ...
    'wmap_onlyFOOD_l2nGM_N99_20220428.img'};    % Koban et al Nat Neurosci 2022

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');  % loads images with spatial basis patterns

end  % function


% ------------------------------------------------------------------------
% WAGER KANG 2015 EMOTION CATEGORY META-ANALYSIS
% ------------------------------------------------------------------------
function [image_obj, networknames, imagenames] = load_emotionmeta

imagenames =     {'Wager_Kang_PlosCB_emometa_2015_anger.nii.gz'  
    'Wager_Kang_PlosCB_emometa_2015_disgust.nii.gz'
    'Wager_Kang_PlosCB_emometa_2015_fear.nii.gz'   
    'Wager_Kang_PlosCB_emometa_2015_happy.nii.gz'  
    'Wager_Kang_PlosCB_emometa_2015_sad.nii.gz'}    ;

networknames = {'Anger' 'Disgust' 'Fear' 'Happy' 'Sad'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function


% ------------------------------------------------------------------------
% Margulies 2016 PNAS principal gradient 1
% ------------------------------------------------------------------------
function [image_obj, networknames, imagenames] = load_marg

imagenames =     {'MNI152NLin2009cAsym_margulies_grad1.nii.gz'};

networknames = {'PrincipalGradient1'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function

% ------------------------------------------------------------------------
% Margulies 2016 PNAS principal gradient 1
% ------------------------------------------------------------------------
function [image_obj, networknames, imagenames] = load_margfsl

imagenames =     {'MNI152NLin6Asym_margulies_grad1.nii.gz'};

networknames = {'PrincipalGradient1'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function


% ------------------------------------------------------------------------
% 2012 Allen Brain Project transcriptomic_gradients
% ------------------------------------------------------------------------
function [image_obj, networknames, imagenames] = load_transcriptomic_gradients

imagenames =     {'transcritptomic_gradients_MNI152NLin6Asym.nii.gz'};

networknames = {'Allen_genePC1' 'Allen_genePC2' 'Allen_genePC3'};

imagenames = check_image_names_get_full_path(imagenames);

image_obj = fmri_data(imagenames, [], 'noverbose');

end % function


% ------------------------------------------------------------------------
% NEUROSYNTH
% ------------------------------------------------------------------------


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

% ------------------------------------------------------------------------
% Neurosynth topic maps - forward inference 
% ------------------------------------------------------------------------
function [image_obj, networknames, imagenames] = load_neurosynth_topics_fi

ns = load(which('neurosynth_data_obj.mat'));

image_obj = ns.topic_obj_forwardinference;

networknames = ns.topic_obj_forwardinference.metadata_table.("Topic name_Summary")';

imagenames = cellstr(image_obj.image_names);

end % function


% ------------------------------------------------------------------------
% Neurosynth topic maps - reverse inference 
% ------------------------------------------------------------------------
function [image_obj, networknames, imagenames] = load_neurosynth_topics_ri

ns = load(which('neurosynth_data_obj.mat'));

image_obj = ns.topic_obj_reverseinference;

networknames = image_obj.metadata_table.("Topic name_Summary")';

imagenames = cellstr(image_obj.image_names);

end % function

            



% ------------------------------------------------------------------------
% BMRK3
% ------------------------------------------------------------------------


function [image_obj, networknames, imagenames] = load_bmrk3 %#ok<*STOUT>

% This code loads a dataset object saved in a mat file, and attempts to
% download it if it cannot be found.

fmri_data_file = which('bmrk3_6levels_pain_dataset.mat');

if isempty(fmri_data_file)
    
    % attempt to download
    disp('Did not find data locally...downloading data file from figshare.com')
    
    fmri_data_file = websave('bmrk3_6levels_pain_dataset.mat', 'https://ndownloader.figshare.com/files/12708989');
    
end

sprintf('Loading: %s\n', fmri_data_file);

load(fmri_data_file); %#ok<*LOAD>

for i = 1:size(image_obj.dat, 2)
    networknames{1, i} = sprintf('Subj%03d_%02dDegreesC', image_obj.additional_info.subject_id(i), image_obj.additional_info.temperatures(i)); %#ok<*AGROW>
end

imagenames = networknames';

descriptives(image_obj);
disp('Temperature data in image_obj.additional_info.temperatures');
disp('Pain ratings in image_obj.Y');

end

% ------------------------------------------------------------------------
% KRAGEL270
% ------------------------------------------------------------------------


function [image_obj, networknames, imagenames] = load_kragel18_alldata

myfile = which('kragel_2018_nat_neurosci_270_subjects_test_images.mat');
if exist(myfile, 'file')
    
    fprintf('Loading %s\n', myfile);
    load(myfile,'data_obj')
    
    image_obj = data_obj;
    networknames = data_obj.additional_info;
    
    if isempty(data_obj.metadata_table)
        error('You have an old file %s\nmetadata_table is empty. Consider removing and rerunning load_image_set to re-download.', myfile);
    end
    
    imagenames = data_obj.metadata_table.imagenames;
    
else
    % Load the files from disk, and clean up downloaded files
    
    fprintf('Did not find %s on path.\nUsing retrieve_neurovault_collection() to download collection 3324\n', 'kragel_2018_nat_neurosci_270_subjects_test_images.mat');
    
    try
        % This is not working for some computers/maybe things have changed
        % with the API. if it fails, try direct download from Dropbox
        
        files_on_disk = retrieve_neurovault_collection(3324);
        data_obj = fmri_data(files_on_disk);
        data_obj = enforce_variable_types(data_obj);
        
        % clean up
        try
            for i = 1:length(files_on_disk), delete(files_on_disk{i}); end
            % remove non-gzipped files just in case
            for i = 1:length(files_on_disk), delete(files_on_disk{i}(1:end-3)); end
            rmdir('3324');
        catch
            disp('Failed to clean up and remove files after download. Check files.')
        end
        
        % resort files/images in order
        labels=regexp(files_on_disk,'Study\d+', 'match');
        labels = cat(1, labels{:});
        labels = strrep(labels, 'Study' ,'');
        for i = 1:length(labels), labels{i} = str2num(labels{i}); end % extract numbers from text
        labels = cat(1, labels{:});
        
        subj=regexp(files_on_disk,'Subject\d+', 'match');
        subj = cat(1, subj{:});
        subj = strrep(subj, 'Subject' ,'');
        for i = 1:length(subj), subj{i} = str2num(subj{i}); end %#ok<*ST2NM> % extract numbers from text
        subj = cat(1, subj{:});
        
        nums = 1000 * labels + subj;
        [~, wh_sort] = sort(nums, 'ascend');
        
        files_on_disk = files_on_disk(wh_sort);
        
        data_obj.dat = data_obj.dat(:, wh_sort);
        data_obj.image_names = data_obj.image_names(wh_sort, :);
        data_obj.fullpath = data_obj.fullpath(wh_sort, :);
        
        % label the images
        sorted_study_labels = labels(wh_sort);
        
        Studynumber = sorted_study_labels;
        Orig_Studynumber = labels;
        
        % imagenames will become text labels
        [imagenames, Domain, Subdomain, StudyCodes] = deal(cell(size(labels)));
        
        networknames = {'ThermalPain1' 'ThermalPain2' 'VisceralPain1' 'VisceralPain2' 'MechanicalPain1' 'MechanicalPain2' ...
            'Cog WM1' 'Cog WM2' 'Cog Inhib1' 'Cog Inhib2' 'Cog RespSel1' 'Cog RespSel2' ...
            'Emotion_Aversiveimages1' 'Emotion_Aversiveimages2' 'Emotion_Rejection1' 'Emotion_VicariousPain2' 'Emotion_AversiveSound1' 'Emotion_AversiveSound2'};
        
        domains = {'Pain' 'Pain' 'Pain' 'Pain' 'Pain' 'Pain' ...
            'Cog_control' 'Cog_control' 'Cog_control' 'Cog_control' 'Cog_control' 'Cog_control' ...
            'Neg_Emotion' 'Neg_Emotion' 'Neg_Emotion' 'Neg_Emotion' 'Neg_Emotion' 'Neg_Emotion'};
        
        subdomains = {'Thermal' 'Thermal' 'Visceral' 'Visceral' 'Mechanical' 'Mechanical' ...
            'WorkingMem' 'WorkingMem' 'Inhibition' 'Inhibition' 'ResponseSelect' 'ResponseSelect' ...
            'Images' 'Images' 'Social' 'Social' 'Sounds' 'Sounds'};
        
        studycodes = {'Atlas_2010_EXP' 'Wager_2013_BMRK3' 'Kano_2017_Rectal' 'Rubio_2015_Rectal' 'Ceko_Woo_MPA1_Mech' 'Ceko_MPA2_Mech' ...
            'DeYoung_2009_WM' 'vanAst_2016_WM' 'Aron_2007_RespSel' 'Xue_2008_RespSel' 'ds101_SimonConflict_NYU' 'Kelly_2008_Flanker' 'Gianaros_2014_IAPS' ...
            'Yarkoni_2011_IAPS' 'Kross_2011_Rejection' 'Krishnan_2016_VicariousPain' 'Losin_Geuter_2018_BMRK5_IADS' 'Kragel_PM01_IADS'};
        
        for i = 1:18
            
            Domain(sorted_study_labels == i) = domains(i);
            Subdomain(sorted_study_labels == i) = subdomains(i);
            imagenames(sorted_study_labels == i) = networknames(i);
            
            StudyCodes(sorted_study_labels == i) = studycodes(i);
        end
        
        data_obj.metadata_table = table(Domain, Subdomain, imagenames, Studynumber, Orig_Studynumber, StudyCodes);
        data_obj.additional_info = networknames;
        
        disp('Downloaded and created object successfully.')
        disp('To save for future use (no re-download), store the object in a')
        disp('variable called data_obj, and save this variable in a file called')
        disp('kragel_2018_nat_neurosci_270_subjects_test_images.mat on your Matlab path.');
        
        % save kragel_2018_nat_neurosci_270_subjects_test_images data_obj
        
        image_obj = data_obj;
        
    catch
        % retrieve_neurovault_collection did not work. Download from
        % Dropbox:
        
        myfile = 'kragel_2018_nat_neurosci_270_subjects_test_images.mat';
        
        disp('Retrieving from Neurovault failed. Trying to download from Google Drive.')
        disp('Saving file with objects: kragel_2018_nat_neurosci_270_subjects_test_images.mat')
        
        %websave(myfile, 'https://www.dropbox.com/s/i88qgg88lgsm0s6/kragel_2018_nat_neurosci_270_subjects_test_images.mat?dl=1');
        
        try
        websave(myfile, 'https://drive.google.com/open?id=1ghDaM55w3dHW2StZ74VTQnJXoE_2uQKm&authuser=tor.d.wager%40dartmouth.edu&usp=drive_fs');
        
        catch
            % figshare
            websave(myfile, 'https://figshare.com/ndownloader/files/42143352')
        end

        fprintf('Loading %s\n', myfile);
        load(myfile)
        
        image_obj = data_obj;
        networknames = data_obj.additional_info;
        imagenames = data_obj.metadata_table.imagenames;
        
    end % try Neurovault...catch
    
end

end % load kragel18_alldata




function table_list = list_signatures



pain =      [1 1 1 1 0 0 0 0 0 0   1 1 0 0 0   0 0 0 0 0 0 0   1 0 0 0   1 0   1 1 1 0 0    0 0 0 0 0 0 0 0 1]';
negemo =    [0 0 0 0 1 1 0 0 0 0   0 0 0 1 1   0 1 0 1 0 1 0   0 0 1 1   0 0   1 0 0 1 1    0 0 0 0 0 0 1 0 0]';
empathy =   [0 0 0 0 0 0 1 1 0 0   0 0 1 1 0   0 0 0 0 0 0 0   0 0 0 0   0 1   0 0 0 0 0    0 0 0 0 0 0 0 0 0]';
physio =    [0 0 0 0 0 0 0 0 1 1   0 0 0 0 0   0 0 0 0 0 0 0   0 0 0 0   0 0   0 0 0 0 0    0 0 0 0 0 0 0 0 0]';
posemo =    [0 0 0 0 0 0 0 0 0 0   0 0 0 0 0   1 0 1 0 0 0 0   0 0 0 0   0 0   0 0 0 0 0    0 0 0 0 0 0 0 1 0]';
cogcontrol =[0 0 0 0 0 0 0 0 0 0   0 0 0 0 0   0 0 0 0 0 0 0   0 1 0 0   0 0   0 0 0 0 0    0 0 1 0 0 0 0 0 0]';
regulation =[0 0 0 0 0 0 0 0 0 0   0 0 0 0 0   0 0 0 0 0 0 0   0 0 0 0   0 0   0 0 0 0 0    1 1 0 0 0 0 0 0 0]';
reward =    [0 0 0 0 0 0 0 0 0 0   0 0 0 0 0   0 0 0 0 0 0 0   0 0 0 0   0 0   0 0 0 0 0    0 0 0 1 1 1 0 1 0]';
other =     [0 0 0 0 0 0 0 0 0 0   0 0 0 0 0   0 0 0 0 1 0 1   0 0 0 0   0 0   0 0 0 0 0    0 0 0 1 1 1 0 0 0]';

keyword = {'NPS' 'NPSpos' 'NPSneg' 'SIIPS' 'PINES' 'Rejection' 'VPS' 'VPS_nooccip' 'GSR' 'Heart' ...
    'FM-Multisens' 'FM-pain' 'Empathic_Care' 'Empathic_Dist' 'Guilt_behavior' ...
    'Amused' 'Angry' 'Content' 'Fearful' 'Neutral' 'Sad' 'Surprised' ...
    'Kragel18Pain' 'Kragel18CogControl' 'Kragel18NegEmotion' 'Reddan18CSplus_vs_CSminus' ...
    'GeuterPaincPDM' 'ZhouVPS' ...
    'General aversive' 'Mech pain' 'Thermal pain' 'Aversive Sound' 'Aversive Visual'  ...
    'PlaceboPvsC_Antic' 'PlaceboPvsC_Pain' 'stroop' 'NCS' 'NCSdrugs' 'NCSfood' 'Painvalue' 'Moneyvalue' 'Shockintensity'}';

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
    'Ashar_2017_empathic_distress_marker.nii' ...       
    'Yu_guilt_SVM_sxpo_sxpx_EmotionForwardmask.nii' ...     % Yu 2019 Cer Ctx Guilt behavior
    'mean_3comp_amused_group_emotion_PLS_beta_BSz_10000it.img' ...  % Kragel 2015 emotion PLS maps
    'mean_3comp_angry_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_content_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_fearful_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_neutral_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_sad_group_emotion_PLS_beta_BSz_10000it.img' ...
    'mean_3comp_surprised_group_emotion_PLS_beta_BSz_10000it.img' ...
    'bPLS_Wholebrain_Pain.nii'  ...                         % Kragel 2018 whole-brain pain cog control neg emotion
    'bPLS_Wholebrain_Cognitive_Control.nii'  ...
    'bPLS_Wholebrain_Negative_Emotion.nii'  ...
    'IE_ImEx_Acq_Threat_SVM_nothresh.nii.gz' ...
    'Geuter_2020_cPDM_combined_pain_map.nii' ...
    'General_vicarious_pain_pattern_unthresholded.nii' ...  % Zhou 2020 eLife vicarious pain
    'General_bplsF_unthr.nii'  ...                          % MPA2 general vs. specific aversiveness
    'Mechanical_bplsF_unthr.nii'  ...
    'Thermal_bplsF_unthr.nii'  ...
    'Sound_bplsF_unthr.nii'  ...
    'Visual_bplsF_unthr.nii'  ...
    'PlaceboPredict_Anticipation.img'   ...                 % Wager 2011 prediction of placebo brain [P - C]->behav [P - C]
    'PlaceboPredict_PainPeriod.img'    ...                  % During pain [P - C]->behav [P - C]
    'stroop_pattern_wani_121416.nii' ...
    'craving_wmapN99_boot10K_02-May-2022.img' ...           % Kober 2022 Craving NCS signature
    'wmap_onlyDRUGS_l2nGM_N99_20220428.img' ...
    'wmap_onlyFOOD_l2nGM_N99_20220428.img' ...
    'painvalue_weights_fdr05.nii.gz' ...                    % Coll 2022 PNAS decision value pain
    'moneyvalue_weights_fdr05.nii.gz' ...
    'shockintensity_weights_fdr05.nii.gz' ...
    }';

table_list = table(keyword, pain, negemo, posemo, empathy, physio, cogcontrol, regulation, reward, other, imagenames);

end % table_list

