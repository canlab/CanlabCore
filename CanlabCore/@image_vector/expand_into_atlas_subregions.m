% expand_into_atlas_subregions Expand a single-image data object into k images, one per atlas region.
%
% Take an image object (one image) and an integer-valued parcellation /
% atlas whose values indicate k regions, and expand the image object into
% a set of k images with the original values only for voxels in the kth
% region, and zero elsewhere. Useful for breaking a pattern mask into
% local atlas-defined subregions.
%
% This is different from fmri_data.apply_parcellation, which returns
% averages or weighted averages per image for a set of images in a data
% object. apply_parcellation can be used to return local pattern
% responses for any parcellation, but it does not return vectors with
% the local weights themselves.
%
% NOTE: This file is a script-style template (no function declaration).
% It currently expects the variable 'nps' to exist in the workspace as
% the target image object, and uses BN_Atlas_274_noCb_uint16.nii as the
% parcellation. Adapt as needed before running.
%
% :Usage:
% ::
%
%     % Set up nps (an fmri_data / image_vector object), then run:
%     expand_into_atlas_subregions
%
% :Inputs (workspace variables):
%
%   **nps:**
%        Single-image image_vector / fmri_data object holding the
%        pattern to be split by region.
%
% :Outputs (workspace variables):
%
%   **target_obj:**
%        The expanded object with one column per non-empty atlas region,
%        each containing the original pattern values restricted to that
%        region (zero elsewhere). target_obj.additional_info holds the
%        retained region names.
%
%   **r:**
%        A region object built from the expanded data, with one element
%        per contiguous region. region_names tracks the parent atlas
%        label for each region in r.
%
% :See also:
%   - fmri_data.apply_parcellation
%   - select_voxels_by_value
%   - region
%   - canlab_pattern_similarity
%
% ..
%    Programmers' notes:
%    fmri_data.subs
% ..



% define target image obj to sample to : nps

nps = replace_empty(nps); % just in case
target_obj = nps;  % this will be expanded to set of images

orig_image_vector = target_obj.dat(:, 1);


parcellation_file = which('BN_Atlas_274_noCb_uint16.nii');
parcel_obj = fmri_data(parcellation_file);

parcel_obj = resample_space(parcel_obj, nps);

names_file = fullfile(fileparts(parcellation_file), 'cluster_names.mat');
clear names
load(names_file, 'names');  
names;      % returns names var. make sure: break if not

% Find which names match
%wh = ~cellfun(@isempty, strfind(names, string_to_find));

% Index values to look for
% indx = find(wh)';

%% Do expansion
% -----------------------------------------------------------------------
for indx = 1:length(names)
    
    output_obj = select_voxels_by_value(parcel_obj, indx);
        
    target_subregion_vector = orig_image_vector; % will become image values only in region, 0 elsewhere
    target_subregion_vector(output_obj.removed_voxels) = 0;
    
    target_obj.dat(:, indx) = target_subregion_vector;
    
    
end

% Add image for everything not already covered
% *******


% Remove empty images, attach names to object

empty_imgs = sum(abs(target_obj.dat)) < 10*eps;
target_obj.dat(:, empty_imgs) = [];

names(empty_imgs) = [];
target_obj.additional_info = names;


%% Get region object, r
% May have different length (n regions) because of contiguity requirement
% -----------------------------------------------------------------------

dat_wh = get_wh_image(target_obj, 1); 
r = region(dat_wh);

region_names = repmat(names(1), 1, length(r));

for i = 2:size(target_obj.dat, 2)
    
    dat_wh = get_wh_image(target_obj, i); 
    
    % threshold at 20% probability, get rid of stray voxels
    %dat_wh = threshold(dat_wh, [.2 Inf], 'raw-between', 'k', 3);

    r_wh = region(dat_wh);
    
    for j = 1:length(r_wh)
        [r_wh(j).title, r_wh(j).shorttitle] = deal(names{i});
    end
    
    r = [r r_wh];
    
    region_names = [region_names repmat(names(i), 1, length(r_wh))];
    
end
%%    

% similarity_output = canlab_pattern_similarity(dat, pattern_weights, varargin)

