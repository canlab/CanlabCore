function [parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg] = extract_data(obj, data_obj, varargin)
% Atlas object-class method that extracts and averages data stored in an
% fmri_data object [data_obj] basd on atlas regions.
%
% - Runs image_vector.apply_parcellation
% - Mask/Atlas image does not have to be in the same space as the images to extract from.
%   Resamples the atlas object to be in the same space as the data object if needed.
% - For low-resolution data images, some parcels may be lost. This will be reported as output
% - For speed (large N), you may want to pre-resample the atlas obj to match the functional data space
%
% - If optional weight_map obj is entered, will apply local pattern weights.
%
% :Usage:
% ::
%
%    data_table = extract_data(atlas_obj, data_obj, [weight_map_obj])
%
% :Inputs:
%
%   - obj:      atlas object 
%   - data_obj: fmri_data object
%
% :Optional inputs:
%
%   **see help for image_vector.apply_parcellation.
%
%   **pattern_expression:**
%       followed by an fmri_data object with multivariate pattern to apply.
%       local patterns in each parcel are applied, and pattern responses
%       returned for each parcel.
%
%   **correlation:**
%        calculate the pearson correlation coefficient of each
%        image in dat and the values in the mask.
%
%   **norm_mask:**
%        normalize the mask weights by L2 norm
%
%   **ignore_missing:**
%        used with pattern expression only. Ignore weights on voxels
%        with zero values in test image. If this is not entered, the function will
%        check for these values and give a warning.
%
%   **cosine_similarity:**
%        used with pattern expression only. scales expression by product of
%        l2 norms (norm(mask)*norm(dat))
%
% :Outputs:
%
%   **parcel_means:**
%       Matrix of mean data values for each parcel, data images x parcels
%       This will always be the parcel means, even if pattern expression is
%       requested. Pattern expression will be returned in a separate
%       output.
%
%   **parcel_pattern_expression:**
%        If 'pattern_expression' is entered followed by a linear
%        multivariate pattern object, parcel_pattern_expression returns the 
%        local expression (dot product) of the pattern for each data image
%        (row), in each parcel region (column)
%
%   **parcel_valence:**
%        A measure of overall activation/deactivation, relative to the region average
%        Define pattern valence as similarity with unit vector (all positive,
%        identical weights)
%        If x is a vector of pattern weights,
%        The cosine similarity with the unit vector is a measure of how uniform
%        the weights are, on a scale of 1 to -1.  1 indicates that the pattern
%        computes the region average (all weights identical and positive). -1
%        indicates that the pattern computes the negative region average (all
%        weights identical and negative).
%
%   **rmsv_pos, rmsv_neg:**
%         Root-mean-square values
%         Signed root mean square values for each image (rows) in each parcel
%         (columns), for positively-valued and negatively-valued voxels,
%         respectively. This is useful for visualizing and assessing pattern
%         weights, i.e., when the input data image is a weight pattern map.
%         It can tell you how strong the contribution of each region is to
%         prediction (high weight values) and how homogenous it is (pos vs.
%         negative energy/rmsvs). 
%         RMSVs are expressed in weights / cm^3 of brain tissue, and volume is
%         regularized by adding 1 cm^3 to avoid unstable results with small
%         regions.
%
% :Examples:
% ::
%
% Extract region averages from an atlas of <500 regions spanning the brain
% ------------------------------------------------------------------------------
% help load_atlas                                       % get a list of atlases readily available
% atlas_obj = load_atlas('canlab2018_2mm');             % Load combined atlas
% test_images = load_image_set('emotionreg');           % Load some test images
% parcel_means = extract_data(atlas_obj, test_images);  % Extract region average data
%
% Also consider:
% DAT = extract_measures_batch(test_images);            % Extract a broad set of data measures, including the regions above
%
%
% Apply local NPS pattern weights to each of a series of atlas regions to
% get local pattern response
% ------------------------------------------------------------------------------
%   [nps, npsnames] = load_image_set('npsplus');
%   siips = get_wh_image(nps, 4);
%   nps = get_wh_image(nps, 1);
%
%    load pain_pathways_atlas_obj.mat
%    test_data_obj = load_image_set('pain');
%    [parcel_means, parcel_pattern_expression, parcel_valence] = extract_data(pain_pathways_finegrained, test_data_obj, 'pattern_expression', nps, 'cosine_similarity');
%
%
% :See also:
%
% For an non-object-oriented alternative, see extract_image_data.m
% fmri_data.extract_roi_averages, apply_parcellation.m

% Programmers' notes:
% This function is different from fmri_data.extract_roi_averages
% Better to have only one function of record in the future...
% Note: 
% cl = extract_roi_averages(imgs, atlas_obj{1});
% accomplishes the same task as apply_parcellation, returns slightly different values due to interpolation.  

% --------------------------------------------
% Resample data to the atlas space if needed
% --------------------------------------------

isdiff = compare_space(data_obj, obj);

if isdiff == 0 || isdiff == 3 % Same space, diff voxels
   % do nothing
elseif isdiff == 2
    error('Invalid object: volInfo structure missing for one or more objects');
else
    % Resample to atlas. 
    disp('Resampling data object to space defined by atlas object.');
    data_obj = resample_space(data_obj, obj);
end

[parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg] = deal([]);

[parcel_means, parcel_pattern_expression, parcel_valence, rmsv_pos, rmsv_neg] = apply_parcellation(data_obj, obj, varargin{:});

end % function


