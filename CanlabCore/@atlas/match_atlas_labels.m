function [matchingLabels, diceMatrix, maxdice, has_match, atlas_subset_match] = match_atlas_labels(obj, other_atlas, varargin)
% match_atlas_labels Calculate Dice coefficients between atlas regions and match labels
%
% :Usage:
% ::
%     [matchingLabels, diceMatrix, maxdice, has_match, atlas_subset_match] = ...
%         atlas_obj.match_atlas_labels(other_atlas, 'dice_threshold', [0 - 1 scalar value])
%
% :Inputs:
%
%   **obj:**
%        An atlas object. Its .dat field is an array with integer-coded region
%        assignments (nonzero for regions, 0 for background), and its .labels field
%        is a cell array of strings describing each region.
%
%   **other_atlas:**
%        A second atlas object with the same structure as obj (i.e., with .dat and 
%        .labels fields). Regions are also integer-coded.
%
%   **'dice_threshold':** [scalar] (Optional)
%        A threshold for the Dice coefficient. For each region in obj, if the maximum 
%        Dice coefficient across regions in other_atlas is not above this threshold, the 
%        corresponding label in matchingLabels is set to 'nomatch'. Default is 0.
%
% :Outputs:
%
%   **matchingLabels:**  
%        A cell array of strings (length = number of regions in obj) where each element 
%        is the label (from other_atlas.labels) of the region in other_atlas that best 
%        matches the corresponding region in obj. Regions with maximum Dice coefficient 
%        below dice_threshold are labeled 'nomatch'.
%
%   **diceMatrix:**  
%        An [n x k] numeric matrix of Dice coefficients, where n is the number of regions 
%        in obj (length(obj.labels)) and k is the number of regions in other_atlas 
%        (length(other_atlas.labels)).
%
%   **maxdice:**  
%        A numeric vector (n x 1) containing the maximum Dice coefficient for each region in obj.
%
%   **has_match:**  
%        A logical vector (n x 1) where each element is true if the corresponding region 
%        in obj has a maximum Dice coefficient above dice_threshold.
%
%   **atlas_subset_match:**  
%        A subset of the atlas (obtained using select_atlas_subset) containing only those 
%        regions from obj that have a match in other_atlas (i.e., where has_match is true).
%
% :Examples:
% ::
%
%     % This example shows how to build a new label field in the canlab2024
%     atlas with large structures, including Schafer/Yeo resting-state networks for cortical regions.
%
%     % Load two atlas objects:
%     atlas_obj = load_atlas('canlab2024');
%     yeo = load_atlas('yeo17networks');
%
%     % Compute the Dice coefficient matrix and matching labels with a dice threshold of 0:
%     % The threshold 0.01 labels all the cortical regions in this example, but also a few others. 
%     [matchingLabels, diceMatrix, maxdice, has_match, atlas_subset_match] = ...
%     match_atlas_labels(atlas_obj, yeo, 'dice_threshold', 0.01);
% 
%     % Make a contingency table of the match between labels and cortical regions
%     contingencyTable = @(vec1, vec2)  [sum(vec1 & vec2), sum(vec1 & ~vec2); sum(~vec1 & vec2), sum(~vec1 & ~vec2)];
%     logicalVector = contains(atlas_obj.labels, 'Ctx')';
%     contingencyTable(logicalVector, has_match)
% 
%     % mask out matches for non-cortical regions
%     matchingLabels(~logicalVector) = {'nomatch'};
%
%     % fill in missing labels from other large-scale structures
%     for labels = {'Thal' 'Cblm' 'BG' 'hypothalamus' 'AMY' 'BStem' 'MTL_Hipp'}
%       logicalVector = contains(atlas_obj.labels, labels{1})';
%       matchingLabels(logicalVector) = {labels{1}};
%     end
% 
%     % check for any non-matched labels
%     atlas_obj.labels(contains(matchingLabels, 'nomatch'))'
% 
%     % assign matchingLabels to labels_5 in atlas_obj 
%     atlas_obj.labels_5 = matchingLabels';  % should be row vector
%
% :References:
%   Dice, L.R. (1945). Measures of the amount of ecologic association between species.
%   Ecology, 26(3), 297-302.
%
% :See also:
%   condf2indic, select_atlas_subset
%
% Author: Tor Wager
% Date: March 2025
% License: GNU General Public License v3 or later

    % Parse optional input: dice_threshold (default = 0)
    p = inputParser;
    addParameter(p, 'dice_threshold', 0, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});
    dice_threshold = p.Results.dice_threshold;

    % Resample the other atlas to match the space of obj (if necessary)
    other_atlas = resample_space(other_atlas, obj);

    % Get unique region codes for each atlas (exclude background code 0)
    regions_obj = unique(obj.dat);
    regions_obj(regions_obj == 0) = [];
    n = numel(regions_obj);
    
    regions_other = unique(other_atlas.dat);
    regions_other(regions_other == 0) = [];
    k = numel(regions_other);
    
    % Preallocate matrix for Dice coefficients.
    diceMatrix = zeros(n, k);
    
    % Calculate Dice coefficients between each region in obj and each region in other_atlas.
    for i = 1:n
        
        mask1 = (obj.dat == regions_obj(i));

        for j = 1:k
            mask2 = (other_atlas.dat == regions_other(j));
            intersection = sum(mask1 & mask2);
            diceMatrix(i,j) = (2 * intersection) / (sum(mask1) + sum(mask2));
        end
    end
    
    % For each region in obj, find the region in other_atlas with the highest Dice coefficient.
    matchingLabels = cell(n, 1);
    for i = 1:n
        [~, bestIdx] = max(diceMatrix(i,:));
        matchingLabels{i} = other_atlas.labels{bestIdx};
    end
    
    % Calculate maximum Dice coefficient per region and determine matches based on threshold.
    maxdice = max(diceMatrix, [], 2);
    has_match = maxdice > dice_threshold;
    matchingLabels(maxdice <= dice_threshold) = {'nomatch'};
    
    % Create atlas_subset_match: subset of obj that has a match in other_atlas.
    atlas_subset_match = select_atlas_subset(obj, find(has_match));
    
end