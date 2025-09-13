function [obj, wh_order] = reorder_atlas_regions(obj, wh_order, varargin)
% reorder_atlas_regions Reorder atlas region indices based on label groups.
%
% :Usage:
% ::
%     [obj, wh_order] = reorder_atlas_regions(obj, wh_order)
%     [obj, wh_order] = reorder_atlas_regions(obj, [], 'labelgroups', labelgroups, 'compare_property', 'labels_5')
%     [obj, wh_order] = reorder_atlas_regions(atlas_obj, [10:-1:1])  % only a subset of the first 10, in reverse order
%
% :Inputs:
%
%   **obj:**  
%        An atlas object containing region labels in one or more properties (e.g., 
%        'labels', 'labels_1', ... 'labels_5').
%
%   **wh_order:**  
%        An n×1 numeric vector indicating the current ordering of atlas regions.
%
%   **'labelgroups':** [1×k cell array of strings] (Optional)  
%        A cell array of label group names in the desired order. Regions whose labels 
%        match these strings (using the function contains) will be reordered accordingly.
%
%   **'compare_property':** [string] (Optional)  
%        The name of the property in the atlas object to compare against labelgroups 
%        (e.g., 'labels'). This property must contain the region labels.
%
% :Outputs:
%
%   **wh_order:**  
%        An n×1 numeric vector of indices representing the new ordering of atlas regions.
%        If the optional inputs 'labelgroups' and 'compare_property' are provided, the input 
%        wh_order is replaced by the ordering derived from reorder_by_labelgroups; otherwise,
%        the original wh_order is returned.
%
% :Examples:
% ::
%     % Suppose atlas_obj is an atlas object with a property 'labels', and an initial ordering wh_order.
%     labelgroups = {'Frontal', 'Parietal', 'Temporal', 'Occipital'};
%     wh_order = reorder_atlas_regions(atlas_obj, wh_order, 'labelgroups', labelgroups, 'compare_property', 'labels');
%
%      atlas_obj = load_atlas('canlab2024');
%      % The group 'labels_5' was added using this script:
%      % canlab2024_add_labels_5_networks_and_large_structures
%      partitionlabels = unique(brainpathway_obj.region_atlas.labels_5);
%      partitionlabels = partitionlabels([5:37 1 38 39 2:4]); % Reorder groups
%      [atlas_obj_reordered, wh_order] = reorder_atlas_regions(atlas_obj, [], 'labelgroups', partitionlabels, 'compare_property', 'labels_5');
%
% :See also:
%   reorder_by_labelgroups
%

% Programmers' notes:
% Created by Tor Wager
% 12/2019: Tor fixed bug in reordering; probability_maps field was not
% preivously being considered, which resulted in wrong ordering when
% present. .dat is rebuilt from probability_maps if available.
% 3/2025 Updated help and added labelgroups option

p = inputParser;
addParameter(p, 'labelgroups', {}, @(x) iscell(x));
addParameter(p, 'compare_property', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
labelgroups = p.Results.labelgroups;
compare_property = p.Results.compare_property;

if ~isempty(labelgroups) && ~isempty(compare_property)

    labels_to_reorder = obj.(compare_property);
    wh_order = reorder_by_labelgroups(labels_to_reorder, labelgroups);

end



n = num_regions(obj);

if any(wh_order > n)
    error('wh_order contains integer values that do not have corresponding region indices.');
end


newdat = int32(zeros(size(obj.dat)));

for i = 1:length(wh_order)
    
    newdat(obj.dat == wh_order(i)) = i;
    
end

obj.dat = newdat;

% Reorder probability maps
if ~isempty(obj.probability_maps)
    
    if ~size(obj.probability_maps, 2) == n
        error('obj.probabilty_maps is not empty, but is the wrong size! illegal object.');
    end
    
    obj.probability_maps = obj.probability_maps(:, wh_order);
    
end

% Reorder labels
myfields = {'labels' 'labels_2' 'labels_3' 'labels_4' 'labels_5' 'property_descriptions', 'label_descriptions'};

for i = 1:length(myfields)
    
    if ~isempty(obj.(myfields{i})) && length(obj.(myfields{i})) == n
        
        obj.(myfields{i}) = obj.(myfields{i})(wh_order);
        
    end
    
end

end % function





function new_order_indx = reorder_by_labelgroups(labels_to_reorder, labelgroups)
% Assume sampleLabels is an n×1 cell array of label strings,
% and partitionlabels is a cell array of group names in the desired order.
%
% This code finds, for each group in partitionlabels, the indices of
% sampleLabels that contain the group string (preserving their original order),
% and concatenates them into a single n×1 vector "new_order_indx".

% labels_to_reorder
% labelgroups, e.g., partitionlabels

new_order_indx = [];  % Initialize an empty vector to hold the new ordering

for i = 1:length(labelgroups)
    currentGroup = labelgroups{i};              % Get the current group label
    groupIdx = find(contains(labels_to_reorder, currentGroup))';  % Find indices of sampleLabels belonging to this group
    new_order_indx = [new_order_indx; groupIdx];                  % Append these indices in order
end

% new_order_indx now is an n×1 vector with the indices reordered by labelgroups.

end % function
