function b = reorder_regions(b, varargin)
% reorder_regions Reorder atlas regions in a brainpathway object.
%
% :Usage:
% ::
%     b = reorder_regions(b)
%     b = reorder_regions(b, 'wh_order', index_vector)
%     b = reorder_regions(b, 'node_clusters', true)
%     b = reorder_regions(b, 'labelgroups', labelgroups, 'compare_property', prop)
%
% :Inputs:
%
%   **b:**  
%        A brainpathway object containing atlas region data. Key properties include:
%           - region_atlas: An atlas object with region labels (e.g., in 'labels', 'labels_1', etc.).
%           - region_dat: A matrix with region data (columns corresponding to regions).
%           - node_clusters: A numeric vector assigning cluster labels to nodes.
%           - node_labels: A cell array of node label strings.
%
%   **'wh_order':** [numeric vector] (Optional)  
%        An arbitrary integer index vector specifying the desired new ordering of regions.
%
%   **'node_clusters':** [logical] (Optional)  
%        If true, the regions will be sorted in ascending order based on the node_clusters property.
%
%   **'labelgroups':** [1×k cell array of strings] (Optional)  
%        A cell array of label group names specifying the desired order. Regions whose labels 
%        (from a specified property) match these strings (using contains) will be grouped 
%        accordingly.
%
%   **'compare_property':** [string] (Optional)  
%        The name of the property in the atlas object to use for comparison with labelgroups 
%        (e.g., 'labels', 'labels_2', etc.). This property must contain the region labels.
%
% :Outputs:
%
%   **b:**  
%        The brainpathway object with its region-related properties reordered. The following 
%        properties are updated:
%           - region_atlas: The atlas object is reordered based on the new index order.
%           - region_dat: Columns are reordered according to the new ordering.
%           - node_labels: If present, reordered to match the new region order.
%           - node_clusters: Reordered accordingly.
%
% :Examples:
% ::
%     % Default: reorder by node_clusters.
%     b = reorder_regions(b, 'node_clusters', true);
%
%     % Custom ordering:
%     customOrder = [5 3 1 2 4];
%     b = reorder_regions(b, 'wh_order', customOrder);
%
%     % Reorder based on label groups using the 'labels' property:
%     labelgroups = {'Frontal', 'Parietal', 'Temporal', 'Occipital'};
%     b = reorder_regions(b, 'labelgroups', labelgroups, 'compare_property', 'labels');
%
% :See also:
%   reorder_atlas_regions, reorder_by_labelgroups
%
% NOTE: May have to check for and reorder more properties as the code is
% developed.

% Author: Tor Wager
% Date: March 2025
% License: GNU General Public License v3 or later

    % Parse optional inputs.
    p = inputParser;
    addParameter(p, 'wh_order', [], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'node_clusters', false, @(x) islogical(x));
    addParameter(p, 'labelgroups', {}, @(x) iscell(x));
    addParameter(p, 'compare_property', '', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    
    wh_order = p.Results.wh_order;
    sortByNodeClusters = p.Results.node_clusters;
    labelgroups = p.Results.labelgroups;
    compare_property = p.Results.compare_property;
    
    % Determine the new ordering.
    if ~isempty(wh_order)

        indx = wh_order;

        % Reorder the atlas object.
        b.region_atlas = reorder_atlas_regions(b.region_atlas, indx);

    elseif sortByNodeClusters && ~isempty(b.node_clusters)
        [~, indx] = sort(b.node_clusters, 'ascend');

        % Reorder the atlas object.
        b.region_atlas = reorder_atlas_regions(b.region_atlas, indx);

    elseif ~isempty(labelgroups) && ~isempty(compare_property)

        % Use the reorder_atlas_regions function on the region_atlas property.
        [b.region_atlas, indx] = reorder_atlas_regions(b.region_atlas, [], 'labelgroups', labelgroups, 'compare_property', compare_property);

    else
        % do nothing
        % indx = (1:numel(b.region_atlas.labels))';
        return

    end

    % Reorder region data columns.
    b.region_dat = b.region_dat(:, indx);
    
    % If node_labels exists and matches the number of regions, reorder it.
    if isfield(b, 'node_labels') && ~isempty(b.node_labels) && numel(b.node_labels) == numel(indx)

        b.node_labels = b.node_labels(indx);

    end
    
    % Similarly, reorder node_clusters if present.

    if isfield(b, 'node_clusters') && ~isempty(b.node_clusters)

        b.node_clusters = b.node_clusters(indx);

    end

end % main function