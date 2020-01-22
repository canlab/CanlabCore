function b = reorder_regions_by_node_cluster(b)
% b = reorder_regions_by_node_cluster(b)

[b.node_clusters, indx] = sort(b.node_clusters, 'ascend');

b.region_atlas = reorder_atlas_regions(b.region_atlas, indx);

b.region_dat = b.region_dat(:, indx);

n = length(b.node_clusters);

if length(b.node_labels) == n
    b.node_labels = b.node_labels(indx);
end

end

