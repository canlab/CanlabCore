function [close_enough,outside_range,nearest_distance,closest_cluster] = cluster_close_enough(cl_match_to,cl_match,mind)
% Finds whether each cluster center in cl_match is within mind mm of a cluster
% center in cl_match_to.
%
% Useful for selecting a list of clusters that are not close to another
% list to, e.g., make a table of.  or this could be used to find clusters
% in a set of correlated clusters that are close to centers in activated
% clsuters.

centers1 = cat(1,cl_match_to.mm_center);
centers2 = cat(1,cl_match.mm_center);

nmatchto = size(centers1,1);
nmatch = size(centers2,1);

% get matrix of distances, centers1 x centers2
d = pdist([centers1; centers2]);
d = squareform(d);
d = d(1:nmatchto,nmatchto+1:end);

[nearest_distance,closest_cluster] = min(d);


d = d < mind;

close_enough = any(d);
outside_range = ~close_enough;

return

