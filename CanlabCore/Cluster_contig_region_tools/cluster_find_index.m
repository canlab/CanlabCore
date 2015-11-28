function [wh_cluster, min_distance] = cluster_find_index(cl, varargin)
% Ever see an interesting blob when visualizing a clusters structure, but
% don't know which index number in the clusters structure vector it
% corresponds to?
%
% With this function, find the index number of the closest cluster to one you specify
% graphically by clicking on.
%
% :Usage:
% ::
%
%    function [wh_cluster, min_distance] = cluster_find_index(cl, [keep display flag, 1/0])
%
keepdisplay = 0;
if length(varargin) > 0
    keepdisplay = varargin{1};
end

if ~keepdisplay
    cluster_orthviews(cl, {[1 0 0]});
end

input('Click near the center of the cluster you want to find the index number of and press return.')

centers = cat(1, cl.mm_center);
pos = spm_orthviews('Pos')';

d = sqrt(sum((centers - repmat(pos, length(cl), 1)) .^ 2, 2));

[min_distance, wh_cluster] = min(d);

cluster_orthviews(cl(wh_cluster), {[0 1 0]}, 'add');

disp(['Closest cluster is number ' num2str(wh_cluster(1)), ' with a center of mass ' num2str(min_distance(1)) ' mm away.'])

end
