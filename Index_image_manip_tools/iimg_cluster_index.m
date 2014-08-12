function [clindx,nvox] = iimg_cluster_index(dat,xyz,k)
% [clindx,nvox] = iimg_cluster_index(dat,xyz,[k])
% xyz is 3 x n list of voxel coords, volInfo.xyzlist'
% dat is index vector of image values
%
% Returns: cluster index and cluster sizes for non-zero, non-nan voxels

%% initialize
nimgs = size(dat,2);
clindx = zeros(size(dat));
wh_good_data = cell(1,nimgs);
nvox = cell(1,nimgs);

%% Get index values
for i = 1:size(dat,2)   % for each image index

    % return cluster index for each voxel, or zero if no data
    [clindx(:,i),wh_good_data{i}] = cluster_index(dat(:,i),xyz);

    % CLuster extent, if entered
    if nargin > 2 && k > 1
        [nvox{i},clindx(:,i),wh_omit] = iimg_cluster_sizes(clindx(:,i),clindx(:,i),k);
    else
        nvox{i} = iimg_cluster_sizes(clindx(:,i));
    end

end

%%

return






%% eliminate zero or NaN data values
% And make sure we get spm cluster indices
function [clindx,wh_good_data] = cluster_index(dat,xyz)

clindx = zeros(size(dat));
wh_good_data = find(dat & ~isnan(dat));

if isempty(wh_good_data), return, end

xyz = xyz(:,wh_good_data);  % in-image coordinates

clusters = get_cluster_index(xyz);
clindx(wh_good_data) = clusters;


%% use spm clusters
function clusters = get_cluster_index(xyz)
nvox = size(xyz,2);
if nvox == 0, clusters = []; return, end
if nvox < 50000
    clusters = spm_clusters(xyz)';
else
    disp('Too many voxels for spm_cluster. Blobs will not be correctly broken up into contiguous clusters.')
    clusters = ones(nvox,1);
end
return


%% Count voxels in each contiguous cluster and threshold, if additional
% args are entered
function [nvox,dat,wh_omit] = iimg_cluster_sizes(clindx,dat,k)

nclust = max(clindx);
nvox = zeros(1,nclust);
if nargin > 2
    wh_omit = false(size(dat,1),1);
end

for i = 1:nclust
    wh = find(clindx == i);
    nvox(i) = length(wh);

    if nargin > 2
        if nvox(i) < k
            wh_omit(wh) = 1;
        end
    end
end

if nargin > 2
    nvox(nvox < k) = [];
    dat(wh_omit) = 0;
end

return
