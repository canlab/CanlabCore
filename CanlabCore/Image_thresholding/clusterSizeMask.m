function [mask,numClusters,XYZ] = clusterSizeMask(sizeThresh,height_mask)
% :Usage:
% ::
%
%     function [mask,numClusters,XYZ] = clusterSizeMask(sizeThresh,height_mask)
%
% ..
%    Tor Wager, 10/27/01
% ..

mask = []; numClusters = 0;, XYZ = [];

% Get point list of activated voxels
% ------------------------------------------------------
voxels = mask2voxel(height_mask);	% returns n x 3
voxels = voxels';					% put xyz in a single column

if isempty(voxels)
	disp('No voxels meet height threshold.')
	mask = zeros(size(height_mask));
	return
end

% Get cluster indices of voxels
% ------------------------------------------------------
[cl_index] = spm_clusters(voxels);


% Find index of voxels of sufficient size
% ------------------------------------------------------
if ~isempty(sizeThresh) & sizeThresh > 0
    for i = 1:max(cl_index)
	    a(cl_index == i) = sum(cl_index == i);
    end
else
    sizeThresh = 0;
    a = ones(size(cl_index));
end

which_vox = (a >= sizeThresh);
numClusters = sum(length(unique(cl_index(find(a >= sizeThresh)))));

if numClusters == 0
	disp('No clusters meet extent threshold.')
	mask = zeros(size(height_mask));
    XYZ = [];
	return
end

voxels = voxels(:,which_vox);


% Voxel point list output - 3 vector
% ------------------------------------------------------
XYZ = voxels;


% Convert back to mask for mask output
% ------------------------------------------------------
voxels = voxels';					% convert back to row vectors
mask = voxel2mask(voxels,size(height_mask));

return
