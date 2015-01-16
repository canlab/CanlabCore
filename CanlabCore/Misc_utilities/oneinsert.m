function yout = oneinsert(removed_voxels, data)

% yout = oneinsert(removed_voxels, data)
%
% Insert ones (for p values) back into data series (i.e., obj.p)
% from which they were removed.
%
% removed_voxels is indicator for the locations of removed voxels
%
% This function is a modified version of naninsert.m, but different.
%
% Wani Woo, Aug 2012
%
% See naninsert.m, nanremove.m, zeroinsert.m

yout = zeros(size(removed_voxels,1), size(data,2));
yout(removed_voxels,:) = 1;
yout(~removed_voxels,:) = data;

return
