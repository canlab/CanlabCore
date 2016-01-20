function [cl1,cl2,dat1,dat2] = iimg_cluster_intersect(dat1,dat2,volInfo)
% Prunes two image data vectors (dat1 and dat2) assumed to contain suprathreshold
% contiguous clusters ("blobs") by saving only those blobs that have one or
% more significant (nonzero) elements in both images
%
% :Usage:
% ::
%
%     [cl1,cl2,dat1,dat2] = iimg_cluster_intersect(dat1,dat2,volInfo)
%
% :Inputs:
%
%   **volInfo.xyzlist:**
%        must contain xyz coordinates corresponding to dat1 and dat2
%
%   **dat1 and dat2:**
%        can be either image-length or in-mask length vectorized images
%
% Try iimg_threshold to create dat vectors from Analyze images (e.g.,
% statistic images)
%
% The outputs cl1 and cl2 are the overlap (intersection) clusters, in the
% space of dat1 (cl1) and dat2 (cl2).  The outputs dat1 and dat2 are the
% 'pruned' intersection data vectors.
%
% ..
%    tor wager, july 06
% ..

% prune dat1 with sig values in dat2
[dat1] = iimg_cluster_prune(dat1,dat2,volInfo);

% do the reverse
[dat2] = iimg_cluster_prune(dat2,dat1,volInfo);

% make clusters
cl1 = iimg_indx2clusters(dat1,volInfo);
cl2 = iimg_indx2clusters(dat2,volInfo);

return
