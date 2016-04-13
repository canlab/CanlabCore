function [dat_out,clindx,keepit] = iimg_cluster_prune(dat,datsig,volInfo)
% Prunes an image data vector (dat) assumed to contain suprathreshold
% contiguous clusters ("blobs") by saving only those blobs that have one or
% more significant (nonzero) elements in another image, datsig
%
% :Usage:
% ::
%
%     [dat_out,clindx,keepit] = iimg_cluster_prune(dat,datsig,volInfo)
%
% One intended use is to define an FWE-corrected significance map in
% datsig, and report blobs at some lower threshold (in dat) that have at
% least one corrected voxel.
%
% Try iimg_threshold to create dat vectors from Analyze images (e.g.,
% statistic images)
%
% ..
%    tor wager, july 06
% ..

% check data type and reduce to in-mask only if necessary
[dattype,dat] = iimg_check_indx(dat,volInfo,'masked');
[dattype2,datsig] = iimg_check_indx(datsig,volInfo,'masked');

% output in same format as dat input
switch dattype
    case 'full'
        dat_out = zeros(volInfo.nvox,1);          % in full image space
    case 'masked'
        dat_out = zeros(volInfo.n_inmask,1);      % in-mask space
    otherwise
        error('Internal error.  dattype should have been checked earlier.')
end

clindx = iimg_cluster_index(dat,volInfo.xyzlist');

n = max(clindx);
keepit = zeros(1,n);

for i = 1:n
    wh = clindx == i;
    keepit(i) = any(datsig(wh));

    if keepit(i)
        switch dattype
            case 'full'
                dat_out(volInfo.wh_inmask(wh)) = dat(wh);
            case 'masked'
                dat_out(wh) = dat(wh);
        end
    end
end

return
