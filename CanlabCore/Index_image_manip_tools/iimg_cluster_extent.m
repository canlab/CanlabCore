function [dat,nvox] = iimg_cluster_extent(dat,volInfo,k)
% Apply a cluster size threshold to a series of index image data vectors
% (each img is one column)
% given volInfo (must be created with iimg_read_img, with extended output)
%
% :Usage:
% ::
%
%     [dat,nvox] = iimg_cluster_extent(dat,volInfo,k)
%
% :Inputs:
%
%   **dat:**
%        may be an indexed image of all image values or only those in mask
%        defined by volInfo


if isempty(k) || k < 2, nvox = []; return, end % do nothing if extent threshold

if size(dat,1) == volInfo.nvox
    dattype = 'full';
elseif size(dat,1) == volInfo.n_inmask
    dattype = 'masked';
else error('data vector does not match size of image in volInfo!')
end

switch dattype
    case 'full', [clindx,nvox] = iimg_cluster_index(dat(volInfo.wh_inmask, :), volInfo.xyzlist', k);
    case 'masked', [clindx,nvox] = iimg_cluster_index(dat,volInfo.xyzlist',k);
end


switch dattype
    case 'full'
        for i = 1:size(dat,2)
            clindx2 = zeros(volInfo.nvox,1);    % put in original index dims
            clindx2(volInfo.wh_inmask) = clindx(:,i);
            dat(:,i) = dat(:,i) .* (clindx2 > 0);
        end

    case 'masked'
        dat = dat .* (clindx > 0);
end

return
