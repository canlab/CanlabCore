function [imgvec,maskvec] = iimg_clusters2indx(cl,volInfo)
% Take a clusters structure and turn it into an indexed image of dims
% volInfo.dim
%
% :Usage:
% ::
%
%     [imgvec,maskvec] = iimg_clusters2indx(cl,volInfo)
%     [imgvec,maskvec] = iimg_clusters2indx(cl,my_image_name)
%
% :Outputs:
%
%   **imgvec:**
%        vector of all image voxels
%
%   **maskvec:**
%        vector of in-mask voxels
%
% Uses:  volInfo.nvox, .wh_inmask, .dim

% convert image to volinfo struct, if necessary
if isstr(volInfo)
    volInfo = iimg_read_img(volInfo,2);
end

imgvec = false(volInfo.nvox,1);

if isempty(cl)
    maskvec = false(volInfo.n_inmask,1);
    return
end

xyz = cat(2,cl.XYZ);

wh = sub2ind(volInfo.dim(1:3),xyz(1,:)',xyz(2,:)',xyz(3,:)');

imgvec(wh) = 1;

maskvec = imgvec(volInfo.wh_inmask);

return
