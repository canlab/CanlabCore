function [ivecobj, orig_cluster_indx]  = region2imagevec2tmp(cl)
% Convert a region object to an image_vector object, replacing the voxels
% and reconstructing as much info as possible.
%
% The .dat field of the new "ivecobj" is made from the cl.all_data field.
% if this is empty, uses cl.val field, then cl.Z as a backup.
% Mask information is available in ivecobj.volInfo.
%
% :Usage:
% ::
%
%    ivecobj = region2imagevec(cl)
%
% ..
%    NEEDS SOME ADDITIONAL WORK/CHECKING
% ..

ivecobj = image_vector;
ivecobj.volInfo.mat = cl(1).M;
ivecobj.volInfo.dim = cl(1).dim;

% no, will be reordered
ivecobj.volInfo.xyzlist = cat(2, cl.XYZ)';

mask = clusters2mask2011(cl, cl(1).dim);

n = sum(mask(:) ~= 0);

% tor changed april 28 2011 to be all voxels
ivecobj.removed_voxels = mask(:) == 0 | isnan(mask(:)); %false(n, 1);

ivecobj.volInfo.image_indx = mask(:) ~= 0;
ivecobj.volInfo.n_inmask = n;

ivecobj.volInfo.wh_inmask = find(ivecobj.volInfo.image_indx);

% re-get continguity; don't just assume
orig_cluster_indx = mask(:);
orig_cluster_indx = orig_cluster_indx(ivecobj.volInfo.wh_inmask);

%ivecobj.volInfo.cluster = mask(:);
%ivecobj.volInfo.cluster = ivecobj.volInfo.cluster(ivecobj.volInfo.wh_inmask);

ivecobj.volInfo.nvox = prod(cl(1).dim);

[i, j, k] = ind2sub(cl(1).dim, ivecobj.volInfo.wh_inmask);
ivecobj.volInfo.xyzlist = [i j k];

if ivecobj.volInfo.n_inmask < 50000
    ivecobj.volInfo.cluster = spm_clusters(ivecobj.volInfo.xyzlist')';
else
    ivecobj.volInfo.cluster = ones(ivecobj.volInfo.n_inmask, 1);
end

ivecobj.volInfo.dt = [16 0]; % for reslicing compatibility
ivecobj.volInfo.fname = 'Reconstructed from clusters';

ivecobj.dat = ones(length(ivecobj.volInfo.wh_inmask), 1);

% region2imagevec creates illegal list of removed_voxels (doesn't match
% .volInfo.wh_inmask).  Fix...
ivecobj.removed_voxels = ivecobj.removed_voxels(ivecobj.volInfo.wh_inmask);




end
