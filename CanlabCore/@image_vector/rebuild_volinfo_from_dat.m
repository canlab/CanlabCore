function dat = rebuild_volinfo_from_dat(dat, newdat)
% Will rebuild volInfo (the image space, or sometimes "mask") from a vectorized image. 
% In other words, will rebuild dat.volInfo from newdat.
%
% Also resets all voxels to be significant, if a statistic image
%
% INPUT
%  dat:  an image_vector
%  newdat: a vector that MUST be size of ENTIRE image (dat.volInfo.nvox)
%
% OUTPUT
%  dat:  dat.dat contains the non-zero values of newdat, and dat.volInfo is
%  correctly defining the image space

if length(newdat) ~= dat.volInfo.nvox
    error('newdat MUST be size of ENTIRE image (dat.volInfo.nvox).  Image_vector.reconstruct_image may be helpful');
end

% rebuild volInfo from newdat
dat.volInfo.image_indx = newdat ~=0 & ~isnan(newdat);
dat.volInfo.wh_inmask = find(dat.volInfo.image_indx);
dat.volInfo.n_inmask = length(dat.volInfo.wh_inmask);
[i, j, k] = ind2sub(dat.volInfo(1).dim(1:3), dat.volInfo(1).wh_inmask);
dat.volInfo(1).xyzlist = [i j k];
        
% rebuild dat.dat with new data
dat.dat = newdat(dat.volInfo.image_indx);

% rebuild volInfo.cluster
dat.volInfo.cluster = dat.volInfo.wh_inmask; % this forces reparse_contiguous to rebuild the cluster correctly.  maybe not ideal solution, but maybe yes -- don't know the code well enough.  Yoni 7/14
dat = reparse_contiguous(dat);  

% rebuild sig field if needed
if isprop(dat, 'sig'), dat.sig = true(size(dat.dat)); end
