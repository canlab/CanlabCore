function dat = rebuild_volinfo_from_dat(dat, newdat)
% Will rebuild volInfo (the image space, or sometimes "mask") from a vectorized image. 
% In other words, will rebuild dat.volInfo from newdat.
%
% Also resets all voxels to be significant, if a statistic image
%
% :Input:
%
%   **dat:**
%        an image_vector
%
%   **newdat:**
%        a vector that MUST be size of ENTIRE image (dat.volInfo.nvox)
%
% :Output:
%
%   **dat:**
%        dat.dat contains the non-zero values of newdat, and dat.volInfo is
%        correctly defining the image space


if length(newdat) ~= dat.volInfo.nvox
    error('newdat MUST be size of ENTIRE image (dat.volInfo.nvox).  Image_vector.reconstruct_image may be helpful');
end

wh_new = newdat ~=0 & ~isnan(newdat);

if isa(dat, 'statistic_image')
    % Rebuild fields specific to statistic_images
    
    dat = replace_empty(dat);
    k = size(dat.dat, 2);
    
    if isempty(dat.ste), dat.ste = single(zeros(size(dat.dat))); end
    if isempty(dat.sig), dat.sig = single(zeros(size(dat.dat))); end
    
    for i = 1:k
        
        p = ones(dat.volInfo.nvox, k);
        p(dat.volInfo.wh_inmask, i) = dat.p(:, i);
        
        ste = Inf .* ones(dat.volInfo.nvox, k);
        ste(dat.volInfo.wh_inmask, i) = dat.ste(:, i);
        
        sig = zeros(dat.volInfo.nvox, k);
        sig(dat.volInfo.wh_inmask, i) = dat.sig(:, i);
        
    end
    
    dat.p = p(wh_new, :);
    dat.ste = ste(wh_new, :);
    dat.sig = sig(wh_new, :);
    
end

% rebuild volInfo from newdat
dat.volInfo.image_indx = wh_new;
dat.volInfo.wh_inmask = find(dat.volInfo.image_indx);
dat.volInfo.n_inmask = length(dat.volInfo.wh_inmask);
[i, j, k] = ind2sub(dat.volInfo(1).dim(1:3), dat.volInfo(1).wh_inmask);
dat.volInfo(1).xyzlist = [i j k];
        
% rebuild dat.dat with new data
dat.dat = newdat(dat.volInfo.image_indx);

dat.removed_voxels = false(dat.volInfo.n_inmask, 1);

% rebuild volInfo.cluster
dat.volInfo.cluster = dat.volInfo.wh_inmask; 
% this forces reparse_contiguous to rebuild the cluster correctly.  
% maybe not ideal solution, but maybe yes -- don't know the code well enough.  Yoni 7/14

dat = reparse_contiguous(dat);  


end
