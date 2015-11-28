function dat = interpolate(dat, varargin)
% Interpolate over missing values in image_vector object
%
% :Usage:
% ::
%
%    dat = interpolate(dat, varargin)
%
% :Input: 
%    image_vector object (dat; e.g., an fmri_data object)
%
% Use when there are some missing values in the mask image
% Performs 3-D linear interpolation to fill in all values in the original
% mask.
%
% e.g., For a standard brain image space that is 91 x 109 x 91, you may
% have 300,000 in-mask values. Only 150,000 of these may be defined in the
% image, however, and the rest are missing (0 or NaN).
% This function will return a dat image with non-missing values for all
% 300,000 voxels (the "in-mask" space). 
% It will not return values for all voxels in the 91 x 109 x 91 space,
% however.
%
% :Note:
% This function does not upsample the data now, but could be extended
% to do so fairly easily.
% 

upsamplevalue = 1; % values > 1 would upsample the data

SPACE = define_sampling_space(dat.volInfo, upsamplevalue);  

% voldata = reconstruct_image(dat);
% voldata(voldata(:) == 0 | isnan(voldata(:))) = 0;

xyzfull = dat.volInfo.xyzlist;

dat = remove_empty(dat);
wh = ~dat.removed_voxels;
x = xyzfull(wh, 1);
y = xyzfull(wh, 2);
z = xyzfull(wh, 3);

% interpolation of scattered data
% note: y and x are reversed in some Matlab functions, such as this one - OK
Vq = griddata(y, x, z, double(dat.dat),  SPACE.X, SPACE.Y, SPACE.Z, 'linear');

datnew = Vq(:);
datnew = datnew(dat.volInfo.wh_inmask);

dat.dat = datnew;
dat.removed_voxels = false;

if upsamplevalue > 1
    % Need to re-create image space and dims from scratch
    % up-sample mask as well
    
end

end % function


% view
%slice(SPACE.X, SPACE.Y, SPACE.Z, Vq,[6 30 60],2,[10 40 70]), shading flat

% resampled_dat = interp3(x, y, z, dat.dat, SPACE.X, SPACE.Y, SPACE.Z);
% resampled_dat = interp3(SPACE.Xo, SPACE.Yo, SPACE.Zo, voldata, SPACE.X, SPACE.Y, SPACE.Z);

