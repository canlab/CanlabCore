function dat = flip(dat, varargin)
% Flips images stored in an object left-to-right
%
% Operates on a single-image (3-D) object only. The slice loop below
% iterates over Z planes of a 3-D volume; calling it on a multi-image
% (4-D after reconstruction) object would silently scramble data, so
% this method now errors instead. To flip a multi-image object, reduce
% it to a single image first (e.g. mean(obj), or get_wh_image(obj, k)).
%
% :Optional Inputs:
%
%   Input 'mirror' to make a symmetrical image, averaging the left
%   and right hemispheres
%
% :Examples:
% ::
%
%    dat = flip(dat)
%    dat = flip(dat, 'mirror')
%    dat = flip(mean(obj))         % flip the group-mean of a multi-image object
%
% ..
%    Tor. may 2012
% ..

% Ensure padded state so reconstruct_image / rebuild_volinfo_from_dat see
% the full voxel space; otherwise the latter errors with
% "newdat MUST be size of ENTIRE image".
dat = replace_empty(dat);

% Single-image guard: the slice loop below assumes a 3-D reconstruction.
% Without this, a multi-image object reconstructs to 4-D, the loop
% iterates only the first image's Z planes, and rebuild_volinfo_from_dat
% errors on a vector of the wrong length.
n_images = size(dat.dat, 2);
if n_images > 1
    error('image_vector:flip:multiImage', ...
        ['flip() operates on a single-image (3-D) object, but this ' ...
         'object has %d images. Reduce to one image first, e.g. ' ...
         'mean(obj) or get_wh_image(obj, k).'], n_images);
end

vdat = reconstruct_image(dat);
orig_dat = vdat(:);

for i = 1:size(vdat, 3)
    slice = vdat(:, :, i);
    slice = flipdim(slice, 1);
    vdat(:, :, i) = slice;
end

vdat = vdat(:);

 % mirror
if ~isempty(varargin) && strcmp(varargin{1}, 'mirror')
    vdat = (orig_dat + vdat) / 2;
end

% rebuild mask around new non-zero values
dat = rebuild_volinfo_from_dat(dat, vdat);

end
