function dat = flip(dat, varargin)
% Flips an image_vector object left to right
%
% Optional: input 'mirror' to make a symmetrical image, averaging the left
% and right hemispheres
%
% dat = flip(dat, ['mirror'])
%
% tor. may 2012


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

