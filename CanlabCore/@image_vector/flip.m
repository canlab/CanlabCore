function dat = flip(dat, varargin)
% Flips images stored in an object left-to-right
%
% :Optional Inputs:
%
%   Input 'mirror' to make a symmetrical image, averaging the left
%   and right hemispheres
%
% :Examples:
% ::
%
%    dat = flip(dat, ['mirror'])
%
% ..
%    Tor. may 2012
% ..

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

