function [mode, info] = canlab_color_mode(dat, varargin)
% canlab_color_mode Single source of truth for choosing a blob color mode.
%
% Decides whether an image/blob set should be rendered with a continuous
% COLORMAP, one solid color per region ('unique'), or a single SOLID color,
% based on how many distinct values it has. Used by montage / surface / orthviews
% and the fmridisplay controller so the default is uniform across methods.
%
% :Usage:
% ::
%     mode = canlab_color_mode(dat)
%     mode = canlab_color_mode(dat, 'max_unique', 300)
%
% :Inputs:
%   **dat:** a numeric data vector/matrix, an image_vector/atlas/region object,
%            or an fmri_surface_data. Zeros and NaNs are treated as background.
%
% :Optional Inputs:
%   **'max_unique':** the cutoff (default 1000) for a *generic* integer map.
%                     Integer-valued images with more than 2 and fewer than this
%                     many distinct nonzero values use 'unique'; at or above it (or
%                     non-integer data) use 'colormap'. Atlas objects and
%                     .dlabel/.label parcellations are always 'unique', ignoring
%                     this cutoff. (image_vector.orthviews uses a 300 parcel cutoff
%                     for its own legacy auto-detection.)
%
% :Outputs:
%   **mode:** one of
%       'solid'    - <= 2 distinct nonzero values (a binary mask / constant map):
%                    render in one solid color.
%       'unique'   - integer-valued with 3..(max_unique-1) distinct values, or an
%                    atlas: one distinct solid color per region.
%       'colormap' - many values or non-integer (continuous) data: a graded colormap.
%   **info:** struct with .n_unique, .is_integer, .max_unique, .is_atlas.
%
% :Examples:
% ::
%     canlab_color_mode(load_atlas('julich_fmriprep20'))     % -> 'unique'
%     canlab_color_mode(mask_obj)                            % -> 'solid'  (binary)
%     canlab_color_mode(ttest(load_image_set('emotionreg'))) % -> 'colormap'
%
% :See also: canlab_colormap, render_blobs, fmri_surface_data.surface, montage

max_unique = 1000;                         % shared cutoff for a *generic* integer map
wh = find(strcmpi(varargin, 'max_unique'), 1);
if ~isempty(wh), max_unique = varargin{wh + 1}; end

% Atlases and label parcellations (.dlabel/.label) are region-labeled: always
% 'unique', regardless of how many regions (parcellations can exceed the cutoff).
is_atlas = isa(dat, 'atlas') || ...
    (isobject(dat) && isprop(dat, 'imagetype') && ...
     any(strcmpi(dat.imagetype, {'dlabel', 'label'})));

% Extract a numeric data vector from whatever was passed.
if isa(dat, 'atlas')
    v = double(dat.dat(:));
elseif isa(dat, 'region')
    v = cat(1, dat.val);  if isempty(v), v = (1:numel(dat))'; end
elseif isobject(dat) && isprop(dat, 'dat')
    v = double(dat.dat(:));
else
    v = double(dat(:));
end

v = v(isfinite(v) & v ~= 0);               % drop background (0) and NaN
u = unique(v);
n = numel(u);
is_integer = ~isempty(u) && all(u == round(u));

info = struct('n_unique', n, 'is_integer', is_integer, ...
    'max_unique', max_unique, 'is_atlas', is_atlas);

if is_atlas
    mode = 'unique';
elseif n <= 2
    mode = 'solid';                        % binary mask / (near-)constant map
elseif is_integer && n < max_unique
    mode = 'unique';                       % few integer labels -> one color each
else
    mode = 'colormap';                     % continuous / many values
end
end
