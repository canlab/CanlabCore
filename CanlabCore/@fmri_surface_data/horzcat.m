function c = horzcat(varargin)
% horzcat Concatenate fmri_surface_data objects with [a, b, ...] syntax.
%
% :Usage:
% ::
%     combined = [obj1, obj2, obj3];
%
% Thin wrapper over cat; concatenates along the image (map) dimension. All
% objects must share the same grayordinate space.
%
% :See also: cat, fmri_surface_data

if nargin == 1
    c = varargin{1};
    return
end
c = cat(varargin{1}, varargin{2:end});
end
