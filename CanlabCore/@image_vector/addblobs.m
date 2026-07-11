function varargout = addblobs(obj, varargin)
% addblobs is a method of fmridisplay, not of image data objects.
%
% This stub exists only to replace MATLAB's cryptic "Undefined function
% 'addblobs' for input arguments of type 'statistic_image'" with a message
% that tells you what to do instead. addblobs adds blobs to a *display*
% object (fmridisplay); you first create that display from your image, then
% add blobs to it.
%
% :Usage:
% ::
%
%     o2 = canlab_results_fmridisplay(img);   % or:  o2 = montage(img)
%     o2 = addblobs(o2, region(img));
%
% :See also:
%   - fmridisplay, fmridisplay.addblobs, canlab_results_fmridisplay,
%     image_vector.montage, image_vector.surface
%
% ..
%    2026 visualization overhaul
% ..

error('image_vector:addblobs:notADisplay', ...
    ['addblobs adds blobs to an fmridisplay display object, not directly to a %s.\n' ...
     'Create a display from your image first, then add blobs to it:\n\n' ...
     '   o2 = canlab_results_fmridisplay(img);   %% or:  o2 = montage(img)\n' ...
     '   o2 = addblobs(o2, region(img));\n\n' ...
     'To just threshold-and-show an image, montage(img) or orthviews(img) are usually enough.'], ...
    class(obj));

varargout = {};   % never reached (error above); present to satisfy the signature

end
