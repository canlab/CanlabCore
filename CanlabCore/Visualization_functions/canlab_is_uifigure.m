function tf = canlab_is_uifigure(f)
% Return true if F is a uifigure (App Designer / uifigure) rather than a
% traditional figure.
%
% :Usage:
% ::
%
%     tf = canlab_is_uifigure(f)          % f defaults to gcf if omitted
%
% Detected via AutoResizeChildren='on', which is the uifigure default (a
% traditional figure is 'off') AND is exactly the property that makes
% subplot()/axes() fail with "Adding subplots to a container with the
% 'AutoResizeChildren' property set to 'on' is not supported". This is used by
% the CANlab display methods (fmridisplay.montage/surface, image_vector.montage/
% surface, canlab_results_fmridisplay) to avoid rendering slices/surfaces into
% the interactive controller window: they open a fresh traditional figure
% instead, so figure-creation behaviour stays consistent across all of them.
%
% NOTE: matlab.ui.internal.isUIFigure returns true for BOTH traditional and
% uifigures, so it cannot be used to tell them apart.
%
% :See also: fmridisplay, canlab_results_fmridisplay

if nargin < 1, f = get(groot, 'CurrentFigure'); end

tf = false;
if isempty(f) || ~isgraphics(f), return, end
try
    tf = isprop(f, 'AutoResizeChildren') && strcmp(f.AutoResizeChildren, 'on');
catch
end
end
