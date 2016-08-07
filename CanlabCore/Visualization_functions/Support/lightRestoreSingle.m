function [h,az,el] = lightRestoreSingle(axishandle)
% [h,az,el] = lightRestoreSingle(axishandle)
% delete all lights from axis and install only a single one
% return handle in h

if nargin == 0, axishandle = gca; end

lh = lighthandles(axishandle);
delete(lh)

[az,el]=view; h = lightangle(az,el);

return
