function obj = update_controller(obj)
% Rebuild an open controller panel so it reflects the current layers.
%
% Called by addblobs / removeblobs after the layer set changes, so a controller
% that is already open redraws itself (new layers appear, removed layers vanish,
% and each control reflects current state). No-op when no controller is open.
%
% :See also:
%   - controller, addblobs, removeblobs
%
% ..
%    2026 visualization overhaul
% ..

if ~isempty(obj.controller_handle) && isgraphics(obj.controller_handle) && isvalid(obj.controller_handle)
    controller(obj);   % rebuilds contents in place
end

end
