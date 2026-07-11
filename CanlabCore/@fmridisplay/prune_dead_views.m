function obj = prune_dead_views(obj, doverbose)
% Drop montage/surface views whose figures the user has closed.
%
% When a montage or surface figure window is closed, its axes/patch handles
% become invalid. Any later addblobs / removeblobs / refresh that tries to draw
% on or erase those handles would error ("Invalid or deleted object"). This
% method detects views with no live graphics and removes them from the object,
% with a short note, so the display keeps working after windows are closed.
%
% :Usage:
% ::
%
%     obj = prune_dead_views(obj)            % verbose note when pruning
%     obj = prune_dead_views(obj, false)     % silent
%
% Called automatically at the top of addblobs, removeblobs, and refresh.
%
% :See also:
%   - addblobs, removeblobs, refresh, render_layer_surfaces
%
% ..
%    2026 visualization overhaul
% ..

if nargin < 2, doverbose = true; end

% --- Surfaces: keep those with at least one live patch handle ---
if ~isempty(obj.surface)
    keep = true(1, numel(obj.surface));
    for i = 1:numel(obj.surface)
        h = obj.surface{i}.object_handle;
        keep(i) = ~isempty(h) && any(ishandle(h));
    end
    if any(~keep) && doverbose
        fprintf(['Note: %d surface figure(s) were closed; removing those ' ...
            'view(s) from the fmridisplay object.\n'], sum(~keep));
    end
    obj.surface = obj.surface(keep);
end

% --- Montages: keep those with at least one live axis handle ---
if ~isempty(obj.montage)
    keep = true(1, numel(obj.montage));
    for i = 1:numel(obj.montage)
        ah = obj.montage{i}.axis_handles;
        keep(i) = ~isempty(ah) && any(ishandle(ah));
    end
    if any(~keep) && doverbose
        fprintf(['Note: %d montage figure(s) were closed; removing those ' ...
            'view(s) from the fmridisplay object.\n'], sum(~keep));
    end
    obj.montage = obj.montage(keep);
end

end
