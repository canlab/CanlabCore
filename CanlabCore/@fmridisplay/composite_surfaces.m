function obj = composite_surfaces(obj, wh_surface)
% Re-render all blob layers onto the surface(s) as a true-colour RGB composite.
%
% Resets each surface to its saved anatomy gray, then paints every blob layer in
% order (bottom to top), compositing in RGB so the top layer wins per vertex and
% lower layers show through where the top layer has no colour. Used by refresh
% and remove_layer (a full recompute); addblobs composites a single new layer on
% top incrementally via render_layer_surfaces.
%
% :Usage:
% ::
%
%     obj = composite_surfaces(obj)              % recomposite all surfaces
%     obj = composite_surfaces(obj, wh_surface)  % only the listed surface indices
%
% :See also:
%   - render_layer_surfaces, refresh, removeblobs, remove_layer, render_on_surface
%
% ..
%    2026 visualization overhaul — central colour pipeline
% ..

if isempty(obj.surface), return, end
if nargin < 2 || isempty(wh_surface), wh_surface = 1:numel(obj.surface); end

% Reset each target surface to its saved anatomy gray
for s = wh_surface
    if s < 1 || s > numel(obj.surface), continue, end
    h = obj.surface{s}.object_handle;
    h = h(ishandle(h));
    if isempty(h), continue, end
    if ~isempty(get(h(1), 'UserData'))
        addbrain('eraseblobs', h);
    end
end

% Paint every layer in order; each composites onto the running result
for k = 1:numel(obj.activation_maps)
    obj = render_layer_surfaces(obj, k, wh_surface);
end

end
