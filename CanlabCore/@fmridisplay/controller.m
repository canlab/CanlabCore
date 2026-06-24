function fig = controller(obj, varargin)
% Open a uifigure controller bound to this fmridisplay instance.
%
% The controller is the MATLAB-side analog of the NiiVue control panel
% (VISUALIZATION_OVERHAUL_NOTES.md, section 6). It lists the object's blob
% layers and gives each a live opacity slider, colormap dropdown, threshold
% control, and a visibility toggle. Because fmridisplay is a handle class,
% every control mutates the SAME object the montages were drawn from and
% calls the corresponding in-place method (set_opacity / set_colormap /
% rethreshold / refresh), so the montage figure updates immediately.
%
% :Usage:
% ::
%
%     o2  = canlab_results_fmridisplay(t);
%     fig = controller(o2);          % open the control panel for o2
%
% :Inputs:
%
%   **obj:**
%        An fmridisplay object (handle) with one or more blob layers
%        (activation_maps). Open a montage and add blobs first.
%
% :Outputs:
%
%   **fig:**
%        The uifigure handle of the controller. The owning fmridisplay
%        instance is stored in its appdata under 'fmridisplay_obj' (a
%        reference, not a copy), so callbacks route back to the one object.
%
% :Examples:
% ::
%
%     imgs = load_image_set('emotionreg');
%     t    = threshold(ttest(imgs), .005, 'unc');
%     o2   = canlab_results_fmridisplay(t);
%     controller(o2);   % drag opacity / change colors / re-threshold live
%
% :See also:
%   - addblobs, refresh, rethreshold, set_colormap, set_opacity
%
% ..
%    2026 visualization overhaul
% ..

nlayers = numel(obj.activation_maps);
if nlayers == 0
    warning('fmridisplay:controller:noLayers', ...
        'No blob layers to control. Add blobs first (addblobs / canlab_results_fmridisplay).');
end

% ---- Figure + scrolling grid -------------------------------------------
fig = uifigure('Name', 'CANlab display controller', 'Position', [80 80 360 120 + 150 * max(nlayers, 1)]);
setappdata(fig, 'fmridisplay_obj', obj);   % back-pointer (reference) to the instance

outer = uigridlayout(fig, [nlayers + 1, 1]);
outer.RowHeight = [repmat({150}, 1, nlayers), {'1x'}];
outer.Scrollable = 'on';

cmap_options = {'split (hot/cool)', 'solid red', 'solid green', 'solid blue', 'orange-yellow'};

for k = 1:nlayers
    build_layer_panel(outer, obj, k, cmap_options);
end

% ---- Footer: refresh-all button ----------------------------------------
btn = uibutton(outer, 'Text', 'Re-render all layers', ...
    'ButtonPushedFcn', @(~, ~) refresh(obj));
btn.Layout.Row = nlayers + 1;

end


function build_layer_panel(parent, obj, k, cmap_options)
% One control panel (opacity / colormap / threshold / visibility) for layer k.

panel = uipanel(parent, 'Title', sprintf('Layer %d', k));
panel.Layout.Row = k;

g = uigridlayout(panel, [4, 2]);
g.RowHeight   = {22, 22, 22, 22};
g.ColumnWidth = {110, '1x'};

% --- Opacity slider ---
uilabel(g, 'Text', 'Opacity');
sld = uislider(g, 'Limits', [0 1], 'Value', 1, ...
    'ValueChangedFcn', @(s, ~) set_opacity(obj, s.Value, 'layers', k));
sld.Layout.Column = 2;

% --- Colormap dropdown ---
uilabel(g, 'Text', 'Colors');
dd = uidropdown(g, 'Items', cmap_options, 'Value', cmap_options{1}, ...
    'ValueChangedFcn', @(d, ~) apply_colormap_choice(obj, k, d.Value));
dd.Layout.Column = 2;

% --- Threshold control (only meaningful for statistic_image / image_vector sources) ---
src = [];
if isfield(obj.activation_maps{k}, 'source_object'), src = obj.activation_maps{k}.source_object; end
can_rethresh = isa(src, 'statistic_image') || isa(src, 'image_vector');

uilabel(g, 'Text', 'p-threshold (unc)');
thr = uieditfield(g, 'numeric', 'Value', 0.005, 'Limits', [0 1], ...
    'Enable', matlab.lang.OnOffSwitchState(can_rethresh), ...
    'ValueChangedFcn', @(e, ~) rethreshold(obj, e.Value, 'unc', 'layers', k));
thr.Layout.Column = 2;

% --- Visibility toggle ---
uilabel(g, 'Text', 'Visible');
cb = uicheckbox(g, 'Text', '', 'Value', true, ...
    'ValueChangedFcn', @(c, ~) set_layer_visible(obj, k, c.Value));
cb.Layout.Column = 2;

end


function apply_colormap_choice(obj, k, choice)
% Map a friendly dropdown label to set_colormap options.
switch choice
    case 'split (hot/cool)'
        set_colormap(obj, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, 'layers', k);
    case 'solid red'
        set_colormap(obj, 'color', [1 0 0], 'layers', k);
    case 'solid green'
        set_colormap(obj, 'color', [0 1 0], 'layers', k);
    case 'solid blue'
        set_colormap(obj, 'color', [0 0 1], 'layers', k);
    case 'orange-yellow'
        set_colormap(obj, 'maxcolor', [1 1 0], 'mincolor', [1 .3 0], 'layers', k);
end
end


function set_layer_visible(obj, k, tf)
% Toggle the visibility of a layer's blob graphics without deleting them.
if k < 1 || k > numel(obj.activation_maps), return, end
bh = obj.activation_maps{k}.blobhandles;
if isempty(bh), return, end
bh = bh(ishandle(bh));
set(bh, 'Visible', matlab.lang.OnOffSwitchState(tf));
end
