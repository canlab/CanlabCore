function fig = controller(obj, varargin)
% Open (or refresh) a uifigure controller bound to this fmridisplay instance.
%
% The controller is the MATLAB-side analog of the NiiVue control panel
% (VISUALIZATION_OVERHAUL_NOTES.md, section 6). It lists the object's blob
% layers and gives each a live opacity slider, colormap dropdown, threshold
% control, and a visibility toggle, all acting across montages AND surfaces
% (fmridisplay is a handle class). The panel:
%   - reflects the current state of each layer (colormap, opacity, threshold);
%   - is rebuilt in place (not duplicated) on a repeat call, and is refreshed
%     automatically by addblobs / removeblobs as layers are added or removed;
%   - is type-aware for thresholding: a statistic_image layer thresholds by
%     p-value, an fmri_data / mask layer by raw value;
%   - closes with `close all` (its HandleVisibility is 'on').
%
% :Usage:
% ::
%
%     o2  = canlab_results_fmridisplay(t);
%     fig = controller(o2);          % open / refresh the control panel for o2
%
% :Inputs:
%
%   **obj:** an fmridisplay object (handle) with one or more blob layers.
%
% :Outputs:
%
%   **fig:** the uifigure handle of the controller (also stored on
%            obj.controller_handle and in the figure appdata 'fmridisplay_obj').
%
% :See also:
%   - addblobs, refresh, rethreshold, set_colormap, set_opacity, update_controller
%
% ..
%    2026 visualization overhaul
% ..

nlayers = numel(obj.activation_maps);

% Reuse an already-open controller for this instance; otherwise create one.
if ~isempty(obj.controller_handle) && isgraphics(obj.controller_handle) && isvalid(obj.controller_handle)
    fig = obj.controller_handle;
    delete(allchild(fig));      % clear to rebuild from current state
else
    fig = uifigure('Name', 'CANlab display controller', 'HandleVisibility', 'on');
    setappdata(fig, 'fmridisplay_obj', obj);   % back-pointer (reference)
    obj.controller_handle = fig;
end

fig.Position(3:4) = [380, 110 + 165 * max(nlayers, 1)];

if nlayers == 0
    uilabel(fig, 'Position', [20 fig.Position(4)-60 340 40], ...
        'Text', 'No blob layers yet. Add some with addblobs(obj, ...).');
    return
end

cmap_options = colormap_options();

outer = uigridlayout(fig, [nlayers + 1, 1]);
outer.RowHeight = [repmat({155}, 1, nlayers), {'1x'}];
outer.Scrollable = 'on';

for k = 1:nlayers
    build_layer_panel(outer, obj, k, cmap_options);
end

btn = uibutton(outer, 'Text', 'Re-render all layers', ...
    'ButtonPushedFcn', @(~, ~) refresh(obj));
btn.Layout.Row = nlayers + 1;

end


% ------------------------------------------------------------------------

function opts = colormap_options()
opts = {'split (hot/cool)', 'solid red', 'solid green', 'solid blue', 'orange-yellow'};
end


function build_layer_panel(parent, obj, k, cmap_options)
% One control panel for layer k, with controls initialized from its state.

layer = obj.activation_maps{k};
args  = {};
if isfield(layer, 'render_args') && ~isempty(layer.render_args), args = layer.render_args; end

src = [];
if isfield(layer, 'source_object'), src = layer.source_object; end

panel = uipanel(parent, 'Title', sprintf('Layer %d  (%s)', k, source_kind(src)));
panel.Layout.Row = k;

g = uigridlayout(panel, [4, 2]);
g.RowHeight   = {22, 22, 22, 22};
g.ColumnWidth = {120, '1x'};

% --- Opacity slider (initialized from stored transvalue) ---
uilabel(g, 'Text', 'Opacity');
sld = uislider(g, 'Limits', [0 1], 'Value', current_opacity(args), ...
    'ValueChangedFcn', @(s, ~) set_opacity(obj, s.Value, 'layers', k));
sld.Layout.Column = 2;

% --- Colormap dropdown (initialized to the current scheme) ---
uilabel(g, 'Text', 'Colors');
dd = uidropdown(g, 'Items', cmap_options, 'Value', current_colormap_label(args, cmap_options), ...
    'ValueChangedFcn', @(d, ~) apply_colormap_choice(obj, k, d.Value));
dd.Layout.Column = 2;

% --- Threshold control (type-aware) ---
[thr_label, thr_value, thr_enable, is_pval] = threshold_spec(src, layer);
uilabel(g, 'Text', thr_label);
thr = uieditfield(g, 'numeric', 'Value', thr_value, ...
    'Enable', matlab.lang.OnOffSwitchState(thr_enable), ...
    'ValueChangedFcn', @(e, ~) apply_threshold(obj, k, e.Value, is_pval));
thr.Layout.Column = 2;

% --- Visibility toggle ---
uilabel(g, 'Text', 'Visible');
cb = uicheckbox(g, 'Text', '', 'Value', true, ...
    'ValueChangedFcn', @(c, ~) set_layer_visible(obj, k, c.Value));
cb.Layout.Column = 2;

end


function kind = source_kind(src)
if isa(src, 'statistic_image'),  kind = 'statistic_image';
elseif isa(src, 'image_vector'), kind = class(src);
elseif isa(src, 'region'),       kind = 'region';
else,                            kind = 'blobs';
end
end


function v = current_opacity(args)
v = 1;
wh = find(strcmp(args, 'transvalue'), 1);
if ~isempty(wh) && isnumeric(args{wh + 1}), v = max(0, min(1, args{wh + 1})); end
end


function lbl = current_colormap_label(args, opts)
% Map the layer's render_args back to the closest dropdown label.
if any(strcmp(args, 'splitcolor'))
    lbl = opts{1};
elseif any(strcmp(args, 'color'))
    c = args{find(strcmp(args, 'color'), 1) + 1};
    if     isequal(c, [1 0 0]), lbl = opts{2};
    elseif isequal(c, [0 1 0]), lbl = opts{3};
    elseif isequal(c, [0 0 1]), lbl = opts{4};
    else,                       lbl = opts{2};
    end
elseif any(strcmp(args, 'maxcolor'))
    lbl = opts{5};
else
    lbl = opts{1};
end
end


function [lbl, val, enable, is_pval] = threshold_spec(src, layer)
% Decide the threshold control's label / default / enabled-state by source type.
applied = [];
if isfield(layer, 'applied_threshold'), applied = layer.applied_threshold; end

if isa(src, 'statistic_image')
    lbl = 'p-threshold (unc)'; is_pval = true; enable = true;
    val = 0.005; if ~isempty(applied) && isscalar(applied), val = applied; end
elseif isa(src, 'image_vector')
    lbl = 'value cutoff (|x|>)'; is_pval = false; enable = true;
    val = 0; if ~isempty(applied) && isscalar(applied), val = applied; end
elseif isa(src, 'region')
    lbl = 'value cutoff (|x|>)'; is_pval = false; enable = true;
    val = 0; if ~isempty(applied) && isscalar(applied), val = applied; end
else
    lbl = 'threshold (n/a)'; is_pval = false; enable = false; val = 0;
end
end


function apply_threshold(obj, k, value, is_pval)
if is_pval
    rethreshold(obj, value, 'unc', 'layers', k);
else
    rethreshold(obj, value, 'layers', k);
end
end


function apply_colormap_choice(obj, k, choice)
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
% Toggle the visibility of a layer's montage blob graphics.
if k < 1 || k > numel(obj.activation_maps), return, end
bh = obj.activation_maps{k}.blobhandles;
if isempty(bh), return, end
bh = bh(ishandle(bh));
set(bh, 'Visible', matlab.lang.OnOffSwitchState(tf));
end
