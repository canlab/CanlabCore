function fig = controller(obj, varargin)
% Open (or refresh) a uifigure controller bound to this fmridisplay instance.
%
% The controller is the MATLAB-side analog of the NiiVue control panel
% (VISUALIZATION_OVERHAUL_NOTES.md, section 6). It lists the object's blob
% layers and gives each a live opacity slider, colormap dropdown (with a
% graphical solid-colour picker), threshold field, and visibility toggle, all
% acting across montages AND surfaces (fmridisplay is a handle class). The panel:
%   - reflects the current state of each layer (colormap, opacity, threshold),
%     and updates in place when you change a layer from the command line
%     (rethreshold / set_colormap / set_opacity);
%   - is rebuilt only when layers are added or removed (addblobs / removeblobs);
%   - is type-aware for thresholding (statistic_image -> p-value, fmri_data/mask
%     -> raw value);
%   - echoes the equivalent line of code for each action to the command window,
%     e.g. `han = set_colormap(han, 'color', [1 0 0]);`;
%   - is light green, names the bound variable in its title, and closes with
%     `close all` (HandleVisibility 'on').
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
%   **fig:** the uifigure handle of the controller.
%
% :See also:
%   - addblobs, refresh, rethreshold, set_colormap, set_opacity, update_controller
%
% ..
%    2026 visualization overhaul
% ..

vname        = inputname(1);     % caller's variable name, for echoed code + title
nlayers      = numel(obj.activation_maps);
cmap_options = colormap_options();

% Reuse an already-open controller; update in place if the layer count is
% unchanged (safe to call from inside a control's own callback), else rebuild.
if ~isempty(obj.controller_handle) && isgraphics(obj.controller_handle) && isvalid(obj.controller_handle)
    fig   = obj.controller_handle;
    vname = getappdata(fig, 'objname');     % keep the original bound name across rebuilds
    existing = findobj(fig, 'Type', 'uipanel');
    if nlayers > 0 && numel(existing) == nlayers
        update_controls_in_place(fig, obj, cmap_options);
        return
    end
    delete(allchild(fig));                  % layer count changed -> rebuild
else
    fig = uifigure('Name', 'CANlab display controller', 'HandleVisibility', 'on', ...
        'Color', [0.90 0.98 0.90]);         % light green, distinct among figures
    obj.controller_handle = fig;
    if isempty(vname), vname = 'obj'; end    % fallback for echoed code / title
end

setappdata(fig, 'objname', vname);
setappdata(fig, 'fmridisplay_obj', obj);
fig.Name = sprintf('CANlab display controller  [ %s ]', vname);   % which object this controls
fig.Position(3:4) = [760, 70 + 120 * max(nlayers, 1)];            % ~2x wider, compact rows

if nlayers == 0
    uilabel(fig, 'Position', [20 fig.Position(4)-60 700 40], 'FontSize', 13, ...
        'Text', 'No blob layers yet. Add some with addblobs(obj, ...).');
    return
end

outer = uigridlayout(fig, [nlayers + 2, 1]);
outer.RowHeight = [repmat({112}, 1, nlayers), {30}, {'1x'}];      % panels, short button, spacer
outer.Scrollable = 'on';

for k = 1:nlayers
    build_layer_panel(outer, obj, k, cmap_options, vname);
end

btn = uibutton(outer, 'Text', 'Re-render layers', 'FontSize', 12, ...
    'ButtonPushedFcn', @(~, ~) refresh(obj));
btn.Layout.Row = nlayers + 1;

end


% ------------------------------------------------------------------------

function opts = colormap_options()
% Dropdown choices. The max/min-colour ramps render reliably on montages and
% surfaces; 'solid colour…' opens the built-in MATLAB colour picker. (True
% perceptual maps like inferno/viridis are deferred — see VISUALIZATION_OVERHAUL_NOTES.md.)
opts = {'split (hot/cool)', 'split (mango)', 'seafire', 'warm (red-yellow)', ...
        'cool (blue-cyan)', 'winter (blue-green)', 'solid colour…'};
end


function build_layer_panel(parent, obj, k, cmap_options, vname)
% One control panel for layer k, controls initialized from its state and tagged
% so they can be updated in place.

layer = obj.activation_maps{k};
args  = render_args_of(layer);
src   = source_of(layer);

panel = uipanel(parent, 'Title', sprintf('Layer %d  (%s)', k, source_kind(src)), 'FontSize', 13);
panel.Layout.Row = k;

% Two compact rows (saves vertical space):
%   Opacity  [==== fixed-width slider ====]        Visible [x]
%   p (unc)  [== slider with p-value ticks ==]     Colors  [dropdown]
g = uigridlayout(panel, [2, 4]);
g.RowHeight     = {'1x', '1x'};
g.ColumnWidth   = {86, 360, 64, 150};   % label | fixed wide slider | label | fixed control
g.ColumnSpacing = 8;
g.Padding       = [8 2 8 2];

fs = 12;

% Row 1: Opacity (fixed-width slider, no ticks) + Visible toggle
l1 = uilabel(g, 'Text', 'Opacity', 'FontSize', fs);   l1.Layout.Row = 1; l1.Layout.Column = 1;
sld = uislider(g, 'Limits', [0 1], 'Value', current_opacity(args), ...
    'MajorTicks', [], 'MinorTicks', [], 'FontSize', fs, 'Tag', sprintf('opacity_%d', k), ...
    'ValueChangedFcn', @(s, ~) on_opacity(obj, k, s.Value, vname));
sld.Layout.Row = 1; sld.Layout.Column = 2;
l2 = uilabel(g, 'Text', 'Visible', 'FontSize', fs);   l2.Layout.Row = 1; l2.Layout.Column = 3;
cb = uicheckbox(g, 'Text', '', 'Value', true, 'Tag', sprintf('visible_%d', k), ...
    'ValueChangedFcn', @(c, ~) set_layer_visible(obj, k, c.Value));
cb.Layout.Row = 1; cb.Layout.Column = 4;

% Row 2: threshold slider (type-aware ticks) + Colors dropdown
[thr_label, thr_value, ~, is_pval] = threshold_spec(src, layer);
l3 = uilabel(g, 'Text', thr_label, 'FontSize', fs);   l3.Layout.Row = 2; l3.Layout.Column = 1;
ts = build_threshold_slider(g, obj, k, src, is_pval, thr_value, vname, fs);
ts.Layout.Row = 2; ts.Layout.Column = 2;
l4 = uilabel(g, 'Text', 'Colors', 'FontSize', fs);    l4.Layout.Row = 2; l4.Layout.Column = 3;
dd = uidropdown(g, 'Items', cmap_options, 'Value', current_colormap_label(args, cmap_options), ...
    'FontSize', fs, 'Tag', sprintf('colormap_%d', k), ...
    'ValueChangedFcn', @(d, ~) on_colormap(obj, k, d.Value, vname));
dd.Layout.Row = 2; dd.Layout.Column = 4;

end


function sld = build_threshold_slider(g, obj, k, src, is_pval, v0, vname, fs)
% Type-aware threshold slider. statistic_image -> p-value slider with ticks at
% common thresholds; fmri_data/region -> raw |x| slider anchored at 0 and the
% 99.9th percentile of |data|.
if is_pval
    lims = [0 0.1]; ticks = [.001 .005 .01 .05 .1]; labs = {'.001','.005','.01','.05','.1'};
elseif isa(src, 'image_vector') || isa(src, 'region')
    b = raw_bound(src);
    lims = [0 b]; ticks = [0 b]; labs = {'0', sprintf('%.3g', b)};
else
    sld = uislider(g, 'Limits', [0 1], 'Value', 0, 'Enable', 'off', 'FontSize', fs, ...
        'Tag', sprintf('threshold_%d', k));
    return
end
v = min(max(v0, lims(1)), lims(2));
sld = uislider(g, 'Limits', lims, 'Value', v, 'MajorTicks', ticks, 'MajorTickLabels', labs, ...
    'FontSize', max(8, fs - 3), 'Tag', sprintf('threshold_%d', k), ...
    'ValueChangedFcn', @(s, ~) on_threshold(obj, k, s.Value, is_pval, vname));
end


function b = raw_bound(src)
% 99.9th percentile of |data| for a raw-valued layer's threshold-slider range.
b = 1;
try
    if isa(src, 'region'), iv = region2imagevec(src); else, iv = src; end
    d = double(iv.dat(:));
    d = d(d ~= 0 & ~isnan(d) & ~isinf(d));
    if ~isempty(d), b = prctile(abs(d), 99.9); end
catch
end
if ~(b > 0), b = 1; end
end


function update_controls_in_place(fig, obj, cmap_options)
% Re-sync each control's value to its layer's current state (no rebuild).
for k = 1:numel(obj.activation_maps)
    layer = obj.activation_maps{k};
    args  = render_args_of(layer);
    src   = source_of(layer);

    sld = findobj(fig, 'Tag', sprintf('opacity_%d', k));
    if ~isempty(sld), sld.Value = current_opacity(args); end

    dd = findobj(fig, 'Tag', sprintf('colormap_%d', k));
    if ~isempty(dd), dd.Value = current_colormap_label(args, cmap_options); end

    thr = findobj(fig, 'Tag', sprintf('threshold_%d', k));
    if ~isempty(thr) && isprop(thr, 'Limits')
        [~, val] = threshold_spec(src, layer);
        thr.Value = min(max(val, thr.Limits(1)), thr.Limits(2));   % clamp into slider range
    end
end
end


% ---- callbacks (apply + echo the equivalent code line) -----------------

function on_opacity(obj, k, value, vname)
set_opacity(obj, value, 'layers', k);
echo_code(vname, sprintf('set_opacity(%s, %g, ''layers'', %d)', vname, value, k));
end

function on_colormap(obj, k, choice, vname)
switch choice
    case 'split (hot/cool)'
        set_colormap(obj, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''splitcolor'', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, ''layers'', %d)', vname, k));
    case 'split (mango)'
        set_colormap(obj, 'splitcolor', {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]}, 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''splitcolor'', {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]}, ''layers'', %d)', vname, k));
    case 'seafire'
        set_colormap(obj, 'splitcolor', {[0 0 1] [.2 .8 .6] [1 .4 .1] [1 1 0]}, 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''splitcolor'', {[0 0 1] [.2 .8 .6] [1 .4 .1] [1 1 0]}, ''layers'', %d)', vname, k));
    case 'warm (red-yellow)'
        set_colormap(obj, 'maxcolor', [1 1 0], 'mincolor', [1 0 0], 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''maxcolor'', [1 1 0], ''mincolor'', [1 0 0], ''layers'', %d)', vname, k));
    case 'cool (blue-cyan)'
        set_colormap(obj, 'maxcolor', [0 1 1], 'mincolor', [0 0 1], 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''maxcolor'', [0 1 1], ''mincolor'', [0 0 1], ''layers'', %d)', vname, k));
    case 'winter (blue-green)'
        set_colormap(obj, 'maxcolor', [0 1 .5], 'mincolor', [0 0 1], 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''maxcolor'', [0 1 0.5], ''mincolor'', [0 0 1], ''layers'', %d)', vname, k));
    case 'solid colour…'
        c = uisetcolor([1 0 0], 'Choose a blob colour');
        if numel(c) ~= 3, return, end          % cancelled (uisetcolor returns 0)
        set_colormap(obj, 'color', c, 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''color'', [%g %g %g], ''layers'', %d)', vname, c(1), c(2), c(3), k));
end
end

function on_threshold(obj, k, value, is_pval, vname)
if is_pval
    rethreshold(obj, value, 'unc', 'layers', k);
    echo_code(vname, sprintf('rethreshold(%s, %g, ''unc'', ''layers'', %d)', vname, value, k));
else
    rethreshold(obj, value, 'layers', k);
    echo_code(vname, sprintf('rethreshold(%s, %g, ''layers'', %d)', vname, value, k));
end
end

function echo_code(vname, callstr)
% Print the runnable line of code that the GUI action just executed.
fprintf('%s = %s;\n', vname, callstr);
end


% ---- small state-readers -----------------------------------------------

function args = render_args_of(layer)
args = {};
if isfield(layer, 'render_args') && ~isempty(layer.render_args), args = layer.render_args; end
end

function src = source_of(layer)
src = [];
if isfield(layer, 'source_object'), src = layer.source_object; end
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

function lbl = current_colormap_label(args, opts) %#ok<INUSD>
% Map the layer's render_args back to a dropdown label (always a valid item).
if any(strcmp(args, 'splitcolor'))
    sc = args{find(strcmp(args, 'splitcolor'), 1) + 1};
    if iscell(sc) && numel(sc) == 4 && isequal(sc, {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]})
        lbl = 'split (mango)';
    elseif iscell(sc) && numel(sc) == 4 && isequal(sc, {[0 0 1] [.2 .8 .6] [1 .4 .1] [1 1 0]})
        lbl = 'seafire';
    else
        lbl = 'split (hot/cool)';
    end
elseif any(strcmp(args, 'color'))
    lbl = 'solid colour…';
elseif any(strcmp(args, 'maxcolor'))
    mx = args{find(strcmp(args, 'maxcolor'), 1) + 1};
    if     isequal(mx, [1 1 0]), lbl = 'warm (red-yellow)';
    elseif isequal(mx, [0 1 1]), lbl = 'cool (blue-cyan)';
    elseif isequal(mx, [0 1 .5]), lbl = 'winter (blue-green)';
    else,  lbl = 'warm (red-yellow)';
    end
else
    lbl = 'split (hot/cool)';
end
end

function [lbl, val, enable, is_pval] = threshold_spec(src, layer)
applied = [];
if isfield(layer, 'applied_threshold'), applied = layer.applied_threshold; end
if isa(src, 'statistic_image')
    lbl = 'p (unc)'; is_pval = true; enable = true;
    val = 0.005; if ~isempty(applied) && isscalar(applied), val = applied; end
elseif isa(src, 'image_vector') || isa(src, 'region')
    lbl = '|x| >'; is_pval = false; enable = true;
    val = 0; if ~isempty(applied) && isscalar(applied), val = applied; end
else
    lbl = 'thresh'; is_pval = false; enable = false; val = 0;
end
end

function set_layer_visible(obj, k, tf)
% Toggle the visibility of a layer's montage blob graphics. (Surface coloring is
% not yet per-layer; see VISUALIZATION_OVERHAUL_NOTES.md.)
if k < 1 || k > numel(obj.activation_maps), return, end
bh = obj.activation_maps{k}.blobhandles;
if isempty(bh), return, end
bh = bh(ishandle(bh));
set(bh, 'Visible', matlab.lang.OnOffSwitchState(tf));
end
