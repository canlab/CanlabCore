function fig = controller(obj, varargin)
% Open (or refresh) a uifigure controller bound to this fmridisplay instance.
%
% The controller is the MATLAB-side analog of the NiiVue control panel
% (VISUALIZATION_OVERHAUL_NOTES.md, section 6). It lists the object's blob
% layers; each layer panel has a colormap stripe, an opacity slider, a
% type-aware threshold slider, a colormap dropdown with a live preview swatch
% (and a graphical solid-colour picker), and a visibility toggle — all acting
% across montages AND surfaces (fmridisplay is a handle class). A footer groups
% global Re-render / Remove legend / Close buttons. The panel:
%   - reflects the current state of each layer and updates in place when you
%     change a layer from the command line (rethreshold/set_colormap/set_opacity);
%   - is rebuilt only when layers are added or removed (addblobs/removeblobs);
%   - is type-aware for thresholding (statistic_image -> p-value on a LOG slider,
%     fmri_data/mask -> raw |x| slider);
%   - echoes the equivalent line of code for each action to the command window;
%   - names the bound variable in its title and closes with `close all`.
%
% To change the controller window colour, edit FIG_COLOR just below.
%
% :Usage:
% ::
%
%     o2  = canlab_results_fmridisplay(t);
%     fig = controller(o2);
%
% :See also:
%   - addblobs, refresh, rethreshold, set_colormap, set_opacity, remove_legend,
%     update_controller
%
% ..
%    2026 visualization overhaul
% ..

% ===> Controller window background colour. Tweak this RGB triplet manually. <===
FIG_COLOR = [1 0.5 0];

vname        = inputname(1);     % caller's variable name, for echoed code + title
nlayers      = numel(obj.activation_maps);
cmap_options = colormap_options();

% Reuse an already-open controller; update in place if the layer count is
% unchanged (safe to call from inside a control's own callback), else rebuild.
if ~isempty(obj.controller_handle) && isgraphics(obj.controller_handle) && isvalid(obj.controller_handle)
    fig   = obj.controller_handle;
    vname = getappdata(fig, 'objname');
    is_new_fig = false;
    existing = findobj(fig, 'Type', 'uipanel');
    if nlayers > 0 && numel(existing) == nlayers
        update_controls_in_place(fig, obj, cmap_options);
        return
    end
    delete(allchild(fig));
else
    fig = uifigure('Name', 'CANlab display controller', 'HandleVisibility', 'on', ...
        'Color', FIG_COLOR);
    obj.controller_handle = fig;
    is_new_fig = true;
    if isempty(vname), vname = 'obj'; end
end

% Prefer the caller's real workspace variable name over 'obj' (for the title and
% the echoed code lines).
vn = obj_varname(obj); if ~isempty(vn), vname = vn; end

setappdata(fig, 'objname', vname);
setappdata(fig, 'fmridisplay_obj', obj);
fig.Name = sprintf('CANlab display controller  [ %s ]', vname);
fig.Position(3:4) = [460, 60 + 262 * max(nlayers, 1)];

% Position a NEW controller in the top-right of the screen so it does not overlap
% the montage / surface figures (which open toward the centre/left). Only on
% creation, so a controller the user has moved stays put on later refreshes.
if is_new_fig
    ss = get(0, 'ScreenSize');
    fig.Position(1:2) = [max(10, ss(3) - fig.Position(3) - 30), ...
                         max(60, ss(4) - fig.Position(4) - 90)];
end

if nlayers == 0
    uilabel(fig, 'Position', [20 fig.Position(4)-70 420 50], 'FontSize', 17, ...
        'Text', 'No blob layers yet. Add some with addblobs(obj, ...).');
    return
end

outer = uigridlayout(fig, [nlayers + 1, 1]);
outer.RowHeight = [repmat({250}, 1, nlayers), {40}];
outer.Scrollable = 'on';

for k = 1:nlayers
    build_layer_panel(outer, obj, k, cmap_options, vname);
end

build_footer(outer, obj, fig, nlayers + 1);

end


% ------------------------------------------------------------------------

function opts = colormap_options()
opts = {'split (hot/cool)', 'split (mango)', 'seafire', 'warm (red-yellow)', ...
        'cool (blue-cyan)', 'winter (blue-green)', ...
        'viridis', 'inferno', 'magma', 'plasma', 'turbo', 'parula', ...
        'indexed (atlas)', 'solid colour…'};
end

function names = perceptual_names()
names = {'viridis', 'inferno', 'magma', 'plasma', 'turbo', 'parula'};
end

function lbl = match_perceptual_name(cmap)
% Identify which named perceptual colormap a stored LUT matrix is (for the
% dropdown). Falls back to 'viridis' for an unrecognised custom matrix.
lbl = 'viridis';
if ~isnumeric(cmap) || size(cmap, 2) ~= 3, return, end
names = perceptual_names();
for i = 1:numel(names)
    ref = canlab_perceptual_colormap(names{i}, size(cmap, 1));
    if isequal(size(ref), size(cmap)) && max(abs(ref(:) - cmap(:))) < 1e-9
        lbl = names{i}; return
    end
end
end

function fs = base_fontsize(), fs = 16; end   % control font size (was 12; +4)


function build_layer_panel(parent, obj, k, cmap_options, vname)
% One single-column panel per layer.

layer = obj.activation_maps{k};
args  = render_args_of(layer);
src   = source_of(layer);
fs    = base_fontsize();

panel = uipanel(parent, 'Title', sprintf('Layer %d  (%s)', k, source_kind(src)), 'FontSize', fs + 1);
panel.Layout.Row = k;

g = uigridlayout(panel, [6, 2]);
% colour bar (~text height), legend value labels, opacity, threshold(+ticks), colors, visible
g.RowHeight     = {26, 22, 34, 44, 34, 30};
g.ColumnWidth   = {96, '1x'};
g.RowSpacing    = 4;
g.ColumnSpacing = 6;
g.Padding       = [8 4 8 2];

% Row 1: colormap "title stripe" doubling as the legend colour bar (full width).
% For a split map this has a grey GAP in the middle (the un-thresholded near-zero
% band that is NOT coloured on the brain), so it reads: extreme-neg .. near-0-neg
% [gap] near-0-pos .. extreme-pos.
stripe = make_colormap_strip(g, stripe_colormap(args, layer));
stripe.Layout.Row = 1; stripe.Layout.Column = [1 2];
stripe.Tag = sprintf('stripe_%d', k);

% Row 2: numeric legend labels under the stripe. Split maps show FOUR values
% aligned to the bar (neg extreme / neg near-0 [gap] pos near-0 / pos extreme);
% single ramps show the two ends. Reflects the layer's current cmaprange.
lg = make_legend_labels(g, layer, k);
lg.Layout.Row = 2; lg.Layout.Column = [1 2];

% Row 3: Opacity
uilabel(g, 'Text', 'Opacity', 'FontSize', fs);
sld = uislider(g, 'Limits', [0 1], 'Value', current_opacity(args), ...
    'MajorTicks', [], 'MinorTicks', [], 'FontSize', fs, 'Tag', sprintf('opacity_%d', k), ...
    'ValueChangedFcn', @(s, ~) on_opacity(obj, k, s.Value, vname));
sld.Layout.Row = 3; sld.Layout.Column = 2;

% Row 4: threshold control (type-aware; p-values on a log scale). A numeric
% entry field sits beside the slider so a value can be typed directly — the
% slider only responds to dragging its thumb, not to clicking the track.
[thr_label, thr_value, ~, is_pval] = threshold_spec(src, layer);
uilabel(g, 'Text', thr_label, 'FontSize', fs);
tgrid = uigridlayout(g, [1 2]); tgrid.Layout.Row = 4; tgrid.Layout.Column = 2;
tgrid.ColumnWidth = {'1x', 74}; tgrid.Padding = [0 0 0 0]; tgrid.ColumnSpacing = 6;
build_threshold_controls(tgrid, obj, k, src, is_pval, thr_value, vname, fs);

% Row 5: Colors dropdown + live preview swatch
uilabel(g, 'Text', 'Colors', 'FontSize', fs);
cgrid = uigridlayout(g, [1 2]); cgrid.Layout.Row = 5; cgrid.Layout.Column = 2;
cgrid.ColumnWidth = {'1x', 54}; cgrid.Padding = [0 0 0 0]; cgrid.ColumnSpacing = 6;
dd = uidropdown(cgrid, 'Items', cmap_options, 'Value', current_colormap_label(args, cmap_options), ...
    'FontSize', fs, 'Tag', sprintf('colormap_%d', k), ...
    'ValueChangedFcn', @(d, ~) on_colormap(obj, k, d.Value, vname));
dd.Layout.Column = 1;
sw = make_colormap_strip(cgrid, swatch_colormap(args, layer));
sw.Layout.Column = 2; sw.Tag = sprintf('swatch_%d', k);
% Clicking the colour swatch opens the colour picker (also re-picks a solid colour)
sw.ImageClickedFcn = @(~, ~) pick_solid_colour(obj, k, vname);

% Row 6: Visible + per-layer Remove
uilabel(g, 'Text', 'Visible', 'FontSize', fs);
vg = uigridlayout(g, [1 2]); vg.Layout.Row = 6; vg.Layout.Column = 2;
vg.ColumnWidth = {40, '1x'}; vg.Padding = [0 0 0 0]; vg.ColumnSpacing = 6;
cb = uicheckbox(vg, 'Text', '', 'Value', true, 'Tag', sprintf('visible_%d', k), ...
    'ValueChangedFcn', @(c, ~) set_layer_visible(obj, k, c.Value));
cb.Layout.Column = 1;
rb = uibutton(vg, 'Text', 'Remove layer', 'FontSize', fs - 2, ...
    'ButtonPushedFcn', @(~, ~) remove_layer(obj, k));
rb.Layout.Column = 2;

end


function build_footer(parent, obj, fig, row) %#ok<INUSD>
fg = uigridlayout(parent, [1 2]); fg.Layout.Row = row;
fg.Padding = [4 2 4 2]; fg.ColumnSpacing = 6;
fs = base_fontsize() - 2;
uibutton(fg, 'Text', 'Re-render',     'FontSize', fs, 'ButtonPushedFcn', @(~, ~) refresh(obj));
uibutton(fg, 'Text', 'Toggle legend', 'FontSize', fs, 'ButtonPushedFcn', @(~, ~) toggle_legend(obj));
end


function toggle_legend(obj)
% Toggle colorbar legends on the montage / surface FIGURES.
%
% Figure legends are OFF by default now that the controller shows a per-layer
% legend (colour bar + numeric end labels) in its own panel. This button puts
% colorbars back on the actual figures (e.g. for a publication export) and
% removes them again. ON targets the montage figure explicitly (legend() would
% otherwise draw into the controller uifigure via gcf — the old "toggle won't
% turn back on" bug) and re-renders surfaces with colorbars; OFF removes both.
fig = obj.controller_handle;
on  = false;
if ~isempty(fig) && isgraphics(fig)
    v = getappdata(fig, 'fig_legends_on'); if ~isempty(v), on = v; end
end

if on
    remove_legend(obj);                                    % surface colorbars + tracked axes
    delete(findall(groot, 'Tag', 'fmridisp_fig_legend'));  % montage legend axes
    new_state = false;
else
    draw_montage_legend(obj);
    if ~isempty(obj.surface), composite_surfaces(obj, [], true); end   % surfaces WITH colorbars
    new_state = true;
end

if ~isempty(fig) && isgraphics(fig), setappdata(fig, 'fig_legends_on', new_state); end
end


function draw_montage_legend(obj)
% Draw the colorbar legend onto the montage FIGURE. legend() places its axes in
% the current figure (gcf); when called from a controller callback gcf is the
% controller uifigure, so we make the montage figure current first. The created
% legend axes are tagged so toggle-off can find and delete them.
montfig = [];
for i = 1:numel(obj.montage)
    ah = obj.montage{i}.axis_handles; ah = ah(ishandle(ah));
    if ~isempty(ah), montfig = ancestor(ah(1), 'figure'); break; end
end
if isempty(montfig) || ~isvalid(montfig), return, end

prev = get(groot, 'CurrentFigure');
set(groot, 'CurrentFigure', montfig);
legend(obj, 'noverbose');
if ~isempty(prev) && isvalid(prev), set(groot, 'CurrentFigure', prev); end

for k = 1:numel(obj.activation_maps)
    if isfield(obj.activation_maps{k}, 'legendhandle')
        lh = obj.activation_maps{k}.legendhandle;
        set(lh(ishandle(lh)), 'Tag', 'fmridisp_fig_legend');
    end
end
end


function pick_solid_colour(obj, k, vname)
% Open the colour picker and apply a solid colour to layer k (used by the
% 'solid colour…' dropdown item and by clicking the colour swatch).
c = uisetcolor([1 0 0], 'Choose a blob colour');
if numel(c) ~= 3, return, end
vn = obj_varname(obj); if ~isempty(vn), vname = vn; end   % echo the real variable name
set_colormap(obj, 'color', c, 'layers', k);
echo_code(vname, sprintf('set_colormap(%s, ''color'', [%g %g %g], ''layers'', %d)', vname, c(1), c(2), c(3), k));
end


% ---- threshold slider --------------------------------------------------

function build_threshold_controls(g, obj, k, src, is_pval, v0, vname, fs)
% Threshold slider (col 1) + numeric entry field (col 2), kept in sync.
% statistic_image -> LOG-scale p-value slider (slider value = log10(p)), ticks at
% .001/.005/.01/.05/.1 so .05–.1 sit close and .001–.005 spread out.
% fmri_data/region -> linear raw |x| slider anchored at 0 and the 99.9th pct of |data|.
% Both controls call on_threshold; update_controller then re-syncs them, so
% dragging the slider updates the field and typing in the field updates the slider.
tickfs = max(9, fs - 5);
sl_tag = sprintf('threshold_%d', k);
ef_tag = sprintf('threshedit_%d', k);
v0 = double(v0);   % uislider/uieditfield require double; data-derived thresholds may be single

if is_pval
    pfloor = 1e-6; ptop = 0.1;                    % extends below .001 down to ~0
    ticks  = [pfloor .001 .005 .01 .05 .1];
    labs   = {'~0', '.001', '.005', '.01', '.05', '.1'};
    lims   = log10([pfloor ptop]);
    v      = log10(min(max(v0, pfloor), ptop));
    sld = uislider(g, 'Limits', lims, 'Value', v, 'MajorTicks', log10(ticks), ...
        'MajorTickLabels', labs, 'FontSize', tickfs, 'Tag', sl_tag, ...
        'ValueChangedFcn',  @(s, ~) on_threshold(obj, k, 10 .^ s.Value, true, vname));
    sld.Layout.Column = 1;
    ef = uieditfield(g, 'numeric', 'Value', min(max(v0, pfloor), ptop), ...
        'Limits', [pfloor ptop], 'ValueDisplayFormat', '%.4g', 'FontSize', fs, ...
        'Tag', ef_tag, 'ValueChangedFcn', @(e, ~) on_threshold(obj, k, e.Value, true, vname));
    ef.Layout.Column = 2;

elseif isa(src, 'image_vector') || isa(src, 'region')
    % Extend the range to cover the current threshold (out-of-range values still
    % sit on the slider and show an empty map).
    b = max(double(raw_bound(src)), v0);
    if ~(b > 0), b = 1; end
    v = min(max(v0, 0), b);
    % NORMALISED [0 1] slider. Raw |x| ranges are tiny (e.g. [0 7e-4]); a uislider
    % with such a small Limits span does NOT reliably fire its callback on
    % drag-release (unlike the O(1)-range p-value and opacity sliders, which work).
    % So the slider runs on [0 1] and maps back to raw |x| via its UserData scale
    % (raw = position * b). The numeric field stays in raw units.
    sld = uislider(g, 'Limits', [0 1], 'Value', v / b, 'MajorTicks', [0 1], ...
        'MajorTickLabels', {'0', sprintf('%.3g', b)}, 'FontSize', tickfs, 'Tag', sl_tag, ...
        'ValueChangedFcn', @(s, ~) on_threshold(obj, k, s.Value * s.UserData, false, vname));
    sld.UserData = b;                                  % raw |x| value at slider position = 1
    sld.Layout.Column = 1;
    ef = uieditfield(g, 'numeric', 'Value', v, 'Limits', [0 Inf], ...
        'ValueDisplayFormat', '%.4g', 'FontSize', fs, 'Tag', ef_tag, ...
        'ValueChangedFcn', @(e, ~) on_threshold(obj, k, e.Value, false, vname));
    ef.Layout.Column = 2;

else
    sld = uislider(g, 'Limits', [0 1], 'Value', 0, 'Enable', 'off', 'FontSize', tickfs, 'Tag', sl_tag);
    sld.Layout.Column = 1;
    ef = uieditfield(g, 'numeric', 'Value', 0, 'Enable', 'off', 'FontSize', fs, 'Tag', ef_tag);
    ef.Layout.Column = 2;
end
end

function b = raw_bound(src)
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


% ---- colormap preview --------------------------------------------------

function h = make_colormap_strip(parent, cm)
% A uiimage (NOT uiaxes) for the colour bar / swatch: it fills its grid cell with
% no reserved tick/label margins, so the bar shows its full height instead of
% collapsing to a hairline the way an invisible uiaxes does.
h = uiimage(parent);
h.ScaleMethod = 'fill';                   % stretch the ramp to fill the cell
draw_colormap_strip(h, cm);
end

function draw_colormap_strip(h, cm)
if isempty(cm), cm = repmat([.5 .5 .5], 8, 1); end
n = size(cm, 1);
img = repmat(reshape(cm, [1 n 3]), [16 1 1]);     % a few rows tall; 'fill' stretches it
h.ImageSource = uint8(round(min(max(img, 0), 1) * 255));
end

function cm = swatch_colormap(args, layer)
% Continuous colour ramp (no gap) for the small preview swatch, taken from the
% central canlab_colormap so the swatch, the stripe, and the rendered blobs agree.
if nargin < 2, layer = []; end
cm = display_colormap(args, layer).colorbar_ramp(64);
end

function c = display_colormap(args, layer)
% The canlab_colormap the controller should PREVIEW for a layer. Usually just
% from_render_args, but a single-region indexed/atlas layer (one index value)
% previews as that region's SOLID colour rather than the whole atlas palette.
c = canlab_colormap.from_render_args(args, []);
if strcmp(c.type, 'indexed') && ~isempty(layer) && isfield(layer, 'cmaprange') ...
        && ~isempty(layer.cmaprange)
    cr = double(layer.cmaprange);
    idx = round(mean(cr));
    n = size(c.colors, 1);
    if (max(cr) - min(cr)) < 1 && idx >= 1 && idx <= n
        c = canlab_colormap.solid(c.colors(idx, :));
    end
end
end

function cm = stripe_colormap(args, layer)
% Colour bar for the legend stripe. Split maps get a grey GAP in the middle (the
% un-thresholded near-zero band that shows as gray on the brain), matching the
% 5-column label layout below: neg ramp (2 units) | gap (1) | pos ramp (2).
if nargin < 2, layer = []; end
c = display_colormap(args, layer);
if strcmp(c.type, 'split')
    nramp = 32; ngap = 8;                         % 32 | 8 | 32  ==  2 : 0.5 : 2 units
    full  = c.colorbar_ramp(2 * nramp);           % [neg(nramp); pos(nramp)]
    graycol = repmat([0.75 0.75 0.75], ngap, 1);
    cm = [full(1:nramp, :); graycol; full(nramp+1:end, :)];
else
    cm = c.colorbar_ramp(72);
end
end


% ---- in-controller legend (numeric labels under the colour stripe) -----

function lg = make_legend_labels(parent, layer, k)
% Numeric labels aligned to the colour stripe above. Five columns (2:2:1:2:2)
% match the split stripe (neg ramp | gap | pos ramp): col1 neg-extreme, col2
% neg-near-0, col3 the gap (blank), col4 pos-near-0, col5 pos-extreme. A single
% ramp uses only col1 (low) and col5 (high).
lg = uigridlayout(parent, [1 5]);
lg.ColumnWidth = {'2x', '2x', '1x', '2x', '2x'};
lg.Padding = [2 0 2 0]; lg.ColumnSpacing = 1; lg.RowSpacing = 0;
fs = max(12, base_fontsize() - 2);             % numeric legend labels (a bit larger)
aligns = {'left', 'right', 'center', 'left', 'right'};
strs = legend_label_set(layer);
for c = 1:5
    u = uilabel(lg, 'Text', strs{c}, 'FontSize', fs, 'HorizontalAlignment', aligns{c}, ...
        'Tag', sprintf('leg%d_%d', c, k));
    u.Layout.Column = c;
end
end

function strs = legend_label_set(layer)
% Five label strings (col3 is the gap, always blank) from the layer's current
% cmaprange (the mapped value range). Split cmaprange = [negExtreme negNear0
% posNear0 posExtreme] -> those four values in cols 1,2,4,5. A 2-element single
% ramp -> ends in cols 1 and 5. Rounded to 2 significant figures.
strs = {'', '', '', '', ''};
% Indexed / atlas layers colour by category, not by a continuous value range,
% so numeric end-labels (e.g. a single index +-0.1) are meaningless -> leave blank.
if isfield(layer, 'render_args') && any(strcmp(layer.render_args, 'indexmap'))
    return
end
cr = [];
if isfield(layer, 'cmaprange'), cr = layer.cmaprange; end
cr = double(cr(~isnan(cr) & ~isinf(cr)));
if numel(cr) < 2, return, end
f = @(x) num2str(round(x, 2, 'significant'));
if numel(cr) >= 4
    cr = sort(cr);                                  % [negExtreme negNear0 posNear0 posExtreme]
    strs{1} = f(cr(1)); strs{2} = f(cr(2)); strs{4} = f(cr(3)); strs{5} = f(cr(4));
else
    strs{1} = f(cr(1)); strs{5} = f(cr(end));
end
end


function update_controls_in_place(fig, obj, cmap_options)
% Re-sync each control + preview to its layer's current state (no rebuild).

% Keep the title / echoed name in sync with the real workspace variable (it may
% only become resolvable after the constructor returns and assigns it).
vn = obj_varname(obj);
if ~isempty(vn)
    setappdata(fig, 'objname', vn);
    fig.Name = sprintf('CANlab display controller  [ %s ]', vn);
end

for k = 1:numel(obj.activation_maps)
    layer = obj.activation_maps{k};
    args  = render_args_of(layer);
    src   = source_of(layer);

    sld = findobj(fig, 'Tag', sprintf('opacity_%d', k));
    if ~isempty(sld), sld.Value = current_opacity(args); end

    dd = findobj(fig, 'Tag', sprintf('colormap_%d', k));
    if ~isempty(dd), dd.Value = current_colormap_label(args, cmap_options); end

    [~, val, ~, is_pval] = threshold_spec(src, layer);
    thr = findobj(fig, 'Tag', sprintf('threshold_%d', k));
    if ~isempty(thr) && isprop(thr, 'Limits')
        if is_pval
            sv = log10(max(val, 1e-6));
            thr.Value = min(max(sv, thr.Limits(1)), thr.Limits(2));
        else
            % Normalised [0 1] |x| slider: convert the raw threshold to a slider
            % position via the stored scale, growing the scale if the threshold
            % now exceeds it (so the thumb sits at the true value).
            b = thr.UserData; if isempty(b) || ~(b > 0), b = 1; end
            if val > b
                b = val; thr.UserData = b;
                thr.MajorTickLabels = {'0', sprintf('%.3g', b)};
            end
            thr.Value = min(max(val / b, 0), 1);
        end
    end
    ef = findobj(fig, 'Tag', sprintf('threshedit_%d', k));
    if ~isempty(ef) && isprop(ef, 'Value')
        lim = [0 Inf]; if isprop(ef, 'Limits'), lim = ef.Limits; end
        ef.Value = min(max(val, lim(1)), lim(2));
    end

    a = findobj(fig, 'Tag', sprintf('stripe_%d', k));
    if ~isempty(a), draw_colormap_strip(a, stripe_colormap(args, layer)); end   % gapped bar for split
    a = findobj(fig, 'Tag', sprintf('swatch_%d', k));
    if ~isempty(a), draw_colormap_strip(a, swatch_colormap(args, layer)); end

    % Refresh the numeric legend labels from the current cmaprange (which
    % rethreshold/refresh may have changed), so labels always match the colours.
    strs = legend_label_set(layer);
    for c = 1:5, set_label_text(fig, sprintf('leg%d_%d', c, k), strs{c}); end
end
end

function set_label_text(fig, tag, txt)
h = findobj(fig, 'Tag', tag);
if ~isempty(h), h.Text = txt; end
end


% ---- callbacks (apply + echo the equivalent code line) -----------------

function on_opacity(obj, k, value, vname)
vn = obj_varname(obj); if ~isempty(vn), vname = vn; end   % echo the real variable name
set_opacity(obj, value, 'layers', k);
echo_code(vname, sprintf('set_opacity(%s, %g, ''layers'', %d)', vname, value, k));
end

function on_colormap(obj, k, choice, vname)
vn = obj_varname(obj); if ~isempty(vn), vname = vn; end   % echo the real variable name
switch choice
    case 'split (hot/cool)'
        set_colormap(obj, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''splitcolor'', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, ''layers'', %d)', vname, k));
    case 'split (mango)'
        set_colormap(obj, 'splitcolor', {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]}, 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''splitcolor'', {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]}, ''layers'', %d)', vname, k));
    case 'seafire'
        set_colormap(obj, 'splitcolor', {[.4 0 1] [.2 .9 .5] [1 .3 .2] [1 .8 .4]}, 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''splitcolor'', {[.4 0 1] [.2 .9 .5] [1 .3 .2] [1 .8 .4]}, ''layers'', %d)', vname, k));
    case 'warm (red-yellow)'
        set_colormap(obj, 'maxcolor', [1 1 0], 'mincolor', [1 0 0], 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''maxcolor'', [1 1 0], ''mincolor'', [1 0 0], ''layers'', %d)', vname, k));
    case 'cool (blue-cyan)'
        set_colormap(obj, 'maxcolor', [0 1 1], 'mincolor', [0 0 1], 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''maxcolor'', [0 1 1], ''mincolor'', [0 0 1], ''layers'', %d)', vname, k));
    case 'winter (blue-green)'
        set_colormap(obj, 'maxcolor', [0 1 .5], 'mincolor', [0 0 1], 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''maxcolor'', [0 1 0.5], ''mincolor'', [0 0 1], ''layers'', %d)', vname, k));
    case perceptual_names()
        % Perceptual / continuous LUT colormaps (viridis, inferno, turbo, ...).
        set_colormap(obj, 'colormap', canlab_perceptual_colormap(choice), 'layers', k);
        echo_code(vname, sprintf('set_colormap(%s, ''colormap'', canlab_perceptual_colormap(''%s''), ''layers'', %d)', vname, choice, k));
    case 'solid colour…'
        pick_solid_colour(obj, k, vname);
end
end

function on_threshold(obj, k, value, is_pval, vname)
% Apply the threshold on slider RELEASE / field entry (the expensive step:
% rethreshold + refresh). Wrapped so a bad value (e.g. one that leaves no voxels)
% cannot break the control. The numeric field re-syncs afterwards via
% update_controls_in_place. (Slider uses ValueChangedFcn only, matching the
% opacity slider, which fires reliably on first use.)
vn = obj_varname(obj); if ~isempty(vn), vname = vn; end   % echo the real variable name
try
    if is_pval
        rethreshold(obj, value, 'unc', 'layers', k);
        echo_code(vname, sprintf('rethreshold(%s, %.4g, ''unc'', ''layers'', %d)', vname, value, k));
    else
        rethreshold(obj, value, 'layers', k);
        echo_code(vname, sprintf('rethreshold(%s, %.4g, ''layers'', %d)', vname, value, k));
    end
catch ME
    warning('fmridisplay:controller:threshold', 'Re-threshold failed: %s', ME.message);
end
end


function echo_code(vname, callstr)
fprintf('%s = %s;\n', vname, callstr);
end

function nm = obj_varname(obj)
% Best-effort lookup of the base-workspace variable name that currently holds
% this fmridisplay handle, so echoed code (and the title) use the user's real
% variable name instead of 'obj'. Returns '' if not found (e.g. the controller
% was auto-opened from the constructor before the variable was assigned, or the
% object lives in a function workspace) — callers keep their fallback then.
nm = '';
try
    vars = evalin('base', 'who');
catch
    return
end
for i = 1:numel(vars)
    try
        val = evalin('base', vars{i});
        if isa(val, 'fmridisplay') && isscalar(val) && val == obj
            nm = vars{i};
            return
        end
    catch
        % skip variables that error on evaluation / comparison
    end
end
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
if any(strcmp(args, 'indexmap'))
    lbl = 'indexed (atlas)';
elseif any(strcmp(args, 'splitcolor'))
    sc = args{find(strcmp(args, 'splitcolor'), 1) + 1};
    if iscell(sc) && numel(sc) == 4 && isequal(sc, {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]})
        lbl = 'split (mango)';
    elseif iscell(sc) && numel(sc) == 4 && isequal(sc, {[.4 0 1] [.2 .9 .5] [1 .3 .2] [1 .8 .4]})
        lbl = 'seafire';
    else
        lbl = 'split (hot/cool)';
    end
elseif any(strcmp(args, 'color'))
    lbl = 'solid colour…';
elseif any(strcmp(args, 'colormap'))
    lbl = match_perceptual_name(args{find(strcmp(args, 'colormap'), 1) + 1});
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
val = double(val);   % uislider/uieditfield require double; data-derived thresholds may be single
end

function set_layer_visible(obj, k, tf)
% Toggle a layer's visibility on montages AND surfaces. Montage blobs hide via
% their graphics Visible; surfaces recomposite, skipping the hidden layer so
% lower layers show through.
if k < 1 || k > numel(obj.activation_maps), return, end
obj.activation_maps{k}.visible = logical(tf);
bh = obj.activation_maps{k}.blobhandles;
if ~isempty(bh)
    bh = bh(ishandle(bh));
    set(bh, 'Visible', matlab.lang.OnOffSwitchState(tf));
end
if ~isempty(obj.surface)
    composite_surfaces(obj);
end
end
