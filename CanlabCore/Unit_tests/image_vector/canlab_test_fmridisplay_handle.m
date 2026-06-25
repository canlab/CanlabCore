function tests = canlab_test_fmridisplay_handle
%CANLAB_TEST_FMRIDISPLAY_HANDLE fmridisplay handle-class + live-layer behavior.
%
% Locks in the 2026 visualization overhaul (VISUALIZATION_OVERHAUL_NOTES.md):
%   - fmridisplay is now a handle class, but the value-style call contract
%     (o2 = addblobs(o2, ...)) is preserved (the returned handle == input).
%   - blob layers retain their source data + render options, so they can be
%     re-rendered in place via refresh / rethreshold / set_colormap / set_opacity.
%
% Smoke + behavior tests: pass if the API runs and the documented state
% changes hold (handle identity, source retention, monotone re-thresholding).

tests = functiontests(localfunctions);
end


function setup(tc) %#ok<DEFNU>
close all force;
end


function teardown(tc) %#ok<DEFNU>
close all force;
end


% ------------------------------------------------------------------------
% Helpers
% ------------------------------------------------------------------------

function o2 = build_display_with_statimage_layer(p_thresh)
% fmridisplay with one axial montage and a blob layer whose SOURCE is a
% statistic_image (so rethreshold can move the threshold up and down).
if nargin < 1, p_thresh = 0.01; end
t  = canlab_get_sample_thresholded_t(p_thresh);
o2 = fmridisplay;
o2 = montage(o2);
o2 = addblobs(o2, t, 'noverbose');
end

function n = displayed_voxels(o2, k)
n = sum(o2.activation_maps{k}.mapdata(:) ~= 0);
end


% ------------------------------------------------------------------------
% Tests
% ------------------------------------------------------------------------

function test_public_vs_internal_method_surface(tc)
% Internal/helper methods are Hidden (off the user-facing surface) but still
% callable; the user-facing methods stay public.
m = methods('fmridisplay');
internal = {'activate_figures', 'prune_dead_views', 'refresh', ...
            'render_layer_surfaces', 'update_controller'};
for i = 1:numel(internal)
    tc.verifyFalse(ismember(internal{i}, m), ...
        sprintf('%s should be Hidden from the public method surface', internal{i}));
end
publicm = {'addblobs', 'controller', 'montage', 'removeblobs', 'rethreshold', ...
           'set_colormap', 'set_opacity', 'surface', 'remove_legend', 'addpoints', 'legend'};
for i = 1:numel(publicm)
    tc.verifyTrue(ismember(publicm{i}, m), ...
        sprintf('%s should be a public method', publicm{i}));
end
% Hidden (not private): the internal methods must remain callable.
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2); o2 = addblobs(o2, t, 'noverbose');
tc.verifyWarningFree(@() refresh(o2));
end


function test_is_handle_class(tc)
% fmridisplay must inherit from handle for stateful, single-source-of-truth use.
mc = meta.class.fromName('fmridisplay');
super = {mc.SuperclassList.Name};
tc.verifyTrue(any(strcmp(super, 'handle')), 'fmridisplay should be a handle subclass');
end


function test_handle_identity_and_backcompat(tc)
% The value-style call contract still works, and the return IS the same handle.
o2  = build_display_with_statimage_layer();
o2b = addblobs(o2, region(canlab_get_sample_thresholded_t(0.01)), 'noverbose');
tc.verifyTrue(o2b == o2, 'addblobs must return the same handle (value-style call preserved)');
tc.verifyGreaterThanOrEqual(numel(o2.activation_maps), 2, 'both addblobs calls registered layers');
end


function test_layer_retains_source_and_options(tc)
% Each layer keeps its source object + the render options + montage targeting.
o2 = build_display_with_statimage_layer();
L  = o2.activation_maps{end};
tc.verifyTrue(isfield(L, 'source_object') && ~isempty(L.source_object), 'source retained');
tc.verifyTrue(isa(L.source_object, 'statistic_image'), 'statistic_image source retained');
tc.verifyTrue(isfield(L, 'render_args'), 'render_args field present');
tc.verifyNotEmpty(L.wh_montage, 'montage targeting retained');
end


function test_refresh_redraws_layer(tc)
% refresh deletes and redraws the layer; valid blob handles exist afterward.
o2 = build_display_with_statimage_layer();
o2 = refresh(o2);
bh = o2.activation_maps{end}.blobhandles;
tc.verifyTrue(~isempty(bh) && all(ishandle(bh)), 'refresh produced valid blob handles');
end


function test_rethreshold_monotone(tc)
% rethreshold from a statistic_image source: tighter => fewer voxels,
% looser => more voxels (re-thresholding works in both directions).
o2 = build_display_with_statimage_layer(0.01);
k  = numel(o2.activation_maps);
v_init = displayed_voxels(o2, k);

o2 = rethreshold(o2, 0.001, 'unc');
v_tight = displayed_voxels(o2, k);

o2 = rethreshold(o2, 0.05, 'unc');
v_loose = displayed_voxels(o2, k);

tc.verifyLessThanOrEqual(v_tight, v_init,  'tightening should not increase displayed voxels');
tc.verifyGreaterThanOrEqual(v_loose, v_init, 'loosening should not decrease displayed voxels');
tc.verifyLessThan(v_tight, v_loose, 'tighter threshold shows strictly fewer voxels than looser');
% and the layer is actually redrawn
tc.verifyTrue(all(ishandle(o2.activation_maps{k}.blobhandles)), 'rethreshold left valid handles');
end


function test_set_opacity_runs_and_records(tc)
% set_opacity runs warning-free and records the opacity in render_args.
o2 = build_display_with_statimage_layer();
tc.verifyWarningFree(@() set_opacity(o2, 0.4));
args = o2.activation_maps{end}.render_args;
whv  = find(strcmp(args, 'transvalue'));
tc.verifyNotEmpty(whv, 'transvalue recorded by set_opacity');
tc.verifyEqual(args{whv(1) + 1}, 0.4, 'recorded opacity value');
end


function test_set_colormap_runs(tc)
% set_colormap runs warning-free and re-renders with valid handles.
o2 = build_display_with_statimage_layer();
tc.verifyWarningFree(@() set_colormap(o2, 'color', [0 1 0]));
tc.verifyTrue(all(ishandle(o2.activation_maps{end}.blobhandles)), 'colors re-rendered with valid handles');
end


function test_controller_opens_and_binds(tc)
% The controller widget opens, returns a uifigure, and stores a back-pointer
% (reference) to the SAME instance in appdata. Skip cleanly when no display.
tc.assumeTrue(usejava('desktop') || (usejava('jvm') && feature('ShowFigureWindows')), ...
    'controller widget requires an interactive figure window');
o2  = build_display_with_statimage_layer();
fig = controller(o2);
tc.addTeardown(@() close(fig));
tc.verifyClass(fig, 'matlab.ui.Figure');
tc.verifyTrue(getappdata(fig, 'fmridisplay_obj') == o2, 'controller binds the same handle');
end


function test_canonical_callsite_compat(tc)
% The canonical, most-common workflow must keep working unchanged.
t  = canlab_get_sample_thresholded_t(0.005);
o2 = canlab_results_fmridisplay(t, 'noverbose');
tc.verifyClass(o2, 'fmridisplay');
o2 = addblobs(o2, region(t), 'noverbose');
o2 = removeblobs(o2);
tc.verifyClass(o2, 'fmridisplay');
end


% ------------------------------------------------------------------------
% Surface integration (unified views): surfaces register on the same object,
% blobs render across montages + surfaces together, later-added views pull in
% existing blobs, and refresh/remove propagate to surfaces.
% ------------------------------------------------------------------------

function o2 = build_montage_surface_blobs(tc, with_legend)
% montage + one surface (added BEFORE blobs) + blobs. Skips cleanly if no GL.
% Surface colorbars are OFF by default now; pass with_legend=true to force them
% on (composite_surfaces with show_legend) for the legend-specific tests.
if nargin < 2, with_legend = false; end
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2);
try
    o2 = surface(o2);
catch ME
    tc.assumeFail(['surface needs a graphics environment: ' ME.message]);
end
o2 = addblobs(o2, t, 'noverbose');
if with_legend, o2 = composite_surfaces(o2, [], true); end   % draw surface colorbars
end


function test_surface_keeps_same_handle(tc)
% Adding a surface must register on the SAME object, not replace it.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
o2 = fmridisplay; o2 = montage(o2);
try, o2b = surface(o2); catch ME, tc.assumeFail(ME.message); end
tc.verifyTrue(o2b == o2, 'surface must return the same handle');
tc.verifyEqual(numel(o2.surface), 1, 'surface registered on the object');
end


function test_surface_addbrain_passthrough(tc)
% surface() is a pass-through to addbrain: a bare direction token works (not
% just 'direction', X), and ANY addbrain surface/composite keyword is usable,
% including the brainstem/caudate composites now centralized in addbrain.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
o2 = fmridisplay; o2 = montage(o2);

% (1) bare token == 'direction', token
try, o2 = surface(o2, 'thalamus'); catch ME, tc.assumeFail(ME.message); end
tc.verifyEqual(numel(o2.surface), 1, 'bare-token surface registered a view');
tc.verifyEqual(o2.surface{end}.direction, 'thalamus', 'bare token recorded as direction');

% (2) addbrain composite keyword passes through and yields a multi-patch handle
try, o2 = surface(o2, 'brainstem left'); catch ME, tc.assumeFail(ME.message); end
tc.verifyEqual(numel(o2.surface), 2, 'composite surface registered a second view');
tc.verifyGreaterThan(numel(o2.surface{end}.object_handle), 1, ...
    'brainstem composite produced multiple surface patches');
end


function test_surface_unknown_keyword_informative_error(tc)
% An unrecognized direction must raise an informative error (naming the token
% + pointing to addbrain), not addbrain's terse 'Unknown method.'.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
o2 = fmridisplay; o2 = montage(o2);
err = [];
try, surface(o2, 'cutaway_left'); catch err, end %#ok<CTCH>
tc.verifyNotEmpty(err, 'unknown direction should error');
tc.verifyEqual(err.identifier, 'fmridisplay:surface:unknownDirection');
tc.verifySubstring(err.message, 'cutaway_left', 'error names the offending keyword');
tc.verifySubstring(err.message, 'addbrain', 'error points to addbrain');
end


function test_surface_medial_flips_azimuth(tc)
% 'medial' mirrors the camera azimuth 180 deg relative to the lateral default
% that addbrain sets (medial = lateral azimuth + 180).
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
oL = fmridisplay; oL = montage(oL);
try
    oL = surface(oL, 'direction', 'hires left', 'orientation', 'lateral');
    [azL, ~] = view(oL.surface{end}.axis_handles);
    oL = surface(oL, 'direction', 'hires left', 'orientation', 'medial');
    [azM, ~] = view(oL.surface{end}.axis_handles);
catch ME
    tc.assumeFail(ME.message);
end
tc.verifyEqual(mod(azM - azL, 360), 180, 'medial azimuth is lateral + 180');
end


function test_surface_pulls_in_existing_blobs(tc)
% A surface added AFTER blobs should render those existing blobs onto it.
% Surface colorbars are off by default now, so we verify the surface PATCH got
% true-colour vertex colouring (the actual pull-in), not a legend.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2); o2 = addblobs(o2, t, 'noverbose');
tc.verifyEmpty(o2.surface, 'no surface exists yet');
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
h = o2.surface{end}.object_handle; h = h(ishandle(h));
cdata = get(h(1), 'FaceVertexCData');
tc.verifyEqual(size(cdata, 2), 3, 'surface uses true-colour per-vertex RGB');
tc.verifyGreaterThan(size(unique(cdata, 'rows'), 1), 1, ...
    'blobs coloured some surface vertices (pull-in happened)');
end


function test_controller_threshold_slider_and_new_colormaps(tc)
% The redesigned controller uses a p-value threshold slider and offers the
% split (mango) / seafire colormaps.
tc.assumeTrue(usejava('desktop') || (usejava('jvm') && feature('ShowFigureWindows')), ...
    'controller widget requires an interactive figure window');
t   = canlab_get_sample_thresholded_t(0.01);
o2  = montage(t);
fig = controller(o2);
tc.addTeardown(@() delete(fig(isvalid(fig))));
ts = findobj(fig, 'Tag', 'threshold_1');
tc.verifyClass(ts, 'matlab.ui.control.Slider');
% p-value slider is LOG scale extending below .001 down toward 0 (1e-6).
tc.verifyEqual(ts.Limits, log10([1e-6 .1]), 'AbsTol', 1e-9, 'log p-value range to ~0');
tc.verifyEqual(ts.MajorTickLabels(:)', {'~0','.001','.005','.01','.05','.1'}, 'p-value tick labels');
dd = findobj(fig, 'Tag', 'colormap_1');
tc.verifyTrue(all(ismember({'split (mango)', 'seafire'}, dd.Items)), 'mango/seafire offered');
end


function test_legend_after_set_colormap(tc)
% legend(obj) must not error after set_colormap switches a layer to a single
% colormap (warm/cool/solid). Previously the stale 2-row split mincolor made
% legend take the split path and index past the 2-element cmaprange.
tc.assumeTrue(usejava('jvm'), 'legend rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2); o2 = addblobs(o2, t, 'noverbose');
o2 = set_colormap(o2, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);   % warm (single ramp)
tc.verifyEqual(size(o2.activation_maps{end}.mincolor, 1), 1, 'mincolor kept single-row in sync');
tc.verifyWarningFree(@() legend(o2));
o2 = set_colormap(o2, 'color', [1 .4 .9]);                          % solid
tc.verifyWarningFree(@() legend(o2));
end


function test_surface_multilayer_compositing(tc)
% Multiple layers composite on surfaces (top wins per vertex; lower layer shows
% through). Removing the top layer recomposites and clears its colour.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t_broad  = canlab_get_sample_thresholded_t(0.05);
t_narrow = canlab_get_sample_thresholded_t(0.001);
o2 = fmridisplay; o2 = montage(o2);
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
o2 = addblobs(o2, t_broad, 'color', [0 1 0], 'noverbose');   % layer 1: solid green, broad
o2 = addblobs(o2, t_narrow, 'noverbose');                    % layer 2: split, narrow, on top

fvc   = o2.surface{1}.object_handle(1).FaceVertexCData;
green = fvc(:, 2) > 0.6 & fvc(:, 1) < 0.4 & fvc(:, 3) < 0.4;          % layer 1
split = (max(fvc, [], 2) - min(fvc, [], 2)) > 0.2 & ~green;          % layer 2
tc.verifyGreaterThan(sum(green), 100, 'lower layer shows through (green)');
tc.verifyGreaterThan(sum(split), 100, 'top layer is visible (split colours)');

o2 = remove_layer(o2, 2);                                            % drop the top layer
fvc = o2.surface{1}.object_handle(1).FaceVertexCData;
g2  = fvc(:, 2) > 0.6 & fvc(:, 1) < 0.4 & fvc(:, 3) < 0.4;
split2 = (max(fvc, [], 2) - min(fvc, [], 2)) > 0.2 & ~g2;
tc.verifyEqual(sum(split2), 0, 'removing the top layer clears its surface colour');
tc.verifyGreaterThan(sum(g2), 100, 'the remaining layer is still composited');
end


function test_surface_layer_visibility(tc)
% Hiding a layer recomposites the surface, skipping it so the lower layer shows
% through; re-showing brings it back.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t1 = canlab_get_sample_thresholded_t(0.05);
t2 = canlab_get_sample_thresholded_t(0.001);
o2 = fmridisplay; o2 = montage(o2);
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
o2 = addblobs(o2, t1, 'color', [0 1 0], 'noverbose');   % layer 1 green
o2 = addblobs(o2, t2, 'noverbose');                     % layer 2 split, on top
isgreen = @(f) f(:, 2) > .6 & f(:, 1) < .4 & f(:, 3) < .4;
issplit = @(f) (max(f, [], 2) - min(f, [], 2)) > .2 & ~isgreen(f);

fvc = o2.surface{1}.object_handle(1).FaceVertexCData;
tc.verifyGreaterThan(sum(issplit(fvc)), 100, 'layer 2 visible initially');

o2.activation_maps{2}.visible = false;
o2 = composite_surfaces(o2);
fvc2 = o2.surface{1}.object_handle(1).FaceVertexCData;
tc.verifyEqual(sum(issplit(fvc2)), 0, 'hidden layer 2 is not shown on the surface');
tc.verifyGreaterThan(sum(isgreen(fvc2)), 100, 'layer 1 shows through where layer 2 was');

o2.activation_maps{2}.visible = true;
o2 = composite_surfaces(o2);
tc.verifyGreaterThan(sum(issplit(o2.surface{1}.object_handle(1).FaceVertexCData)), 100, 'layer 2 returns when re-shown');
end


function test_surface_first_render_uncolored_is_gray(tc)
% Regression: on the FIRST true-colour render onto a fresh (solid-FaceColor)
% surface, uncoloured areas must be proper gray, not near-black (the fresh
% surface's FaceVertexCData can be a stale dark value).
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2);
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
o2 = addblobs(o2, t, 'noverbose');                       % FIRST render
fvc = o2.surface{1}.object_handle(1).FaceVertexCData;
unc = fvc((max(fvc, [], 2) - min(fvc, [], 2)) < 0.05, :);  % uncoloured (gray) vertices
tc.verifyGreaterThan(mean(unc(:)), 0.3, 'uncoloured areas are bright gray, not near-black');
end


function test_surface_truecolor_via_central_map(tc)
% Surfaces render in true-colour RGB (N x 3 FaceVertexCData) computed from the
% central canlab_colormap, so colours match it (e.g. warm = red->yellow, no blue).
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2);
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
o2 = addblobs(o2, t, 'noverbose');
s = o2.surface{1}.object_handle(1);
tc.verifyEqual(size(s.FaceVertexCData, 2), 3, 'surface uses true-colour N x 3 RGB');
tc.verifyEqual(s.CDataMapping, 'direct');
o2 = set_colormap(o2, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);   % warm
fvc = o2.surface{1}.object_handle(1).FaceVertexCData;
colored = (max(fvc, [], 2) - min(fvc, [], 2)) > 0.1;
tc.verifyLessThan(mean(fvc(colored, 3)), 0.2, 'warm true-colour surface has little blue');
end


function test_surface_colormap_follows_maxmin_color(tc)
% maxcolor/mincolor (warm/cool/winter) must actually colour the surface, not
% fall back to render_on_surface's default split hot/cool. Warm = red->yellow
% with little blue; switching to cool changes the surface colormap.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2);
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
o2 = addblobs(o2, t, 'noverbose');
ax = ancestor(o2.surface{1}.object_handle(1), 'axes');

o2 = set_colormap(o2, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);   % warm
cm_warm = colormap(ax);
top = cm_warm(end-50:end, :);                                       % blob end of the map
tc.verifyLessThan(mean(top(:, 3)), 0.2, 'warm surface map has little blue (not the default split)');

o2 = set_colormap(o2, 'maxcolor', [0 1 1], 'mincolor', [0 0 1]);   % cool
tc.verifyFalse(isequal(cm_warm, colormap(ax)), 'surface colormap updates when layer colors change');
end


function test_foursurfaces_uses_dedicated_figure(tc)
% A multi-surface layout must not draw over the montage figure when the user
% forgets to open a figure first.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
fm = figure; o2 = montage(t);
tc.addTeardown(@() close(fm(isvalid(fm))));
try, o2 = surface(o2, 'foursurfaces'); catch ME, tc.assumeFail(ME.message); end   % no figure first
sf = ancestor(o2.surface{end}.object_handle(1), 'figure');
tc.addTeardown(@() close(sf(isvalid(sf) & sf ~= fm)));
tc.verifyNotEqual(sf, fm, 'four-surface layout uses its own figure, not the montage');
tc.verifyEqual(numel(o2.surface), 4, 'four surfaces registered');
end


function test_foursurfaces_adds_four_views_same_handle(tc)
% A multi-surface keyword adds four registered views to the same object.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
o2 = fmridisplay; o2 = montage(o2);
n0 = numel(o2.surface);
try, o2b = surface(o2, 'foursurfaces'); catch ME, tc.assumeFail(ME.message); end
tc.verifyTrue(o2b == o2, 'foursurfaces must not replace the object');
tc.verifyEqual(numel(o2.surface), n0 + 4, 'four surface views added');
end


function test_repeated_addblobs_then_removeblobs_clears_surface(tc)
% Regression: re-painting a surface (addblobs twice, or addblobs after adding
% a new surface) must not poison the saved gray, so removeblobs still restores
% anatomy. Previously render_on_surface saved the already-painted blob colors
% as the restore data, and the surface could never be cleared.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2);
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
o2 = addblobs(o2, t, 'noverbose');        % paint #1
o2 = addblobs(o2, t, 'noverbose');        % paint #2 (re-paint same surface)
% With blobs the surface has many distinct vertex colors:
blobs = o2.surface{1}.object_handle(1).FaceVertexCData;
n_blobs = size(unique(blobs, 'rows'), 1);
tc.assumeGreaterThan(n_blobs, 10, 'surface should be multi-colored while blobs are shown');
o2 = removeblobs(o2);
after = o2.surface{1}.object_handle(1).FaceVertexCData;
n_after = size(unique(after, 'rows'), 1);
% After removeblobs the surface must collapse back to ~uniform gray. Before the
% fix it kept the blob colors (n_after stayed large) because the re-paint had
% poisoned the saved restore data.
tc.verifyLessThan(n_after, n_blobs, 'removeblobs cleared the surface blobs');
tc.verifyLessThanOrEqual(n_after, 3, 'surface is essentially uniform gray after removeblobs');
end


function test_set_opacity_blends_surface(tc)
% set_opacity blends a layer's surface colours toward what's underneath (gray
% for a single layer), so the colours desaturate rather than the whole patch
% going translucent.
o2 = build_montage_surface_blobs(tc);
s  = o2.surface{1}.object_handle(1);
fvc0    = s.FaceVertexCData;
sat0    = max(fvc0, [], 2) - min(fvc0, [], 2);
colored = sat0 > 0.2;
tc.assumeGreaterThan(sum(colored), 100, 'layer has coloured surface vertices');
o2  = set_opacity(o2, 0.3);                            % 30% opaque -> blend toward gray
fvc1 = o2.surface{1}.object_handle(1).FaceVertexCData;
sat1 = max(fvc1(colored, :), [], 2) - min(fvc1(colored, :), [], 2);
tc.verifyLessThan(mean(sat1), mean(sat0(colored)), 'opacity blended the surface colours toward gray');
end


function test_removeblobs_survives_closed_surface_figure(tc)
% Closing a surface window must not break removeblobs: the dead view is pruned.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2);
f  = figure;
try, o2 = surface(o2); catch ME, close(f); tc.assumeFail(ME.message); end
o2 = addblobs(o2, t, 'noverbose');
close(f);                                  % user closes the surface window
o2 = removeblobs(o2);                      % must not error
tc.verifyEmpty(o2.surface, 'closed surface view pruned from the object');
end


function test_montage_layer_source_is_statistic_image(tc)
% montage(t) must retain the statistic_image as the layer source so the layer
% can be re-thresholded by p-value later.
t = canlab_get_sample_thresholded_t(0.01);
o2 = montage(t);
tc.verifyClass(o2.activation_maps{1}.source_object, 'statistic_image');
end


function test_rethreshold_raw_for_image_vector_layer(tc)
% A non-statistic_image layer (mask/mean) re-thresholds by raw value.
dat = canlab_get_sample_fmri_data();
m   = mean(dat);                            % fmri_data, no p-values
o2  = montage(m);
k   = numel(o2.activation_maps);
n0  = sum(o2.activation_maps{k}.mapdata(:) ~= 0);
o2  = rethreshold(o2, 0.5);                 % raw magnitude cutoff |x| > 0.5
n1  = sum(o2.activation_maps{k}.mapdata(:) ~= 0);
tc.verifyLessThan(n1, n0, 'raw rethreshold reduced the displayed voxels');
end


function test_controller_autoupdates_and_reflects_state(tc)
% The controller rebuilds when a layer is added, and its colormap dropdown
% reflects the layer's current colors.
tc.assumeTrue(usejava('desktop') || (usejava('jvm') && feature('ShowFigureWindows')), ...
    'controller widget requires an interactive figure window');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = montage(t);
fig = controller(o2);
tc.addTeardown(@() delete(fig(isvalid(fig))));
n_before = numel(findobj(fig, 'Type', 'uipanel'));
o2 = addblobs(o2, region(t), 'noverbose');  % add a 2nd layer
n_after = numel(findobj(fig, 'Type', 'uipanel'));
tc.verifyEqual(n_after, n_before + 1, 'controller auto-added a panel for the new layer');
% colormap state reflected (a solid color maps to the 'solid colour…' item)
o2 = set_colormap(o2, 'color', [1 0 0], 'layers', 1);
controller(o2);                              % update in place
dd = findobj(o2.controller_handle, 'Type', 'uidropdown');
tc.verifyTrue(any(strcmp({dd.Value}, 'solid colour…')), 'a dropdown reflects the solid-color layer');
end


function test_default_colormap_by_sign(tc)
% Default colormap depends on the data's sign: mixed +/- -> mango split;
% positive-only -> warm; negative-only -> cool; binary mask -> solid colour.
t = canlab_get_sample_thresholded_t(0.01);                      % mixed
o = montage(t);
tc.verifyTrue(any(strcmp(o.activation_maps{end}.render_args, 'splitcolor')), 'mixed -> mango split');

tpos = t; tpos.dat(tpos.dat < 0) = 0;
try, tpos.sig = tpos.sig & tpos.dat > 0; catch, end
a = montage(tpos); a = a.activation_maps{end}.render_args;
tc.verifyEqual(a{find(strcmp(a, 'maxcolor'), 1) + 1}, [1 1 0], 'positive-only -> warm');

tneg = t; tneg.dat(tneg.dat > 0) = 0;
try, tneg.sig = tneg.sig & tneg.dat < 0; catch, end
a = montage(tneg); a = a.activation_maps{end}.render_args;
tc.verifyEqual(a{find(strcmp(a, 'maxcolor'), 1) + 1}, [0 1 1], 'negative-only -> cool');

tbin = t; tbin.dat = single(tbin.dat ~= 0);
o = montage(tbin);
tc.verifyTrue(any(strcmp(o.activation_maps{end}.render_args, 'color')), 'binary mask -> solid colour');
end


function test_default_colormap_is_mango(tc)
% The default blob colormap is now "mango".
t  = canlab_get_sample_thresholded_t(0.01);
o2 = montage(t);
args = o2.activation_maps{1}.render_args;
sc = args{find(strcmp(args, 'splitcolor'), 1) + 1};
tc.verifyEqual(sc, {[.5 0 1] [0 .8 .3] [1 .2 1] [1 1 .3]}, 'default split colormap is mango');
end


function test_remove_layer_removes_one(tc)
% remove_layer drops a single layer (vs removeblobs which removes all).
tc.assumeTrue(usejava('jvm'), 'rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = montage(t); o2 = addblobs(o2, t, 'noverbose');
tc.verifyEqual(numel(o2.activation_maps), 2);
o2 = remove_layer(o2, 1);
tc.verifyEqual(numel(o2.activation_maps), 1, 'remove_layer removed exactly one layer');
end


function test_squeeze_figure_preserves_slice_size(tc)
% squeeze_figure removes vertical white space (shrinks the figure) while keeping
% each slice's on-screen pixel size and horizontal placement.
tc.assumeTrue(usejava('jvm'), 'montage rendering requires Java');
t   = canlab_get_sample_thresholded_t(0.01);
o2  = canlab_results_fmridisplay(t, 'noverbose');
ah  = o2.montage{end}.axis_handles; ah = ah(ishandle(ah));
a   = ah(1);
fig = ancestor(a, 'figure');
tc.addTeardown(@() close(fig(isvalid(fig))));
h0  = fig.Position(4);
pix_h0 = a.Position(4) * fig.Position(4);
pix_w0 = a.Position(3) * fig.Position(3);

o2 = squeeze_figure(o2);

tc.verifyLessThan(fig.Position(4), h0, 'figure got shorter');
tc.verifyEqual(a.Position(4) * fig.Position(4), pix_h0, 'RelTol', 0.03, 'slice pixel height preserved');
tc.verifyEqual(a.Position(3) * fig.Position(3), pix_w0, 'RelTol', 0.03, 'slice pixel width preserved');
end


function test_controller_shows_legend_labels(tc)
% The controller HOSTS the legend: numeric end labels under each layer's colour
% stripe (extremes; 0 in the centre for split maps), read from cmaprange.
tc.assumeTrue(usejava('desktop') || (usejava('jvm') && feature('ShowFigureWindows')), ...
    'controller widget requires an interactive figure window');
t   = canlab_get_sample_thresholded_t(0.01);
o2  = montage(t);
fig = controller(o2);
tc.addTeardown(@() delete(fig(isvalid(fig))));
cr  = o2.activation_maps{1}.cmaprange;
lo  = findobj(fig, 'Tag', 'leglo_1');
mid = findobj(fig, 'Tag', 'legmid_1');
hi  = findobj(fig, 'Tag', 'leghi_1');
tc.verifyNotEmpty(lo, 'low-end legend label exists');
tc.verifyNotEmpty(hi, 'high-end legend label exists');
tc.verifyEqual(lo.Text, num2str(round(min(cr), 1, 'significant')), 'low end labelled from cmaprange');
tc.verifyEqual(hi.Text, num2str(round(max(cr), 1, 'significant')), 'high end labelled from cmaprange');
if numel(cr) >= 4
    tc.verifyEqual(mid.Text, '0', 'split map labels 0 in the centre');
end
end


function test_montage_figure_legend_off_by_default(tc)
% The colorbar legend on the montage FIGURE is off by default (the controller
% carries the legend now); passing 'legend' opts back in.
t = canlab_get_sample_thresholded_t(0.01);
haslh = @(o) isfield(o.activation_maps{1}, 'legendhandle') && ...
             ~isempty(o.activation_maps{1}.legendhandle) && ...
             any(ishandle(o.activation_maps{1}.legendhandle));
o1 = montage(t);
tc.verifyFalse(haslh(o1), 'no montage figure legend by default');
o2 = montage(t, 'legend');
tc.verifyTrue(haslh(o2), '''legend'' draws the montage figure colorbar');
end


function test_surface_legend_not_on_montage(tc)
% When surface colorbars ARE drawn (legend forced on), they must be parented to
% the surface figure, not the montage figure.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
fm = figure; o2 = montage(t);
tc.addTeardown(@() close(fm(isvalid(fm))));
fs = figure;
tc.addTeardown(@() close(fs(isvalid(fs))));
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
o2 = addblobs(o2, t, 'noverbose');
o2 = composite_surfaces(o2, [], true);     % force surface colorbars on (default is off)
tc.verifyEmpty(findobj(fm, 'Type', 'colorbar'), 'no surface colorbar on the montage figure');
tc.verifyNotEmpty(findobj(fs, 'Type', 'colorbar'), 'surface colorbar lives on the surface figure');
end


function test_remove_legend_clears_surface_colorbars(tc)
% remove_legend deletes the surface colorbar legends but keeps the blobs.
o2  = build_montage_surface_blobs(tc, true);   % legend forced on
fig = ancestor(o2.surface{1}.object_handle(1), 'figure');
o2  = remove_legend(o2);
tc.verifyEmpty(findobj(fig, 'Type', 'colorbar'), 'remove_legend cleared surface colorbars');
tc.verifyNotEmpty(o2.activation_maps, 'blob layers are kept');
end


function test_rethreshold_syncs_open_controller(tc)
% A command-line rethreshold updates the open controller's threshold field.
tc.assumeTrue(usejava('desktop') || (usejava('jvm') && feature('ShowFigureWindows')), ...
    'controller widget requires an interactive figure window');
t   = canlab_get_sample_thresholded_t(0.01);
o2  = montage(t);
fig = controller(o2);
tc.addTeardown(@() delete(fig(isvalid(fig))));
o2  = rethreshold(o2, 0.001, 'unc');
thr = findobj(o2.controller_handle, 'Tag', 'threshold_1');
% threshold is a LOG-scale p-value slider, so its Value is log10(p).
tc.verifyEqual(thr.Value, log10(0.001), 'AbsTol', 1e-9, ...
    'controller threshold slider synced to command-line rethreshold');
end


function test_removeblobs_clears_layers_and_surface_legend(tc)
% removeblobs must act on all views together: drop the blob layers and remove
% the surface colorbar(s) it created. (Surface vertex colors are restored to
% gray by addbrain('eraseblobs'); that mechanism is exercised here too.)
o2 = build_montage_surface_blobs(tc, true);   % legend forced on
leg = [];
if isfield(o2.activation_maps{end}, 'legendhandle'), leg = o2.activation_maps{end}.legendhandle; end
o2 = removeblobs(o2);
tc.verifyEmpty(o2.activation_maps, 'activation maps cleared by removeblobs');
tc.verifyNotEmpty(o2.surface, 'surface views themselves are kept');
if ~isempty(leg)
    tc.verifyFalse(any(ishandle(leg)), 'surface colorbar(s) removed by removeblobs');
end
end
