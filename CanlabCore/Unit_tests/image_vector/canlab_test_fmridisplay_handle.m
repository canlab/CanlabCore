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

function o2 = build_montage_surface_blobs(tc)
% montage + one surface (added BEFORE blobs) + blobs. Skips cleanly if no GL.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2);
try
    o2 = surface(o2);
catch ME
    tc.assumeFail(['surface needs a graphics environment: ' ME.message]);
end
o2 = addblobs(o2, t, 'noverbose');
end


function test_surface_keeps_same_handle(tc)
% Adding a surface must register on the SAME object, not replace it.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
o2 = fmridisplay; o2 = montage(o2);
try, o2b = surface(o2); catch ME, tc.assumeFail(ME.message); end
tc.verifyTrue(o2b == o2, 'surface must return the same handle');
tc.verifyEqual(numel(o2.surface), 1, 'surface registered on the object');
end


function test_surface_pulls_in_existing_blobs(tc)
% A surface added AFTER blobs should render those existing blobs onto it.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
o2 = fmridisplay; o2 = montage(o2); o2 = addblobs(o2, t, 'noverbose');
has_legend = @() isfield(o2.activation_maps{end}, 'legendhandle') && ~isempty(o2.activation_maps{end}.legendhandle);
tc.verifyFalse(has_legend(), 'no surface legend before any surface exists');
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
tc.verifyTrue(has_legend(), 'adding a surface rendered the existing blob layer onto it');
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


function test_set_opacity_dims_surface(tc)
% set_opacity must reach surfaces (controller opacity slider): the surface
% patch FaceAlpha follows the requested value.
o2 = build_montage_surface_blobs(tc);
o2 = set_opacity(o2, 0.4);
tc.verifyEqual(o2.surface{1}.object_handle(1).FaceAlpha, 0.4, 'AbsTol', 1e-9, ...
    'set_opacity dims the surface patch');
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


function test_surface_legend_not_on_montage(tc)
% A surface colorbar must be parented to the surface figure, not the montage.
tc.assumeTrue(usejava('jvm'), 'surface rendering requires Java');
t  = canlab_get_sample_thresholded_t(0.01);
fm = figure; o2 = montage(t);
tc.addTeardown(@() close(fm(isvalid(fm))));
fs = figure;
tc.addTeardown(@() close(fs(isvalid(fs))));
try, o2 = surface(o2); catch ME, tc.assumeFail(ME.message); end
o2 = addblobs(o2, t, 'noverbose');
tc.verifyEmpty(findobj(fm, 'Type', 'colorbar'), 'no surface colorbar on the montage figure');
tc.verifyNotEmpty(findobj(fs, 'Type', 'colorbar'), 'surface colorbar lives on the surface figure');
end


function test_remove_legend_clears_surface_colorbars(tc)
% remove_legend deletes the surface colorbar legends but keeps the blobs.
o2  = build_montage_surface_blobs(tc);
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
tc.verifyEqual(thr.Value, 0.001, 'AbsTol', 1e-9, ...
    'controller threshold field synced to command-line rethreshold');
end


function test_removeblobs_clears_layers_and_surface_legend(tc)
% removeblobs must act on all views together: drop the blob layers and remove
% the surface colorbar(s) it created. (Surface vertex colors are restored to
% gray by addbrain('eraseblobs'); that mechanism is exercised here too.)
o2 = build_montage_surface_blobs(tc);
leg = [];
if isfield(o2.activation_maps{end}, 'legendhandle'), leg = o2.activation_maps{end}.legendhandle; end
o2 = removeblobs(o2);
tc.verifyEmpty(o2.activation_maps, 'activation maps cleared by removeblobs');
tc.verifyNotEmpty(o2.surface, 'surface views themselves are kept');
if ~isempty(leg)
    tc.verifyFalse(any(ishandle(leg)), 'surface colorbar(s) removed by removeblobs');
end
end
