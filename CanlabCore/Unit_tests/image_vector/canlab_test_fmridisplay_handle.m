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
