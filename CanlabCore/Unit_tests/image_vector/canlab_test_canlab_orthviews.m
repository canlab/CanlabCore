function tests = canlab_test_canlab_orthviews
%CANLAB_TEST_CANLAB_ORTHVIEWS Tests for the SPM-free canlab_orthviews viewer.
%
% Covers the canonical use cases promised by canlab_orthviews:
%   - default underlay loads on no-argument call
%   - thresholded t-map from the emotionreg dataset renders with blobs
%   - 'position' / 'Reposition' commands round-trip
%   - simulated click on a panel moves the crosshair
%   - SPM-compatible 'AddBlobs' API works
%   - right-click context menu is installed on each panel
%   - zoom command changes the panel axis limits
%   - figure carries the canlab_orthviews tag and the legacy Graphics marker
%
% Tests assume an interactive figure window is available; they are skipped
% in pure headless / -batch runs without a display.

tests = functiontests(localfunctions);
end


function setup(tc) %#ok<DEFNU>
close all force;
canlab_orthviews('Reset');   % no-op when no figure exists
end

function teardown(tc) %#ok<DEFNU>
canlab_orthviews('Reset');
close all force;
end


% ---------- helpers ----------------------------------------------------

function assume_display(tc)
tc.assumeTrue(usejava('desktop') ...
    || (usejava('jvm') && feature('ShowFigureWindows')), ...
    'canlab_orthviews requires an interactive figure window');
end


function fh = expect_open_figure(tc)
fh = findobj('Type','figure','Tag','canlab_orthviews');
tc.verifyNotEmpty(fh, 'canlab_orthviews figure should be open');
fh = fh(1);
end


% ---------- tests ------------------------------------------------------

function test_default_underlay_opens(tc) %#ok<DEFNU>
assume_display(tc);
tc.verifyWarningFree(@() canlab_orthviews);
fh = expect_open_figure(tc);

% Three axes in a 1x3 layout
ax = findall(fh, 'Type','axes');
tc.verifyEqual(numel(ax), 3, 'Expected 3 panels');

% Default position is the origin
xyz = canlab_orthviews('position');
tc.verifyEqual(xyz, [0 0 0], 'AbsTol', 1e-9, ...
    'Default crosshair should be [0 0 0]');
end


function test_figure_tag_and_legacy_marker(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews;
fh = expect_open_figure(tc);
tc.verifyEqual(get(fh, 'Tag'), 'canlab_orthviews');
% Legacy SPM code searches for the 'Graphics' tag; we mark UserData
% so old utilities can co-exist without us hijacking SPM's actual tag.
tc.verifyEqual(get(fh, 'UserData'), 'Graphics');
end


function test_reposition_roundtrip(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews;
xyz_in = [-12 8 -4];
canlab_orthviews('Reposition', xyz_in);
xyz_out = canlab_orthviews('position');
tc.verifyEqual(xyz_out, xyz_in, 'AbsTol', 1e-9);
end


function test_thresholded_tmap_renders(tc) %#ok<DEFNU>
assume_display(tc);

t = canlab_get_sample_thresholded_t(0.005);
tc.verifyWarningFree(@() canlab_orthviews(t));
fh = expect_open_figure(tc);

% At least one overlay image (RGB) was added on top of the underlay.
imgs = findall(fh, 'Type','image');
% Underlay imagesc objects appear as 'Type','image' too; expect more
% than the three underlay layers if blobs rendered.
tc.verifyGreaterThanOrEqual(numel(imgs), 3, ...
    'Expected at least the 3 underlay images on the figure');

% State should record the blob layer
state = getappdata(fh, 'canlab_orthviews_state');
tc.verifyNotEmpty(state.blobs, ...
    'Expected at least one blob layer after rendering t-map');
end


function test_addblobs_api(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews;
fh = expect_open_figure(tc);

% Build a tiny synthetic blob: a single voxel near (10, -8, 12) mm
XYZ_vox = [50; 60; 40];
Z       = 3.4;
M       = eye(4); M(1,4) = 10 - 50; M(2,4) = -8 - 60; M(3,4) = 12 - 40;
canlab_orthviews('AddBlobs', 1, XYZ_vox, Z, M);

state = getappdata(fh, 'canlab_orthviews_state');
tc.verifyNotEmpty(state.blobs);
mm_back = state.blobs{end}.XYZmm(:,1)';
tc.verifyEqual(mm_back, [10 -8 12], 'AbsTol', 1e-9);
end


function test_remove_blobs(tc) %#ok<DEFNU>
assume_display(tc);
t = canlab_get_sample_thresholded_t(0.005);
canlab_orthviews(t);
fh = expect_open_figure(tc);

canlab_orthviews('RemoveBlobs');
state = getappdata(fh, 'canlab_orthviews_state');
tc.verifyEmpty(state.blobs);
end


function test_zoom_changes_axis_limits(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews;
fh = expect_open_figure(tc);

% Center at origin, then zoom to 40 mm radius
canlab_orthviews('Reposition', [0 0 0]);
canlab_orthviews('zoom', 40);

state = getappdata(fh, 'canlab_orthviews_state');
xl = get(state.ax(3), 'XLim');
tc.verifyEqual(diff(xl), 80, 'AbsTol', 1e-6, ...
    'Axial-panel X span should be 2*zoom mm');

% Reset zoom
canlab_orthviews('zoom', Inf);
state = getappdata(fh, 'canlab_orthviews_state');
xl2 = get(state.ax(3), 'XLim');
tc.verifyGreaterThan(diff(xl2), diff(xl), ...
    'Full-FOV zoom should yield a wider axial X span');
end


function test_xhairs_toggle(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews;
fh = expect_open_figure(tc);

% Crosshairs are line objects on the panels
n1 = numel(findall(fh, 'Type','line'));
tc.verifyGreaterThan(n1, 0, 'Crosshairs should be drawn by default');

canlab_orthviews('xhairs', 'off');
n2 = numel(findall(fh, 'Type','line'));
tc.verifyLessThan(n2, n1, 'Turning xhairs off should remove lines');

canlab_orthviews('xhairs', 'on');
n3 = numel(findall(fh, 'Type','line'));
tc.verifyGreaterThanOrEqual(n3, n1 - 1, 'xhairs on should restore them');
end


function test_context_menu_attached(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews;
fh = expect_open_figure(tc);

state = getappdata(fh, 'canlab_orthviews_state');
for k = 1:numel(state.ax)
    cm = get(state.ax(k), 'UIContextMenu');
    tc.verifyNotEmpty(cm, sprintf('Panel %d missing context menu', k));
    items = findall(cm, 'Type','uimenu');
    tc.verifyGreaterThanOrEqual(numel(items), 4, ...
        sprintf('Panel %d context menu should have several items', k));
end
end


function test_simulated_click_navigates(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews;
fh = expect_open_figure(tc);

state = getappdata(fh, 'canlab_orthviews_state');
axi = state.ax(3);   % axial panel: (X, Y) in mm

% Invoke the same callback the axis would call on a left click,
% injecting a synthetic IntersectionPoint (CurrentPoint is read-only
% on Axes in modern MATLAB, so we can't set it from outside).
target_x = 14;
target_y = -22;
synth_evt = struct('IntersectionPoint', [target_x target_y 0]);

cb = get(axi, 'ButtonDownFcn');
tc.verifyNotEmpty(cb, 'Axis must have a ButtonDownFcn for click nav');
feval(cb, axi, synth_evt);

xyz = canlab_orthviews('position');
tc.verifyEqual(xyz(1), target_x, 'AbsTol', 1e-6);
tc.verifyEqual(xyz(2), target_y, 'AbsTol', 1e-6);
end


function test_region_input_renders(tc) %#ok<DEFNU>
assume_display(tc);
t = canlab_get_sample_thresholded_t(0.005);
r = region(t);
tc.assumeFalse(isempty(r), 'No surviving regions for region() input test');

tc.verifyWarningFree(@() canlab_orthviews(r));
fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');
tc.verifyNotEmpty(state.blobs);
end


function test_blobs_fill_voxels_not_dots(tc) %#ok<DEFNU>
% Regression: an earlier rasterizer drew one panel pixel per blob voxel,
% which made overlays look like sparse dots instead of filled regions.
% After volume-based slicing each blob voxel should cover many panel
% pixels (2 mm voxels at ~0.8 mm/pixel ≈ 6 pixels per voxel). Assert
% that the rendered overlay covers a meaningful fraction of the panel.

assume_display(tc);

t = canlab_get_sample_thresholded_t(0.005);
canlab_orthviews(t);
% Crosshair at a known activation cluster (typical right-DLPFC peak)
canlab_orthviews('Reposition', [38 22 30]);
drawnow;

fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');

% Count nonzero alpha pixels in each panel's overlay image (the second
% 'image' child is the blob layer; the first is the underlay).
coverage = zeros(1, 3);
for k = 1:3
    imgs = findall(state.ax(k), 'Type', 'image');
    for ii = 1:numel(imgs)
        ad = get(imgs(ii), 'AlphaData');
        if isempty(ad) || isscalar(ad), continue; end
        coverage(k) = coverage(k) + nnz(ad > 0);
    end
end

% Each blob voxel must cover > 1 panel pixel, otherwise we're back to
% the dots bug. Across the three panels we'd expect easily hundreds of
% lit pixels for the emotionreg cluster at this crosshair.
tc.verifyGreaterThan(max(coverage), 200, ...
    sprintf(['Blob overlay coverage too low (%s pixels per panel) — ' ...
             'voxels are probably being rendered as single dots again.'], ...
            mat2str(coverage)));
end


function test_default_background_is_white(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews;
fh = expect_open_figure(tc);
tc.verifyEqual(get(fh, 'Color'), [1 1 1], 'AbsTol', 1e-9, ...
    'Default figure background should be white');
state = getappdata(fh, 'canlab_orthviews_state');
for k = 1:numel(state.ax)
    tc.verifyEqual(get(state.ax(k), 'Color'), [1 1 1], 'AbsTol', 1e-9, ...
        sprintf('Default axis %d background should be white', k));
end
end


function test_black_keyword_sets_dark_background(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews([], 'black');
fh = expect_open_figure(tc);
tc.verifyEqual(get(fh, 'Color'), [0 0 0], 'AbsTol', 1e-9);
state = getappdata(fh, 'canlab_orthviews_state');
tc.verifyEqual(state.bgcolor, [0 0 0], 'AbsTol', 1e-9);
end


function test_backgroundcolor_keyvalue(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews([], 'backgroundcolor', [0.2 0.4 0.6]);
fh = expect_open_figure(tc);
tc.verifyEqual(get(fh, 'Color'), [0.2 0.4 0.6], 'AbsTol', 1e-6);
end


function test_nocrosshairs_keyword(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews([], 'nocrosshairs');
fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');
tc.verifyFalse(state.show_xhairs);
tc.verifyEmpty(findall(fh, 'Type', 'line'), ...
    'No line objects should remain when crosshairs are disabled');
end


function test_crosshairs_keyvalue(tc) %#ok<DEFNU>
assume_display(tc);
canlab_orthviews([], 'crosshairs', false);
fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');
tc.verifyFalse(state.show_xhairs);

canlab_orthviews([], 'crosshairs', true);
state = getappdata(fh, 'canlab_orthviews_state');
tc.verifyTrue(state.show_xhairs);
tc.verifyNotEmpty(findall(fh, 'Type', 'line'));
end


function test_crosshair_gap_at_center(tc) %#ok<DEFNU>
% Crosshairs are split into two segments per line with a ~5 mm gap.
% Each panel should contribute four short line objects (h-left,
% h-right, v-below, v-above), not two long ones.
assume_display(tc);
canlab_orthviews;
canlab_orthviews('Reposition', [0 0 0]);
fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');

per_panel = zeros(1, 3);
for k = 1:3
    per_panel(k) = numel(findall(state.ax(k), 'Type', 'line'));
end
tc.verifyTrue(all(per_panel >= 4), sprintf( ...
    ['Expected at least 4 line segments per panel (split crosshairs), ' ...
     'got %s'], mat2str(per_panel)));
end


function test_smooth_edges_default_on(tc) %#ok<DEFNU>
assume_display(tc);
t = canlab_get_sample_thresholded_t(0.005);
canlab_orthviews(t);
fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');
tc.assumeNotEmpty(state.blobs);
tc.verifyTrue(state.blobs{1}.smooth_edges, ...
    'smooth_edges should default to true');
tc.verifyEqual(state.blobs{1}.smooth_sigma_mm, 0, ...
    'smooth_sigma_mm should default to 0 (no pre-blur)');
end


function test_no_smooth_edges_keyword(tc) %#ok<DEFNU>
assume_display(tc);
t = canlab_get_sample_thresholded_t(0.005);
canlab_orthviews(t, 'no_smooth_edges');
fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');
tc.assumeNotEmpty(state.blobs);
tc.verifyFalse(state.blobs{1}.smooth_edges);

% NN sampling -> alpha is binary (only 0 or L.alpha). Trilinear sampling
% produces fractional alpha values. Check that the rendered overlay has
% only ~two distinct alpha values when smoothing is off.
imgs_per_panel = arrayfun(@(ax) findall(ax, 'Type', 'image'), ...
                         state.ax, 'UniformOutput', false);
all_unique = [];
for k = 1:numel(imgs_per_panel)
    overlays = imgs_per_panel{k};
    for ii = 1:numel(overlays)
        ad = get(overlays(ii), 'AlphaData');
        if ~isempty(ad) && ~isscalar(ad)
            all_unique = [all_unique; unique(round(double(ad(:)) * 1000) / 1000)]; %#ok<AGROW>
        end
    end
end
n_unique_alpha = numel(unique(all_unique));
tc.verifyLessThanOrEqual(n_unique_alpha, 4, ...
    sprintf(['Without smoothing, alpha should be near-binary (got %d ' ...
             'distinct values across overlays).'], n_unique_alpha));
end


function test_smooth_edges_produces_fractional_alpha(tc) %#ok<DEFNU>
% With smoothing on (the default), edges should have fractional alpha
% values strictly between 0 and 1 — many more distinct alpha levels than
% the binary case above.
assume_display(tc);
t = canlab_get_sample_thresholded_t(0.005);
canlab_orthviews(t);
fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');

n_fractional = 0;
for k = 1:numel(state.ax)
    overlays = findall(state.ax(k), 'Type', 'image');
    for ii = 1:numel(overlays)
        ad = get(overlays(ii), 'AlphaData');
        if ~isempty(ad) && ~isscalar(ad)
            n_fractional = n_fractional + nnz(ad > 0 & ad < 1);
        end
    end
end
tc.verifyGreaterThan(n_fractional, 50, ...
    sprintf(['Smoothed blobs should yield many fractional-alpha edge ' ...
             'pixels; got %d.'], n_fractional));
end


function test_smooth_sigma_prebluUR_grows_blob(tc) %#ok<DEFNU>
% 'smooth', sigma_mm pre-blurs the blob volume, so the rendered
% footprint should grow noticeably compared to the no-smoothing case.
assume_display(tc);

t = canlab_get_sample_thresholded_t(0.005);

canlab_orthviews(t);
fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');
n_nosmooth = nnz(state.blobs{1}.vol ~= 0);

canlab_orthviews(t, 'smooth', 4, 'replaceblobs');
state = getappdata(fh, 'canlab_orthviews_state');
n_smoothed = nnz(state.blobs{1}.vol ~= 0);

tc.verifyGreaterThan(n_smoothed, n_nosmooth, ...
    'Pre-blurred blob should occupy more voxels than the original');
tc.verifyEqual(state.blobs{1}.smooth_sigma_mm, 4, 'AbsTol', 1e-9);
end


function test_callbacks_survive_reused_figure(tc) %#ok<DEFNU>
% Regression: anonymous-closure callbacks referencing local functions
% went stale on reused figures, producing
% "Unable to find function @(src,evt)on_axis_click(src,evt,fh,k)".
% After the fix, callbacks must (a) be non-anonymous handles, and
% (b) actually fire when invoked on a figure that has been reused
% across calls.

assume_display(tc);

canlab_orthviews;                       % first call -> new figure
t = canlab_get_sample_thresholded_t(0.005);
canlab_orthviews(t);                    % second call -> reuses figure
fh = expect_open_figure(tc);

state = getappdata(fh, 'canlab_orthviews_state');

% Verify each panel has a non-anonymous handle. Anonymous closures
% report 'anonymous' in their functions() metadata; named handles
% report 'simple' or 'scopedfunction'.
for k = 1:numel(state.ax)
    cb = get(state.ax(k), 'ButtonDownFcn');
    tc.verifyClass(cb, 'function_handle', ...
        sprintf('Panel %d should have a function handle callback', k));
    info = functions(cb);
    tc.verifyNotEqual(info.type, 'anonymous', sprintf( ...
        'Panel %d callback is anonymous (will go stale across runs)', k));
end

% Actually fire the callback on the reused figure and verify it works.
axi = state.ax(3);
evt = struct('IntersectionPoint', [-6 4 0]);
feval(get(axi, 'ButtonDownFcn'), axi, evt);
xyz = canlab_orthviews('position');
tc.verifyEqual(xyz(1), -6, 'AbsTol', 1e-6);
tc.verifyEqual(xyz(2),  4, 'AbsTol', 1e-6);
end


function test_statistic_image_threshold_respected(tc) %#ok<DEFNU>
% Confirm that .sig from threshold() is honored end-to-end. Compare blob
% voxel count to the number of surviving (.sig == 1) voxels in the
% statistic_image, not to the unthresholded .dat.

assume_display(tc);

t = canlab_get_sample_thresholded_t(0.005);
canlab_orthviews(t);

fh = expect_open_figure(tc);
state = getappdata(fh, 'canlab_orthviews_state');
tc.verifyNotEmpty(state.blobs);

% Total survivors after thresholding (number of nonzero voxels in the
% rendered volume, summed across split layers).
n_sig = sum(t.sig(:) ~= 0);
n_blob = 0;
for k = 1:numel(state.blobs)
    V = state.blobs{k}.vol;
    n_blob = n_blob + nnz(V ~= 0 & isfinite(V));
end

% The blob volume is rebuilt from reconstruct_image, which multiplies
% .dat by .sig — counts should match (allow small slack for any zero
% .dat values that incidentally pass threshold).
tc.verifyGreaterThanOrEqual(n_blob, max(1, n_sig - 5));
tc.verifyLessThanOrEqual(n_blob, n_sig + 5);
tc.verifyLessThan(n_blob, sum(t.sig(:) == 0) + n_sig + 5, ...
    'Blob count should reflect thresholded voxels, not the full image.');
end
