function tests = canlab_test_surface_fmridisplay
% Unit tests for fmri_surface_data <-> fmridisplay integration (visualization
% harmonization, Level 2). A surface object is added to a managed fmridisplay as
% a surface-native layer via addblobs, painted directly on matching cortical
% meshes, and driven by set_colormap / removeblobs like a volume layer.
%
% Run:    runtests('canlab_test_surface_fmridisplay')
%
% :See also: add_surface_blobs, render_layer_surfaces, addblobs, fmridisplay

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'Surface_tools'));
addpath(fullfile(here, '..', '..'));
assert(~isempty(which('fmri_surface_data')), 'fmri_surface_data not on path.');
assert(~isempty(which('add_surface_blobs')), 'add_surface_blobs not on path.');
f = which('S1200.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii');
assert(~isempty(f), 'bundled myelin sample not on path.');
t.TestData.s = fmri_surface_data(f);     % fs_LR-32k cortex
t.TestData.figvis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
end

function teardownOnce(t)
close all force
set(0, 'DefaultFigureVisible', t.TestData.figvis);
end


% -------------------------------------------------------------------------
function test_isempty_fixed(t)
% A populated surface object must not report empty (cortex-only volInfo is empty).
verifyFalse(t, isempty(t.TestData.s), 'A populated cortex object must not be isempty.');
verifyTrue(t, isempty(fmri_surface_data), 'An empty object must be isempty.');
end


% -------------------------------------------------------------------------
function test_addblobs_paints_surface_native(t)
o2 = fmridisplay;
o2 = surface(o2, 'hcp inflated left');
o2 = surface(o2, 'hcp inflated right');
o2 = addblobs(o2, t.TestData.s, 'colormap', 'hot', 'cmaprange', [1 1.8]);

% A surface-native layer was created
verifyEqual(t, numel(o2.activation_maps), 1);
verifyTrue(t, isa(o2.activation_maps{1}.source_surface, 'fmri_surface_data'));

% Both hemispheres painted with per-vertex truecolor (in-data verts non-gray)
cL = get(o2.surface{1}.object_handle(1), 'FaceVertexCData');
cR = get(o2.surface{2}.object_handle(1), 'FaceVertexCData');
verifyEqual(t, size(cL), [32492 3], 'Left mesh must be painted per-vertex.');
verifyEqual(t, size(cR), [32492 3], 'Right mesh must be painted per-vertex.');
verifyGreaterThan(t, nnz(~all(abs(cL-0.5) < 1e-6, 2)), 20000, 'Many left vertices should be colored.');
close all force
end


% -------------------------------------------------------------------------
function test_set_colormap_and_removeblobs(t)
o2 = fmridisplay;
o2 = surface(o2, 'hcp inflated left');
o2 = addblobs(o2, t.TestData.s, 'colormap', 'hot');
c1 = get(o2.surface{1}.object_handle(1), 'FaceVertexCData');

% set_colormap must recolor the surface (managed-display integration)
o2 = set_colormap(o2, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);
c2 = get(o2.surface{1}.object_handle(1), 'FaceVertexCData');
verifyFalse(t, isequal(c1, c2), 'set_colormap must propagate to surface layers.');

% removeblobs must restore the anatomy (gray)
o2 = removeblobs(o2);
c3 = get(o2.surface{1}.object_handle(1), 'FaceVertexCData');
verifyTrue(t, all(all(abs(c3 - 0.5) < 1e-6)) || size(unique(c3,'rows'),1) <= 1, ...
    'removeblobs must restore gray on surfaces.');
close all force
end


% -------------------------------------------------------------------------
function test_surface_returns_managed_and_paints_correct_hemis(t)
% surface(obj) now returns a stateful fmridisplay with the matching native
% surfaces added and the data painted -- and each hemisphere gets ITS OWN data
% (regression: the foursurfaces_hcp patches lose their L/R tag, so hemisphere is
% resolved by x-centroid; a bug there painted left data on both hemispheres).
o2 = surface(t.TestData.s);
verifyTrue(t, isa(o2, 'fmridisplay'), 'surface(obj) must return a managed fmridisplay.');
verifyEqual(t, numel(o2.activation_maps), 1, 'Data added as one managed layer.');

r = reconstruct_image(t.TestData.s);
Ld = double(r.cortex_left(:,1)); Rd = double(r.cortex_right(:,1));
fig = ancestor(o2.surface{1}.object_handle(1), 'figure');
pp = findobj(fig, 'Type', 'patch');
checked = 0;
for hh = pp(:)'
    V = get(hh, 'Vertices'); nv = size(V,1);
    if nv ~= 32492, continue; end
    c = get(hh, 'FaceVertexCData'); if ~isequal(size(c),[nv 3]), continue; end
    bright = mean(c, 2);
    if mean(V(:,1)) < 0, own = Ld; other = Rd; else, own = Rd; other = Ld; end
    m = own ~= 0 & isfinite(own) & ~all(abs(c-0.5) < 1e-6, 2);
    verifyGreaterThan(t, corr(bright(m), own(m)), 0.8, 'Hemisphere must show its OWN data.');
    checked = checked + 1;
end
verifyEqual(t, checked, 4, 'All four cortical panels should be painted.');
close all force
end


% -------------------------------------------------------------------------
function test_surface_on_existing_display(t)
% surface(o2, obj) adds the object's native surfaces to an existing managed
% display and renders onto them (under controller management).
o2 = montage(fmridisplay, 'axial');
n_mont = numel(o2.surface);
o2 = surface(o2, t.TestData.s);
verifyGreaterThan(t, numel(o2.surface), n_mont, 'Native surfaces must be added.');
verifyEqual(t, numel(o2.activation_maps), 1, 'A managed surface-native layer is added.');
verifyNotEmpty(t, cortex_patches(o2, 32492), 'fs_LR cortical patches must exist.');
close all force
end


% -------------------------------------------------------------------------
function test_explicit_space_mismatch_warns(t)
% surface(o2, obj, 'foursurfaces_freesurfer') with fs_LR data cannot map onto the
% fsaverage meshes: a clear one-time warning, and the fsaverage patches are NOT
% painted with fs_LR data.
o2 = fmridisplay;
lastwarn('');
w = warning('off', 'all'); c = onCleanup(@() warning(w));
o2 = surface(o2, t.TestData.s, 'foursurfaces_freesurfer');
[~, wid] = lastwarn;
verifyEqual(t, wid, 'fmridisplay:render_layer_surfaces:spacemismatch', ...
    'A cross-space mismatch must warn clearly.');
% No fsaverage (163842) patch should be painted per-vertex with the fs_LR data.
for hh = cortex_patches(o2, 163842)
    cc = get(hh, 'FaceVertexCData');
    verifyTrue(t, ~isequal(size(cc), [163842 3]) || all(all(abs(cc-0.5) < 1e-6)), ...
        'fsaverage mesh must not be painted with fs_LR data.');
end
close all force
end


% -------------------------------------------------------------------------
function test_nonmatching_and_nosurface(t)
% No surface views -> a clear warning, not an error.
o2 = fmridisplay;
w = warning('off', 'fmridisplay:add_surface_blobs:nosurface');
c = onCleanup(@() warning(w));
o2 = addblobs(o2, t.TestData.s);            %#ok<NASGU>  must not error
verifyEqual(t, numel(o2.activation_maps), 1, 'Layer is still recorded.');

% Non-matching mesh (fsaverage vs fs_LR data) is skipped, not errored.
o3 = fmridisplay;
o3 = surface(o3, 'inflated left');          % fsaverage-164k (163842 verts)
w2 = warning('off', 'fmridisplay:render_layer_surfaces:spacemismatch');
c2 = onCleanup(@() warning(w2));
o3 = addblobs(o3, t.TestData.s);            % fs_LR data -> no matching mesh
% The fsaverage patch keeps its anatomy colour (data not force-painted onto it)
cc = get(o3.surface{1}.object_handle(1), 'FaceVertexCData');
verifyTrue(t, size(cc,1) ~= 32492, 'fsaverage mesh should not be painted with fs_LR data.');
close all force
end


% -------------------------------------------------------------------------
function p = cortex_patches(o2, nv)
% Cortical patches (vertex count nv) across all managed surface views.
p = gobjects(0);
for i = 1:numel(o2.surface)
    h = o2.surface{i}.object_handle;
    h = h(ishandle(h));
    for jj = 1:numel(h)
        hh = h(jj);
        if strcmp(get(hh, 'Type'), 'patch') && size(get(hh, 'Vertices'), 1) == nv
            p(end + 1) = hh; %#ok<AGROW>
        end
    end
end
end
