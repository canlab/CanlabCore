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
o3 = addblobs(o3, t.TestData.s);            % fs_LR data -> no matching mesh
% The fsaverage patch keeps its anatomy colour (data not force-painted onto it)
cc = get(o3.surface{1}.object_handle(1), 'FaceVertexCData');
verifyTrue(t, size(cc,1) ~= 32492, 'fsaverage mesh should not be painted with fs_LR data.');
close all force
end
