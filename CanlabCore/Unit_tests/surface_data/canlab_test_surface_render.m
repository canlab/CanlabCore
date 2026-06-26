function tests = canlab_test_surface_render
% Unit tests for fmri_surface_data rendering + contiguity (M5).
%
% Covers the geom loader, native surface() rendering (4-panel, direct vertex
% coloring), render_on_surface onto existing patches (native + via-volume),
% reparse_contiguous (mesh connected components), and plot QC. Figures are
% rendered offscreen and checked structurally (handles, FaceVertexCData) -- no
% visual inspection needed. Requires the bundled meshes + CBIG warps (in-repo);
% no external toolbox.
%
% Run:    runtests('canlab_test_surface_render')
%
% :See also: surface, render_on_surface, reparse_contiguous, fmri_surface_data

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'Surface_tools'));
addpath(fullfile(here, '..', '..'));
assert(~isempty(which('fmri_surface_data')), 'fmri_surface_data not on path.');
assert(~isempty(which('addbrain')), 'addbrain not on path.');
t.TestData.figvis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
% A real fs_LR dscalar (covers all vertices incl. medial wall)
f = which('transcriptomic_gradients.dscalar.nii');
assert(~isempty(f), 'transcriptomic_gradients.dscalar.nii not on path.');
t.TestData.s = fmri_surface_data(f);
end

function teardownOnce(t)
close all force;
set(0, 'DefaultFigureVisible', t.TestData.figvis);
end


% -------------------------------------------------------------------------
function test_geom_loads_for_both_spaces(t)
% Geometry is loaded (correct vertex counts) indirectly via surface() for both
% the fs_LR-32k and fsaverage-164k bundled meshes. (load_surface_geom is a
% private helper, verified through the public surface() API.)
h1 = surface(t.TestData.s);                 % fs_LR-32k
verifyEqual(t, size(get(h1.surfaces(1), 'Vertices'), 1), 32492);
close(h1.figure);

vol = local_synthetic_volume();
sf = vol2surf(vol);                          % fsaverage_164k object
h2 = surface(sf);
verifyEqual(t, size(get(h2.surfaces(1), 'Vertices'), 1), 163842);
close(h2.figure);
end


% -------------------------------------------------------------------------
function test_surface_native_render(t)
h = surface(t.TestData.s, 'which_image', 1);
verifyTrue(t, isgraphics(h.figure));
verifyEqual(t, numel(h.surfaces), 4, 'Expected 4 panels (L/R x lateral/medial).');
for k = 1:4
    c = get(h.surfaces(k), 'FaceVertexCData');
    V = get(h.surfaces(k), 'Vertices');
    verifyEqual(t, size(c), [size(V,1) 3], 'Each patch must be truecolor-colored per vertex.');
    verifyEqual(t, size(V,1), 32492, 'Native fs_LR mesh should have 32492 vertices.');
end
close(h.figure);
end


% -------------------------------------------------------------------------
function test_medial_wall_is_gray(t)
% For a 91k object (medial wall excluded), medial-wall vertices render gray.
f = which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii');
if isempty(f), t.assumeFail('No 91k dlabel on path.'); end
o = fmri_surface_data(f);
h = surface(o);
% Left hemi panel 1: medial-wall vertices (not in vertlist) should be exactly gray
lh_model = o.brain_model.models{1};
inmask = false(lh_model.numvert, 1); inmask(lh_model.vertlist + 1) = true;
c = get(h.surfaces(1), 'FaceVertexCData');
graymask = all(abs(c - 0.5) < 1e-6, 2);
% Every medial-wall vertex must be gray
verifyTrue(t, all(graymask(~inmask)), 'Medial wall vertices must render gray.');
close(h.figure);
end


% -------------------------------------------------------------------------
function test_render_on_existing_native_patch(t)
% Color an addbrain native fs_LR patch directly (no resampling).
hp = addbrain('hcp inflated left');     % 32492-vertex fs_LR mesh, Tag has 'left'
render_on_surface(t.TestData.s, hp, 'which_image', 1);
c = get(hp, 'FaceVertexCData');
verifyEqual(t, size(c), [32492 3], 'Direct native coloring should set per-vertex truecolor.');
verifyEqual(t, get(hp, 'FaceColor'), 'interp');
close all force;
end


% -------------------------------------------------------------------------
function test_render_on_mni_surface_via_volume(t)
% An fsaverage object on an arbitrary MNI surface routes through a volume.
vol = local_synthetic_volume();
sf = vol2surf(vol);
hp = addbrain('left');                  % MNI pial-ish surface (not fsaverage topology)
nverts = size(get(hp, 'Vertices'), 1);
verifyNotEqual(t, nverts, 163842, 'addbrain left should not be fsaverage topology.');
render_on_surface(sf, hp, 'clim', [-2 2]);
c = get(hp, 'FaceVertexCData');
verifyEqual(t, size(c, 1), nverts, 'Via-volume render should color all patch vertices.');
close all force;
end


% -------------------------------------------------------------------------
function test_reparse_contiguous(t)
st = threshold(t.TestData.s, 1.0, 'positive');
[st, ncl] = reparse_contiguous(st, 'which_image', 1);
verifyGreaterThan(t, ncl, 0, 'Expected at least one cluster.');
verifyEqual(t, numel(st.brain_model.cluster), size(st.dat,1));
active = st.dat(:,1) ~= 0 & ~isnan(st.dat(:,1));
verifyEqual(t, nnz(st.brain_model.cluster > 0), nnz(active), ...
    'Every active grayordinate must get a cluster label.');
verifyEqual(t, nnz(st.brain_model.cluster(~active)), 0, ...
    'Inactive grayordinates must have label 0.');
end


% -------------------------------------------------------------------------
function test_plot_runs(t)
h = plot(t.TestData.s, 'norender');
verifyTrue(t, isgraphics(h.figure));
close(h.figure);
end


% =========================================================================
function vol = local_synthetic_volume()
dims = [91 109 91];
tmat = [-2 0 0 92; 0 2 0 -128; 0 0 2 -74; 0 0 0 1];
[I, J, K] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
xyz = tmat * [I(:)'; J(:)'; K(:)'; ones(1, numel(I))];
f = sin(xyz(1,:)/30) + cos(xyz(2,:)/40);
iv = image_vector;
iv.volInfo = struct('mat', tmat, 'dim', dims, 'dt', [16 0], ...
    'xyzlist', [I(:) J(:) K(:)], 'nvox', prod(dims), ...
    'image_indx', true(prod(dims),1), 'wh_inmask', (1:prod(dims))', ...
    'n_inmask', prod(dims), 'fname', '');
iv.dat = single(f(:));
iv.removed_voxels = false(prod(dims),1); iv.removed_images = false;
vol = fmri_data(iv);
end
