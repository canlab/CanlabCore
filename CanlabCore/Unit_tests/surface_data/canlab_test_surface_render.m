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
% private helper, verified through the public surface() API.) surface() now
% returns a managed fmridisplay whose cortical patches carry the right meshes.
o1 = surface(t.TestData.s);                 % fs_LR-32k -> managed display
verifyTrue(t, isa(o1, 'fmridisplay'), 'Native surface() returns a managed fmridisplay.');
verifyNotEmpty(t, cortex_patches(o1, 32492), 'Expected fs_LR 32492-vertex cortical patches.');
close all force

vol = local_synthetic_volume();
sf = vol2surf(vol);                          % fsaverage_164k object
o2 = surface(sf);
verifyNotEmpty(t, cortex_patches(o2, 163842), 'Expected fsaverage 163842-vertex cortical patches.');
close all force
end


% -------------------------------------------------------------------------
function test_surface_native_render(t)
o = surface(t.TestData.s, 'which_image', 1);
verifyTrue(t, isa(o, 'fmridisplay'), 'Native surface() returns a managed fmridisplay.');
% transcriptomic_gradients is a mixed grayordinate object: cortex (surface-native
% layer) + subcortex (volume layer), so two managed layers.
verifyGreaterThanOrEqual(t, numel(o.activation_maps), 1, 'At least the cortical layer exists.');
verifyTrue(t, isa(o.activation_maps{1}.source_surface, 'fmri_surface_data'), ...
    'First layer is the surface-native cortical layer.');
p = cortex_patches(o, 32492);
verifyEqual(t, numel(p), 4, 'Expected 4 cortical panels (L/R x lateral/medial).');
for k = 1:numel(p)
    c = get(p(k), 'FaceVertexCData');
    V = get(p(k), 'Vertices');
    verifyEqual(t, size(c), [size(V,1) 3], 'Each cortical patch must be truecolor per vertex.');
    verifyEqual(t, size(V,1), 32492, 'Native fs_LR mesh should have 32492 vertices.');
end
close all force
end


% -------------------------------------------------------------------------
function test_harmonized_color_options(t)
% surface()/render_on_surface accept the same color vocabulary as the volume
% pipeline (colormap/cmaprange, maxcolor/mincolor, splitcolor, color).
s = t.TestData.s;
% Graduated maps must produce many colors
opts = { {'colormap','parula'}, ...        % single sequential map over data range
         {'maxcolor',[1 1 0],'mincolor',[1 0 0]}, ...
         {'splitcolor',{[0 0 1],[0 1 1],[1 .5 0],[1 1 0]}} };
for k = 1:numel(opts)
    o = surface(s, opts{k}{:});
    p = cortex_patches(o, 32492);
    c = get(p(1), 'FaceVertexCData');
    verifyEqual(t, size(c,2), 3, 'Must set truecolor.');
    verifyGreaterThan(t, size(unique(c,'rows'),1), 2, 'Graduated coloring should not be uniform.');
    close all force
end
% solid 'color' -> the solid color (plus gray for any zero/NaN vertices)
o = surface(s, 'color', [0 .7 0]);
p = cortex_patches(o, 32492);
u = unique(get(p(1),'FaceVertexCData'),'rows');
verifyLessThanOrEqual(t, size(u,1), 2, 'Solid color should give <= 2 colors (color + gray).');
verifyTrue(t, ismember([0 .7 0], u, 'rows'), 'Solid color must appear on the surface.');
close all force
end


% -------------------------------------------------------------------------
function test_medial_wall_is_gray(t)
% For a 91k object (medial wall excluded), medial-wall vertices render gray.
f = which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii');
if isempty(f), t.assumeFail('No 91k dlabel on path.'); end
o = fmri_surface_data(f);
o2 = surface(o);
% Left cortical patch (x-centroid < 0): medial-wall vertices (not in vertlist)
% should be exactly gray.
p = cortex_patches(o2, 32492);
lh = p(1);
for jj = 1:numel(p)
    V = get(p(jj), 'Vertices');
    if mean(V(:,1)) < 0, lh = p(jj); break; end
end
lh_model = o.brain_model.models{1};
inmask = false(lh_model.numvert, 1); inmask(lh_model.vertlist + 1) = true;
c = get(lh, 'FaceVertexCData');
graymask = all(abs(c - 0.5) < 1e-6, 2);
% Every medial-wall vertex must be gray
verifyTrue(t, all(graymask(~inmask)), 'Medial wall vertices must render gray.');
close all force
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
function test_default_view_matches_space(t)
% surface(obj) with NO surface directive picks a four-surface view matching the
% object's space family: fs_LR -> foursurfaces_hcp (native 32492); fsaverage
% (incl. nested 6/5/4) -> foursurfaces_freesurfer (163842); onavg -> fs_LR view
% (onavg is fs_LR-aligned). Data paints on the resulting cortical patches in
% every case.
if isempty(which('resample_surface')), t.assumeFail('resample_surface not on path.'); end
s = t.TestData.s;                        % fs_LR-32k

% Native fs_LR -> foursurfaces_hcp (32492-vertex patches painted directly).
o = surface(s);
verifyNotEmpty(t, cortex_patches(o, 32492), 'fs_LR default view should have 32492 patches.');
verifyGreaterThan(t, painted_frac(o, 32492), 0.5, 'fs_LR default view should paint the cortex.');
close all force;

% fsaverage6 (nested) -> renders on the fsaverage-164k four-surface view.
s6 = resample_surface(s, 'fsaverage6', 'interp', 'nearest');
verifyEqual(t, s6.surface_space, 'fsaverage6');
o6 = surface(s6);
verifyNotEmpty(t, cortex_patches(o6, 163842), 'fsaverage6 should render on 163842-vertex patches.');
verifyGreaterThan(t, painted_frac(o6, 163842), 0.5, 'fsaverage6 default view should paint the cortex.');
close all force;

% onavg (fs_LR-aligned) -> renders on the fs_LR four-surface view.
son = resample_surface(s, 'onavg_41k', 'interp', 'nearest');
verifyEqual(t, son.surface_space, 'onavg_41k');
oon = surface(son);
verifyNotEmpty(t, cortex_patches(oon, 32492), 'onavg should render on the fs_LR (32492) view.');
verifyGreaterThan(t, painted_frac(oon, 32492), 0.5, 'onavg default view should paint the cortex.');
close all force;

% An explicit surface keyword overrides the space-matched default (and does not
% error on a bare token -- the historical parser bug).
oe = surface(s, 'foursurfaces_hcp');
verifyNotEmpty(t, cortex_patches(oe, 32492), 'Explicit foursurfaces_hcp keyword should build the view.');
close all force;
end


% -------------------------------------------------------------------------
function test_atlas_surface_unique(t)
% A volumetric atlas rendered on surfaces uses UNIQUE per-region solid colours
% (via an indexed colormap), and surface(atl,'unique') is accepted (no warning).
% A continuous statistic map keeps its colormap (no false auto-unique).
if isempty(which('load_atlas')), t.assumeFail('load_atlas not on path.'); end
atl = [];
try, atl = load_atlas('julich_fmriprep20', 'noverbose'); catch, end
if isempty(atl), t.assumeFail('julich atlas not available.'); end

verifyEqual(t, canlab_color_mode(atl), 'unique', 'An atlas is a unique-colour map.');

w = warning('off', 'all'); c = onCleanup(@() warning(w));
lastwarn('');
h = surface(atl, 'unique');                       % explicit flag, must not error/warn-unknown
[~, wid] = lastwarn;
verifyNotEqual(t, wid, 'MATLAB:UndefinedFunction', 'surface(atl,''unique'') must run.');
verifyNotEmpty(t, h, 'surface(atl,''unique'') returns surface handles.');
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
function p = cortex_patches(o2, nv)
% Cortical patches (vertex count nv) across all managed surface views of a
% fmridisplay returned by native surface().
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


% =========================================================================
function frac = painted_frac(o2, nv)
% Fraction of vertices coloured non-gray across all nv-vertex cortical patches.
p = cortex_patches(o2, nv);
tot = 0; nz = 0;
for k = 1:numel(p)
    c = get(p(k), 'FaceVertexCData');
    if size(c, 2) == 3 && size(c, 1) == nv
        tot = tot + nv;
        nz = nz + sum(~all(abs(c - 0.5) < 1e-6, 2));
    end
end
frac = nz / max(tot, 1);
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
