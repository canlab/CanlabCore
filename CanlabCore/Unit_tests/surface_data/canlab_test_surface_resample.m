function tests = canlab_test_surface_resample
% Unit tests for fmri_surface_data/resample_surface (surface-to-surface resampling).
%
% Covers fsaverage<->fs_LR (spherical barycentric / nearest), exact nested
% fsaverage downsampling, the 'list' output, auto method selection (barycentric
% for continuous, nearest for binary/label), multi-map weight reuse, and carrying
% subcortex through unchanged. The core cases use a synthetic smooth fsaverage
% object built from the bundled sphere, so they need only CanlabCore (the sphere
% assets + interpolation engine ship in-repo); the 91k subcortex case needs
% Neuroimaging_Pattern_Masks and assumeFail-skips if it is absent.
%
% Run:    runtests('canlab_test_surface_resample')
%
% :See also: resample_surface, vol2surf, surf2vol, fmri_surface_data

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
assert(~isempty(which('fmri_surface_data')), 'fmri_surface_data not on path.');
assert(~isempty(which('fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.mat')), ...
    'Bundled fsaverage<->fs_LR registration sphere not on path.');
% A synthetic SMOOTH fsaverage-164k object (value = a low-frequency function of
% sphere position), so resampling round-trips faithfully and accuracy is testable.
D = load(which('fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.mat'));
V = double(D.vertices); V = V ./ sqrt(sum(V.^2, 2));      % unit sphere
smooth = V(:, 1) + 0.5 * V(:, 2) - 0.3 * V(:, 3);         % smooth map, both hemis same
t.TestData.obj = local_make_fsaverage(smooth, smooth);
end


% -------------------------------------------------------------------------
function test_list_returns_table(t)
tbl = resample_surface(fmri_surface_data, 'list');
verifyClass(t, tbl, 'table');
verifyGreaterThanOrEqual(t, height(tbl), 5, 'Expected at least the 5 core spaces.');
verifyTrue(t, all(ismember({'fsaverage_164k', 'fs_LR_32k', 'onavg_41k'}, tbl.space)), ...
    'List must include fsaverage_164k, fs_LR_32k, and onavg_41k.');
end


% -------------------------------------------------------------------------
function test_fsaverage_to_fslr(t)
s32 = resample_surface(t.TestData.obj, 'fsLR_32k');
verifyEqual(t, s32.surface_space, 'fs_LR_32k');
verifyEqual(t, size(s32.dat, 1), 2 * 32492, 'fs_LR-32k cortex = 2 x 32492.');
% Values stay in the source range (no wild extrapolation)
src = double(t.TestData.obj.dat(:, 1)); out = double(s32.dat(:, 1));
verifyGreaterThanOrEqual(t, min(out), min(src) - 1e-3);
verifyLessThanOrEqual(t, max(out), max(src) + 1e-3);
end


% -------------------------------------------------------------------------
function test_nested_downsample_is_exact(t)
% fsaverage_164k -> fsaverage6 is the first 40962 vertices per hemisphere (exact).
s6 = resample_surface(t.TestData.obj, 'fsaverage6');
verifyEqual(t, s6.surface_space, 'fsaverage6');
verifyEqual(t, size(s6.dat, 1), 2 * 40962);
verifyEqual(t, double(s6.dat(1:40962, 1)), double(t.TestData.obj.dat(1:40962, 1)), 'AbsTol', 1e-6, ...
    'fsaverage6 must be the exact nested subset of the 164k left hemisphere.');
end


% -------------------------------------------------------------------------
function test_roundtrip_accuracy(t)
% Barycentric fsaverage -> fs_LR -> fsaverage should recover a smooth map closely.
s32   = resample_surface(t.TestData.obj, 'fsLR_32k');       % barycentric (continuous)
sback = resample_surface(s32, 'fsaverage_164k');
a = double(t.TestData.obj.dat(:, 1)); b = double(sback.dat(:, 1));
m = isfinite(a) & isfinite(b);
verifyGreaterThan(t, corr(a(m), b(m)), 0.999, 'Barycentric round-trip should be near-lossless.');
end


% -------------------------------------------------------------------------
function test_binary_defaults_to_nearest(t)
% A binary mask must stay binary (auto nearest-neighbor), not blend to [0 1].
sb = t.TestData.obj;
sb.dat = single(sb.dat(:, 1) > 0);
sb.imagetype = 'dscalar';
out = resample_surface(sb, 'fsLR_32k');
u = unique(double(out.dat(isfinite(out.dat))));
verifyTrue(t, all(ismember(u, [0 1])), 'Binary mask must resample to {0,1} (nearest).');
end


% -------------------------------------------------------------------------
function test_interp_override(t)
% Explicit 'interp','nearest' on continuous data; output differs from barycentric.
s_nn   = resample_surface(t.TestData.obj, 'fsLR_32k', 'interp', 'nearest');
s_bary = resample_surface(t.TestData.obj, 'fsLR_32k', 'interp', 'linear');
verifyEqual(t, size(s_nn.dat), size(s_bary.dat));
verifyFalse(t, isequal(s_nn.dat, s_bary.dat), 'nearest and barycentric should differ on continuous data.');
end


% -------------------------------------------------------------------------
function test_multimap_weights_reused(t)
% A multi-map object resamples every column correctly (weights built once).
sm = t.TestData.obj;
sm.dat = [sm.dat, 2 * sm.dat, -sm.dat];       % 3 maps
out = resample_surface(sm, 'fsLR_32k');
verifyEqual(t, size(out.dat, 2), 3);
% Column 2 = 2x column 1; column 3 = -column 1 (linear operator applied per map)
verifyEqual(t, double(out.dat(:, 2)), 2 * double(out.dat(:, 1)), 'AbsTol', 1e-4);
verifyEqual(t, double(out.dat(:, 3)), -double(out.dat(:, 1)), 'AbsTol', 1e-4);
end


% -------------------------------------------------------------------------
function test_onavg_space(t)
% onavg (equal-area) template: resample to onavg and back via the fs_LR frame.
if isempty(which('onavg_sphere_fsLR_lh_41k.mat'))
    t.assumeFail('onavg registration spheres not on path.');
end
son = resample_surface(t.TestData.obj, 'onavg');       % fsaverage_164k -> onavg_41k
verifyEqual(t, son.surface_space, 'onavg_41k');
verifyEqual(t, size(son.dat, 1), 2 * 40962, 'onavg den-41k cortex = 2 x 40962.');
sback = resample_surface(son, 'fsaverage_164k');
a = double(t.TestData.obj.dat(:, 1)); b = double(sback.dat(:, 1));
m = isfinite(a) & isfinite(b);
verifyGreaterThan(t, corr(a(m), b(m)), 0.99, 'onavg round-trip should recover a smooth map.');
% den-10k alias
s10 = resample_surface(t.TestData.obj, 'onavg_10k');
verifyEqual(t, size(s10.dat, 1), 2 * 10242);
end


% -------------------------------------------------------------------------
function test_subcortex_carried_through(t)
% A 91k grayordinate object: cortex resamples, subcortex passes through unchanged.
f = which('transcriptomic_gradients.dscalar.nii');
if isempty(f), t.assumeFail('transcriptomic_gradients.dscalar.nii not on path.'); end
s91 = fmri_surface_data(f);
nsub = sum(cellfun(@(m) strcmp(m.type, 'vox'), s91.brain_model.models));
out = resample_surface(s91, 'fsaverage6');
verifyEqual(t, out.surface_space, 'fsaverage6');
verifyEqual(t, sum(cellfun(@(m) strcmp(m.type, 'vox'), out.brain_model.models)), nsub, ...
    'All subcortical models must be carried through.');
sub_in  = s91.dat(2 * 32492 + 1:end, :);
sub_out = out.dat(2 * 40962 + 1:end, :);
verifyEqual(t, sub_out, sub_in, 'Subcortical data must be unchanged.');
end


% =========================================================================
function obj = local_make_fsaverage(lh, rh)
% Build a minimal dense fsaverage_164k cortex-only object from per-hemi vectors.
NV = 163842;
mL = struct('struct', 'CORTEX_LEFT',  'type', 'surf', 'start', 1, 'count', NV, ...
    'numvert', NV, 'vertlist', 0:NV - 1, 'voxlist', []);
mR = struct('struct', 'CORTEX_RIGHT', 'type', 'surf', 'start', NV + 1, 'count', NV, ...
    'numvert', NV, 'vertlist', 0:NV - 1, 'voxlist', []);
bm = struct('type', 'dense', 'length', 2 * NV, 'models', {{mL, mR}}, 'vol', []);
bm.grayordinate_type = 'cortex_only'; bm.cluster = [];
obj = fmri_surface_data('dat', single([lh(:); rh(:)]), 'brain_model', bm, ...
    'surface_space', 'fsaverage_164k', 'imagetype', 'dscalar');
end
