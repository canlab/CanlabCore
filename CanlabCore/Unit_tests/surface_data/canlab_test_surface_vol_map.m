function tests = canlab_test_surface_vol_map
% Unit tests for fmri_surface_data volume<->surface mapping (M4).
%
% Covers vol2surf (volume -> fsaverage_164k via CBIG RF), surf2vol (fsaverage ->
% MNI fmri_data), the vol->surf->vol round-trip fidelity, and the mean /
% apply_mask / threshold overrides. Uses a synthetic smooth MNI volume + the
% vendored CBIG warps (in-repo) -- no external toolbox.
%
% Run:    runtests('canlab_test_surface_vol_map')
%
% :See also: vol2surf, surf2vol, fmri_surface_data

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'Surface_tools'));
addpath(fullfile(here, '..', '..'));
assert(~isempty(which('fmri_surface_data')), 'fmri_surface_data not on path.');
assert(~isempty(which('vol2surf')), 'vol2surf not on path.');
t.TestData.vol = local_synthetic_volume();
end


% -------------------------------------------------------------------------
function test_vol2surf_basic(t)
s = vol2surf(t.TestData.vol);
verifyEqual(t, class(s), 'fmri_surface_data');
verifyEqual(t, s.surface_space, 'fsaverage_164k');
verifyEqual(t, size(s.dat), [2*163842 1], 'Expected 2 hemispheres x 163842 vertices.');
verifyEqual(t, numel(s.brain_model.models), 2);
verifyEqual(t, s.brain_model.models{1}.struct, 'CORTEX_LEFT');
verifyEqual(t, s.brain_model.models{2}.numvert, 163842);
verifyFalse(t, any(isnan(s.dat(:))), 'No NaNs expected from interpn with extrapval 0.');
end


% -------------------------------------------------------------------------
function test_vol2surf_matches_interpn(t)
% A surface vertex value must equal the volume sampled at that vertex's RAS.
vol = t.TestData.vol;
s = vol2surf(vol);
L = load(canlab_cbig_warp_path('lh_ras')); ras = L.ras;
V = reconstruct_image(vol);
for vtest = [1 1000 80000 163842]
    voxc = vol.volInfo.mat \ [ras(:, vtest); 1];
    expected = interpn(V, voxc(1), voxc(2), voxc(3), 'linear', 0);
    verifyEqual(t, double(s.dat(vtest, 1)), expected, 'AbsTol', 1e-5, ...
        sprintf('Left vertex %d value != interpn sample.', vtest));
end
end


% -------------------------------------------------------------------------
function test_vol_surf_vol_roundtrip(t)
% Smooth data should round-trip vol->surf->vol with high correlation on the
% cortical voxels touched by the projection.
vol = t.TestData.vol;
s = vol2surf(vol);
vback = surf2vol(s);

dims = vol.volInfo.dim;
v0 = zeros(prod(dims), 1); v0(:) = double(vol.dat);
vb = zeros(prod(dims), 1); vb(vback.volInfo.wh_inmask) = double(vback.dat);
hit = vback.volInfo.image_indx;

verifyGreaterThan(t, nnz(hit), 20000, 'Expected many cortical voxels hit.');
cc = corr(v0(hit), vb(hit));
verifyGreaterThan(t, cc, 0.99, sprintf('vol->surf->vol correlation too low (r=%.4f).', cc));
verifyEqual(t, class(vback), 'fmri_data');
verifyEqual(t, vback.volInfo.dim, dims);
end


% -------------------------------------------------------------------------
function test_surf2vol_requires_fsaverage(t)
% A non-fsaverage object should be rejected by surf2vol.
o = fmri_surface_data(local_synthetic_cifti());     % fsLR_32k synthetic
verifyError(t, @() surf2vol(o), 'surf2vol:space');
end


% -------------------------------------------------------------------------
function test_vol2surf_nearest_preserves_labels(t)
% 'nearest' interpolation must keep integer label values intact.
vol = t.TestData.vol;
vol.dat = single(mod(round(double(vol.dat)*7), 4));   % integer "labels" 0..3
s = vol2surf(vol, 'interp', 'nearest');
u = unique(s.dat(:));
verifyTrue(t, all(ismember(u, [0 1 2 3])), 'Nearest interp introduced non-integer labels.');
end


% -------------------------------------------------------------------------
function test_mean_apply_mask_threshold(t)
s = vol2surf(t.TestData.vol);

% mean across a 3-map version
s3 = s; s3.dat = [s.dat, s.dat*2, s.dat*3];
s3.removed_images = false(3,1);
m = mean(s3);
verifyEqual(t, size(m.dat), [size(s.dat,1) 1]);
verifyEqual(t, double(m.dat), double(s.dat)*2, 'AbsTol', 1e-4, 'mean of [x 2x 3x] should be 2x.');

% apply_mask: keep first 1000 grayordinates
keep = false(size(s.dat,1),1); keep(1:1000) = true;
sm = apply_mask(s, keep);
verifyEqual(t, nnz(any(sm.dat~=0,2)), nnz(s.dat(1:1000)~=0), 'apply_mask zeroed wrong rows.');
verifyEqual(t, size(sm.dat), size(s.dat), 'apply_mask must keep .dat full (D5b).');

% threshold (raw, two-tailed)
st = threshold(s, 0.5);
verifyTrue(t, all(abs(st.dat(st.dat~=0)) >= 0.5), 'threshold left sub-threshold values.');
verifyEqual(t, size(st.dat), size(s.dat));
end


% =========================================================================
function vol = local_synthetic_volume()
% Smooth fmri_data on the standard MNI152 2 mm grid.
dims = [91 109 91];
tmat = [-2 0 0 92; 0 2 0 -128; 0 0 2 -74; 0 0 0 1];
[I, J, K] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
xyz = tmat * [I(:)'; J(:)'; K(:)'; ones(1, numel(I))];
f = 0.01*xyz(1,:) + 0.005*xyz(2,:) + ...
    exp(-((xyz(1,:).^2 + xyz(2,:).^2 + xyz(3,:).^2) / (2*40^2)));
iv = image_vector;
iv.volInfo = struct('mat', tmat, 'dim', dims, 'dt', [16 0], ...
    'xyzlist', [I(:) J(:) K(:)], 'nvox', prod(dims), ...
    'image_indx', true(prod(dims),1), 'wh_inmask', (1:prod(dims))', ...
    'n_inmask', prod(dims), 'fname', '');
iv.dat = single(f(:));
iv.removed_voxels = false(prod(dims), 1);
iv.removed_images = false;
vol = fmri_data(iv);
end


% -------------------------------------------------------------------------
function cii = local_synthetic_cifti()
mL = struct('struct','CORTEX_LEFT','type','surf','start',1,'count',5, ...
    'numvert',10,'vertlist',[0 2 4 6 8],'voxlist',[]);
mR = struct('struct','CORTEX_RIGHT','type','surf','start',6,'count',5, ...
    'numvert',10,'vertlist',[1 3 5 7 9],'voxlist',[]);
sform = [2 0 0 -20; 0 2 0 -22; 0 0 2 -24; 0 0 0 1];
dense = struct('type','dense','length',10,'models',{{mL,mR}}, ...
    'vol',struct('dims',[10 10 10],'sform',sform));
maps = struct('name',{'m1'},'table',{[]});
mapsdim = struct('type','scalars','length',1,'maps',maps);
cii = struct('cdata', single((1:10)'), 'intent','dscalar', 'diminfo', {{dense, mapsdim}});
end
