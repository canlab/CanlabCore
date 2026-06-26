function tests = canlab_test_surface_space_recon
% Unit tests for fmri_surface_data spatial overrides + interop (M3).
%
% Covers compare_space (0/1/2/3 contract), reconstruct_image (cortex dense
% arrays + subcortex volume), to_fmri_data (subcortex -> fmri_data), write
% (native CIFTI round-trip, faithful + regenerated), and rebuild_like. Core
% tests use a synthetic grayordinate object (no external toolbox); one
% opportunistic test uses a real .dlabel.nii if present.
%
% Run:    runtests('canlab_test_surface_space_recon')
%
% :See also: fmri_surface_data, reconstruct_image, to_fmri_data, compare_space

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'Surface_tools'));
addpath(fullfile(here, '..', '..'));
assert(~isempty(which('fmri_surface_data')), 'fmri_surface_data not on path.');
t.TestData.tmp = tempname; mkdir(t.TestData.tmp);
t.TestData.obj = fmri_surface_data(local_synthetic_cifti());
end

function teardownOnce(t)
if isfield(t.TestData,'tmp') && exist(t.TestData.tmp,'dir'), rmdir(t.TestData.tmp,'s'); end
end


% -------------------------------------------------------------------------
function test_volinfo_subblock_populated(t)
o = t.TestData.obj;
verifyNotEmpty(t, o.volInfo, 'volInfo (subcortical sub-block) should be populated.');
verifyEqual(t, o.volInfo.dim, [10 10 10]);
verifyEqual(t, o.volInfo.n_inmask, 4, 'Subcortical sub-block should have 4 voxels.');
% 0-based CIFTI affine -> 1-based SPM .mat conversion
expected_mat = [2 0 0 -20; 0 2 0 -22; 0 0 2 -24; 0 0 0 1];
expected_mat(:,4) = expected_mat(:,4) - sum(expected_mat(:,1:3),2);
verifyEqual(t, o.volInfo.mat, expected_mat, 'Subcortical affine (1-based) incorrect.');
end


% -------------------------------------------------------------------------
function test_to_fmri_data(t)
o = t.TestData.obj;
vol = to_fmri_data(o);
verifyEqual(t, class(vol), 'fmri_data');
verifyEqual(t, size(vol.dat), [4 2], 'Subcortex should be 4 voxels x 2 maps.');
% Values equal the subcortical grayordinate rows (rows 11:14)
verifyEqual(t, double(vol.dat), double(o.dat(11:14, :)), 'AbsTol', 1e-6, ...
    'to_fmri_data values must equal the subcortical grayordinate rows.');
verifyEqual(t, vol.volInfo.dim, [10 10 10]);
end


% -------------------------------------------------------------------------
function test_reconstruct_image(t)
o = t.TestData.obj;
r = reconstruct_image(o);

% Cortex: dense [numvert x nMaps]; medial wall (non in-data verts) -> NaN
verifyEqual(t, size(r.cortex_left), [10 2]);
verifyEqual(t, sum(isnan(r.cortex_left(:,1))), 5, 'Half the vertices are medial wall (NaN).');
% In-data vertices (0-based [0 2 4 6 8] -> 1-based [1 3 5 7 9]) carry the data
verifyEqual(t, r.cortex_left([1 3 5 7 9], 1), double(o.dat(1:5, 1)), 'AbsTol', 1e-6);

% Volume present and a known voxel matches the subcortical data
verifyTrue(t, isfield(r, 'volume'));
verifyEqual(t, size(r.volume), [10 10 10 2]);
% First subcortical voxel is IJK (2,2,2) 0-based -> (3,3,3) 1-based, row 11
verifyEqual(t, r.volume(3,3,3,1), double(o.dat(11,1)), 'AbsTol', 1e-6);
end


% -------------------------------------------------------------------------
function test_compare_space_contract(t)
o = t.TestData.obj;
verifyEqual(t, compare_space(o, o), 0, 'Identical objects -> 0.');

% Selecting maps does not change the grayordinate space -> 0
verifyEqual(t, compare_space(o, get_wh_image(o, 1)), 0);

% Same space tag + layout but different in-data vertices -> 3
o3 = o;
o3.brain_model.models{1}.vertlist = [0 1 2 3 4];   % different selection
verifyEqual(t, compare_space(o, o3), 3, 'Different in-data grayordinates -> 3.');

% Different space tag -> 1
o1 = o; o1.surface_space = 'fsaverage_164k';
verifyEqual(t, compare_space(o, o1), 1, 'Different space tag -> 1.');

% Missing brain_model -> 2
oempty = fmri_surface_data;
verifyEqual(t, compare_space(o, oempty), 2, 'Missing brain_model -> 2.');
end


% -------------------------------------------------------------------------
function test_rebuild_like(t)
o = t.TestData.obj;
rb = rebuild_like(o, o.dat * 3);
verifyEqual(t, class(rb), 'fmri_surface_data');
verifyTrue(t, isequaln(rb.brain_model, o.brain_model), 'brain_model must be preserved.');
verifyEqual(t, double(rb.dat), double(o.dat)*3, 'AbsTol', 1e-3);
verifyEqual(t, numel(rb.removed_voxels), size(rb.dat,1));
% Row-count mismatch must error (geometry is fixed)
verifyError(t, @() rebuild_like(o, o.dat(1:5,:)), 'fmri_surface_data:rebuild_like:rowmismatch');
end


% -------------------------------------------------------------------------
function test_write_roundtrip_cifti(t)
o = t.TestData.obj;

% Faithful path (stashed xml present? synthetic has none -> regenerate path)
f = fullfile(t.TestData.tmp, 'syn.dscalar.nii');
write(o, f);
o2 = fmri_surface_data(f);
verifyEqual(t, size(o2.dat), size(o.dat));
verifyEqual(t, double(o2.dat), double(o.dat), 'AbsTol', 0, 'CIFTI write->read changed data.');
verifyEqual(t, numel(o2.brain_model.models), numel(o.brain_model.models));
verifyEqual(t, o2.brain_model.vol.sform, o.brain_model.vol.sform, 'Affine lost on write.');
verifyEqual(t, reshape(o2.image_names,1,[]), {'mapOne','mapTwo'}, 'Map names lost on write.');
end


% -------------------------------------------------------------------------
function test_real_dlabel_write_if_available(t)
f0 = which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii');
if isempty(f0)
    t.assumeFail('No real .dlabel.nii on path; synthetic tests cover write.');
end
o = fmri_surface_data(f0);
f = fullfile(t.TestData.tmp, 'real.dlabel.nii');
write(o, f);
o2 = fmri_surface_data(f);
verifyEqual(t, double(o2.dat), double(o.dat), 'AbsTol', 0, 'dlabel keys changed on write.');
verifyEqual(t, numel(o2.label_table), numel(o.label_table), 'label table changed on write.');
% subcortex export
vol = to_fmri_data(o);
verifyEqual(t, class(vol), 'fmri_data');
verifyGreaterThan(t, size(vol.dat,1), 1000, 'Expected thousands of subcortical voxels.');
end


% =========================================================================
function cii = local_synthetic_cifti()
mL = struct('struct','CORTEX_LEFT','type','surf','start',1,'count',5, ...
    'numvert',10,'vertlist',[0 2 4 6 8],'voxlist',[]);
mR = struct('struct','CORTEX_RIGHT','type','surf','start',6,'count',5, ...
    'numvert',10,'vertlist',[1 3 5 7 9],'voxlist',[]);
voxijk = [2 2 2; 3 2 2; 2 3 2; 2 2 3]';
mV = struct('struct','THALAMUS_LEFT','type','vox','start',11,'count',4, ...
    'numvert',NaN,'vertlist',[],'voxlist',voxijk);
sform = [2 0 0 -20; 0 2 0 -22; 0 0 2 -24; 0 0 0 1];
dense = struct('type','dense','length',14,'models',{{mL,mR,mV}}, ...
    'vol',struct('dims',[10 10 10],'sform',sform));
maps = struct('name',{'mapOne','mapTwo'},'table',{[],[]});
mapsdim = struct('type','scalars','length',2,'maps',maps);
cii = struct('cdata', single(reshape(1:28,14,2)), 'intent','dscalar', ...
    'diminfo', {{dense, mapsdim}});
end
