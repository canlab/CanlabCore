function tests = canlab_test_surface_io
% Unit tests for the native CIFTI-2 / GIFTI reader+writer (M1 of fmri_surface_data).
%
% Verifies that canlab_read_gifti / canlab_write_gifti / canlab_read_cifti /
% canlab_write_cifti work with NO external toolbox (no gifti, FieldTrip,
% cifti-matlab, or Connectome Workbench) -- only core MATLAB + the JVM.
%
% Run:    runtests('canlab_test_surface_io')
% or it is auto-discovered by canlab_run_all_tests.
%
% Coverage:
%   - GIFTI .surf round-trip (shipped S1200 fs_LR-32k surface): exact verts/faces,
%     GZip and Base64 encodings.
%   - GIFTI label round-trip with a LabelTable (synthetic).
%   - CIFTI-2 dscalar round-trip (synthetic grayordinates: 2 surface models +
%     1 voxel model): cdata, BrainModel start/count/vertlist/voxlist, and the
%     subcortical affine all preserved.
%   - CIFTI-2 dlabel round-trip with a LabelTable (synthetic).
%   - If a real .dscalar.nii is on the path, an end-to-end read->write->read.
%
% :See also: canlab_read_cifti, canlab_write_cifti, canlab_read_gifti, canlab_write_gifti

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
here = fileparts(mfilename('fullpath'));
st = fullfile(here, '..', '..', 'Surface_tools');
if exist(st, 'dir'), addpath(st); end
% Sanity: the codec functions must be on the path
for fn = {'canlab_read_gifti','canlab_write_gifti','canlab_read_cifti','canlab_write_cifti'}
    assert(~isempty(which(fn{1})), 'Missing required function on path: %s', fn{1});
end
t.TestData.tmp = tempname;
mkdir(t.TestData.tmp);
end

function teardownOnce(t)
if isfield(t.TestData, 'tmp') && exist(t.TestData.tmp, 'dir')
    rmdir(t.TestData.tmp, 's');
end
end


% -------------------------------------------------------------------------
function test_gifti_surf_roundtrip(t)
surf = which('S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii');
assert(~isempty(surf), 'Shipped surface S1200.L.midthickness...surf.gii not found on path.');

g = canlab_read_gifti(surf);
verifyEqual(t, size(g.vertices), [32492 3], 'fs_LR-32k left surface should have 32492 vertices.');
verifyEqual(t, size(g.faces),    [64980 3], 'fs_LR-32k left surface should have 64980 faces.');
verifyGreaterThanOrEqual(t, min(g.faces(:)), 1, 'Faces must be 1-based for patch().');
verifyLessThanOrEqual(t, max(g.faces(:)), 32492, 'Face indices must not exceed vertex count.');

for enc = {'GZipBase64Binary', 'Base64Binary', 'ASCII'}
    f = fullfile(t.TestData.tmp, ['rt_' enc{1} '.surf.gii']);
    canlab_write_gifti(f, g, 'encoding', enc{1});
    g2 = canlab_read_gifti(f);
    tol = 0; if strcmp(enc{1}, 'ASCII'), tol = 1e-4; end   % ASCII prints %g
    verifyLessThanOrEqual(t, max(abs(g.vertices(:) - g2.vertices(:))), tol, ...
        sprintf('Vertex round-trip failed for encoding %s', enc{1}));
    verifyEqual(t, g.faces, g2.faces, ...
        sprintf('Face round-trip failed for encoding %s', enc{1}));
end
end


% -------------------------------------------------------------------------
function test_gifti_label_roundtrip(t)
n = 50;
keys = mod((0:n-1)', 4);                 % 0..3
g = struct();
g.vertices = [];  g.faces = [];
g.cdata = keys;
g.intents = {'NIFTI_INTENT_LABEL'};
g.labels = struct('key', {0 1 2 3}, ...
    'name', {'Unknown','A','B','C'}, ...
    'rgba', {[0 0 0 0],[1 0 0 1],[0 1 0 1],[0 0 1 1]});

f = fullfile(t.TestData.tmp, 'lbl.label.gii');
canlab_write_gifti(f, g);
g2 = canlab_read_gifti(f);
verifyEqual(t, g2.cdata, double(keys), 'Label keys did not round-trip.');
verifyEqual(t, numel(g2.labels), 4, 'Label table size changed.');
verifyEqual(t, g2.labels(2).name, 'A', 'Label name did not round-trip.');
verifyEqual(t, g2.labels(2).rgba, [1 0 0 1], 'Label RGBA did not round-trip.');
end


% -------------------------------------------------------------------------
function test_cifti_dscalar_roundtrip(t)
cii = local_make_synthetic_cifti('dscalar');
f = fullfile(t.TestData.tmp, 'syn.dscalar.nii');
canlab_write_cifti(f, cii);
c2 = canlab_read_cifti(f);

verifyEqual(t, c2.intent, 'dscalar');
verifyEqual(t, size(c2.cdata), size(cii.cdata), 'cdata size changed.');
verifyLessThanOrEqual(t, max(abs(cii.cdata(:) - c2.cdata(:))), 0, 'cdata values changed.');

verifyEqual(t, numel(c2.diminfo{1}.models), numel(cii.diminfo{1}.models), 'Model count changed.');
for i = 1:numel(cii.diminfo{1}.models)
    a = cii.diminfo{1}.models{i};  b = c2.diminfo{1}.models{i};
    verifyEqual(t, b.struct, a.struct);
    verifyEqual(t, b.start,  a.start);
    verifyEqual(t, b.count,  a.count);
    if strcmp(a.type, 'surf')
        verifyEqual(t, b.vertlist, a.vertlist, 'Surface vertlist changed.');
    else
        verifyEqual(t, b.voxlist, a.voxlist, 'Voxel IJK list changed.');
    end
end
verifyEqual(t, c2.diminfo{1}.vol.sform, cii.diminfo{1}.vol.sform, 'Subcortical affine changed.');
verifyEqual(t, {c2.diminfo{2}.maps.name}, {cii.diminfo{2}.maps.name}, 'Scalar map names changed.');
end


% -------------------------------------------------------------------------
function test_cifti_dlabel_roundtrip(t)
cii = local_make_synthetic_cifti('dlabel');
f = fullfile(t.TestData.tmp, 'syn.dlabel.nii');
canlab_write_cifti(f, cii);
c2 = canlab_read_cifti(f);

verifyEqual(t, c2.intent, 'dlabel');
verifyLessThanOrEqual(t, max(abs(cii.cdata(:) - c2.cdata(:))), 0, 'Label keys changed.');
tbl = c2.diminfo{2}.maps(1).table;
verifyEqual(t, numel(tbl), 3, 'Label table size changed.');
verifyEqual(t, tbl(2).name, 'RegionA', 'Label name did not round-trip.');
verifyEqual(t, tbl(2).rgba, [1 0 0 1], 'Label RGBA did not round-trip.');
end


% -------------------------------------------------------------------------
function test_real_cifti_if_available(t)
% Opportunistic: if a real dscalar is on the path (Neuroimaging_Pattern_Masks),
% verify an end-to-end native read->write->read with exact cdata.
f0 = which('transcriptomic_gradients.dscalar.nii');
if isempty(f0)
    t.assumeFail('No real .dscalar.nii on path; skipping (synthetic tests cover the codec).');
end
c = canlab_read_cifti(f0);
f = fullfile(t.TestData.tmp, 'real_rt.dscalar.nii');
canlab_write_cifti(f, c);
c2 = canlab_read_cifti(f);
verifyLessThanOrEqual(t, max(abs(double(c.cdata(:)) - double(c2.cdata(:)))), 0, ...
    'Real dscalar cdata did not round-trip exactly.');
verifyEqual(t, numel(c2.diminfo{1}.models), numel(c.diminfo{1}.models));
end


% =========================================================================
function cii = local_make_synthetic_cifti(kind)
% Build a minimal valid grayordinate CIFTI struct: 2 surface models (L/R, 10
% vertices each, 5 in-data) + 1 subcortical voxel model (4 voxels).

mL = struct('struct','CORTEX_LEFT','type','surf','start',1,'count',5, ...
    'numvert',10,'vertlist',[0 2 4 6 8],'voxlist',[]);
mR = struct('struct','CORTEX_RIGHT','type','surf','start',6,'count',5, ...
    'numvert',10,'vertlist',[1 3 5 7 9],'voxlist',[]);
voxijk = [2 2 2; 3 2 2; 2 3 2; 2 2 3]';     % 3x4 IJK (0-based)
mV = struct('struct','THALAMUS_LEFT','type','vox','start',11,'count',4, ...
    'numvert',NaN,'vertlist',[],'voxlist',voxijk);

sform = [2 0 0 -20; 0 2 0 -22; 0 0 2 -24; 0 0 0 1];
dense = struct('type','dense','length',14,'models',{{mL,mR,mV}}, ...
    'vol',struct('dims',[10 10 10],'sform',sform));

G = 14;
switch kind
    case 'dscalar'
        nMaps = 2;
        cdata = single(reshape(1:G*nMaps, G, nMaps));   % deterministic
        maps = struct('name', {'mapOne','mapTwo'}, 'table', {[],[]});
        mapsdim = struct('type','scalars','length',nMaps,'maps',maps);
        intent = 'dscalar';
    case 'dlabel'
        cdata = single(mod((0:G-1)', 3));               % keys 0,1,2
        tbl = struct('key',{0 1 2}, 'name',{'Unknown','RegionA','RegionB'}, ...
            'rgba',{[0 0 0 0],[1 0 0 1],[0 0 1 1]});
        maps = struct('name', {'parc'}, 'table', {tbl});
        mapsdim = struct('type','labels','length',1,'maps',maps);
        intent = 'dlabel';
end

cii = struct('cdata', cdata, 'intent', intent, ...
    'diminfo', {{dense, mapsdim}});
end
