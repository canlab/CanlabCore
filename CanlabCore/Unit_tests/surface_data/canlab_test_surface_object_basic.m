function tests = canlab_test_surface_object_basic
% Unit tests for the fmri_surface_data class skeleton (M2).
%
% Verifies construction, the grayordinate property model, inheritance from
% image_vector, and the D5b "no empty-squeezing" behavior -- all with NO
% external toolbox. The core tests build a synthetic grayordinate object so they
% run anywhere; one opportunistic test uses a real .dscalar.nii if present.
%
% Run:    runtests('canlab_test_surface_object_basic')
%
% :See also: fmri_surface_data, canlab_read_cifti, canlab_test_surface_io

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'Surface_tools'));
addpath(fullfile(here, '..', '..'));        % so @fmri_surface_data resolves
assert(~isempty(which('fmri_surface_data')), 'fmri_surface_data not on path.');
t.TestData.cii = local_synthetic_cifti();
end


% -------------------------------------------------------------------------
function test_empty_construction(t)
o = fmri_surface_data;
verifyTrue(t, isa(o, 'fmri_surface_data'));
verifyTrue(t, isa(o, 'image_vector'), 'Must subclass image_vector.');
verifyEmpty(t, o.dat);
verifyEmpty(t, o.brain_model);
verifyEqual(t, o.intent, '');
end


% -------------------------------------------------------------------------
function test_build_from_grayordinate_struct(t)
o = fmri_surface_data(t.TestData.cii);
verifyEqual(t, class(o), 'fmri_surface_data');
verifyEqual(t, size(o.dat), [14 2], 'Grayordinate dat should be 14 x 2.');
verifyEqual(t, o.intent, 'dscalar');

% brain_model carries the model split, 1:1 with .dat rows
bm = o.brain_model;
verifyEqual(t, numel(bm.models), 3, 'Expected L cortex + R cortex + 1 voxel model.');
total = sum(cellfun(@(m) m.count, bm.models));
verifyEqual(t, total, size(o.dat, 1), 'Sum of model counts must equal grayordinate rows.');
verifyEqual(t, bm.models{1}.struct, 'CORTEX_LEFT');
verifyEqual(t, bm.models{3}.type, 'vox');
verifyEqual(t, bm.vol.sform, t.TestData.cii.diminfo{1}.vol.sform, 'Subcortical affine lost.');

% map names carried
verifyEqual(t, reshape(o.image_names,1,[]), {'mapOne','mapTwo'});
end


% -------------------------------------------------------------------------
function test_removed_vectors_are_vestigial(t)
% D5b: removed_voxels/removed_images are all-false and length-correct.
o = fmri_surface_data(t.TestData.cii);
verifyEqual(t, numel(o.removed_voxels), size(o.dat,1));
verifyEqual(t, numel(o.removed_images), size(o.dat,2));
verifyFalse(t, any(o.removed_voxels), 'removed_voxels must be all-false.');
verifyFalse(t, any(o.removed_images), 'removed_images must be all-false.');
end


% -------------------------------------------------------------------------
function test_remove_replace_empty_are_noops(t)
% Even with all-zero grayordinate rows, nothing is squeezed.
o = fmri_surface_data(t.TestData.cii);
o.dat(7, :) = 0;        % zero out a whole grayordinate row
n0 = size(o.dat, 1);

o1 = remove_empty(o);
verifyEqual(t, size(o1.dat, 1), n0, 'remove_empty must NOT drop rows.');
verifyEqual(t, o1.dat, o.dat, 'remove_empty must return data unchanged.');
verifyFalse(t, any(o1.removed_voxels), 'removed_voxels must stay all-false.');

o2 = replace_empty(o1);
verifyEqual(t, size(o2.dat, 1), n0, 'replace_empty must not change size.');
verifyEqual(t, o2.dat, o.dat, 'replace_empty must return data unchanged.');
end


% -------------------------------------------------------------------------
function test_get_wh_image_subsets_maps_only(t)
% Inherited get_wh_image selects map columns; grayordinate rows are preserved.
o = fmri_surface_data(t.TestData.cii);
o2 = get_wh_image(o, 2);
verifyEqual(t, size(o2.dat), [14 1], 'get_wh_image should keep all rows, 1 map.');
verifyEqual(t, o2.dat, o.dat(:, 2), 'Wrong map selected.');
verifyEqual(t, numel(o2.image_names), 1, 'image_names not subset with maps.');
verifyEqual(t, numel(o2.removed_images), 1, 'removed_images not subset with maps.');
% rows (grayordinates) unchanged
verifyEqual(t, numel(o2.removed_voxels), 14, 'Grayordinate rows must be preserved.');
end


% -------------------------------------------------------------------------
function test_method_surface_is_inherited(t)
% The image_vector method surface is inherited by dispatch (same method names
% as fmri_data/image_vector). Note: methods that dereference volInfo
% (descriptives, montage, orthviews, flip, ...) get surface-aware guards/
% overrides in M3 (Risk #1); here we just confirm the surface is present.
o = fmri_surface_data(t.TestData.cii);
m = methods(o);
for name = {'get_wh_image','remove_empty','replace_empty','apply_parcellation','ica'}
    verifyTrue(t, ismember(name{1}, m), sprintf('Method %s should be on the object.', name{1}));
end
% remove_empty / replace_empty must be overridden in THIS class folder (no-op
% semantics are verified separately in test_remove_replace_empty_are_noops).
clsdir = fileparts(which('fmri_surface_data'));
verifyTrue(t, exist(fullfile(clsdir, 'remove_empty.m'), 'file') == 2, ...
    'remove_empty.m override missing from @fmri_surface_data.');
verifyTrue(t, exist(fullfile(clsdir, 'replace_empty.m'), 'file') == 2, ...
    'replace_empty.m override missing from @fmri_surface_data.');
end


% -------------------------------------------------------------------------
function test_real_dscalar_if_available(t)
f = which('transcriptomic_gradients.dscalar.nii');
if isempty(f)
    t.assumeFail('No real .dscalar.nii on path; synthetic tests cover the class.');
end
o = fmri_surface_data(f);
verifyEqual(t, class(o), 'fmri_surface_data');
verifyGreaterThan(t, size(o.dat,1), 60000, 'Expected tens of thousands of grayordinates.');
verifyEqual(t, sum(cellfun(@(m) m.count, o.brain_model.models)), size(o.dat,1));
verifyEqual(t, o.brain_model.models{1}.struct, 'CORTEX_LEFT');
% .dat is full, no squeezing
verifyEqual(t, numel(o.removed_voxels), size(o.dat,1));
verifyFalse(t, any(o.removed_voxels));
end


% =========================================================================
function cii = local_synthetic_cifti()
% Minimal grayordinate dscalar: 2 surface models (L/R, 5 in-data each) + 1
% subcortical voxel model (4 voxels) = 14 grayordinates, 2 maps.
mL = struct('struct','CORTEX_LEFT','type','surf','start',1,'count',5, ...
    'numvert',32492,'vertlist',[0 2 4 6 8],'voxlist',[]);
mR = struct('struct','CORTEX_RIGHT','type','surf','start',6,'count',5, ...
    'numvert',32492,'vertlist',[1 3 5 7 9],'voxlist',[]);
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
