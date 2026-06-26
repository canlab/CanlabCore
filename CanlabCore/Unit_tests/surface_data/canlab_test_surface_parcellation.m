function tests = canlab_test_surface_parcellation
% Unit tests for fmri_surface_data parcellation / region methods (M7).
%
% Covers apply_parcellation (parcel means, labels, area weighting, space check),
% cluster-extent threshold ('k'), and surface_region. Core tests use a synthetic
% object with a known parcellation (self-contained); mesh-dependent tests use a
% real fs_LR file if present.
%
% Run:    runtests('canlab_test_surface_parcellation')
%
% :See also: apply_parcellation, threshold, surface_region, reparse_contiguous

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'Surface_tools'));
addpath(fullfile(here, '..', '..'));
assert(~isempty(which('fmri_surface_data')), 'fmri_surface_data not on path.');
t.TestData.figvis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
end

function teardownOnce(t)
close all force
set(0, 'DefaultFigureVisible', t.TestData.figvis);
end


% -------------------------------------------------------------------------
function test_apply_parcellation_known(t)
% 20 grayordinates, 3 parcels + background (0). Data == key, so each parcel mean
% equals its key.
keys = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 0 0 0 0 0]';
D = [double(keys), double(keys)*10];
o = local_obj(D);

[pm, labels, tbl] = apply_parcellation(o, keys);
verifyEqual(t, size(pm), [2 3], 'parcel_means should be [nMaps x nParcels].');
verifyEqual(t, pm(1,:), [1 2 3], 'AbsTol', 1e-6, 'map1 parcel means should equal the keys.');
verifyEqual(t, pm(2,:), [10 20 30], 'AbsTol', 1e-5, 'map2 parcel means should equal 10*keys.');
verifyEqual(t, numel(labels), 3);
verifyEqual(t, tbl.n_grayordinates', [5 5 5], 'each parcel has 5 grayordinates.');
% background key 0 must be excluded
verifyFalse(t, any(tbl.key == 0), 'background (key 0) must be excluded.');
end


% -------------------------------------------------------------------------
function test_apply_parcellation_from_object_and_space_check(t)
keys = [ones(10,1); 2*ones(10,1)];
o = local_obj([double(keys), double(keys)]);
parc = local_obj(double(keys));        % a "dlabel"-like object on the same space
parc.intent = 'dlabel';
parc.label_table = struct('key', {1 2}, 'name', {'A','B'}, 'rgba', {[1 0 0 1],[0 0 1 1]});

[pm, labels] = apply_parcellation(o, parc);
verifyEqual(t, pm(1,:), [1 2], 'AbsTol', 1e-6);
verifyEqual(t, labels, {'A','B'}, 'labels should come from the label_table.');

% Different space must error
parc2 = parc; parc2.surface_space = 'fsaverage_164k';
verifyError(t, @() apply_parcellation(o, parc2), 'fmri_surface_data:apply_parcellation:space');
end


% -------------------------------------------------------------------------
function test_real_atlas_parcellation(t)
f = which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii');
if isempty(f), t.assumeFail('No real .dlabel atlas on path.'); end
atl = fmri_surface_data(f);
% Use the atlas keys as the data: each parcel mean must equal its key
data = atl; data.dat = single(double(atl.dat)); data.intent = 'dscalar';
data.removed_images = false(1,1);
[pm, labels] = apply_parcellation(data, atl);
ukeys = unique(round(double(atl.dat(atl.dat > 0))));
verifyEqual(t, numel(labels), numel(ukeys), 'one label per positive key.');
verifyEqual(t, pm(1,:)', ukeys, 'AbsTol', 1e-4, 'parcel mean of key-data must equal the key.');
% area-weighted means are finite
pmA = apply_parcellation(data, atl, 'area');
verifyTrue(t, all(isfinite(pmA(:))), 'area-weighted parcel means must be finite.');
end


% -------------------------------------------------------------------------
function test_cluster_extent_threshold(t)
f = which('transcriptomic_gradients.dscalar.nii');
if isempty(f), t.assumeFail('No fs_LR dscalar on path.'); end
s = fmri_surface_data(f);
raw = threshold(s, 1.0, 'positive');
big = threshold(s, 1.0, 'positive', 'k', 50);
verifyLessThanOrEqual(t, nnz(big.dat(:,1) ~= 0), nnz(raw.dat(:,1) ~= 0), ...
    'cluster-extent threshold cannot increase the surviving set.');
% Every surviving cluster must have >= 50 grayordinates
[chk, ~] = reparse_contiguous(big, 'which_image', 1);
cl = chk.brain_model.cluster;
if any(cl > 0)
    szs = accumarray(cl(cl > 0), 1);
    verifyGreaterThanOrEqual(t, min(szs(szs > 0)), 50, 'a surviving cluster is smaller than k.');
end
end


% -------------------------------------------------------------------------
function test_surface_region(t)
f = which('transcriptomic_gradients.dscalar.nii');
if isempty(f), t.assumeFail('No fs_LR dscalar on path.'); end
s = threshold(fmri_surface_data(f), 1.5, 'positive');
reg = surface_region(s, 'which_image', 1);
verifyGreaterThan(t, numel(reg), 0, 'expected at least one region.');
% Field sanity on the first region
r1 = reg(1);
verifyTrue(t, all(isfield(r1, {'struct','type','grayord_rows','XYZmm','numVox','val'})));
verifyEqual(t, r1.numVox, numel(r1.grayord_rows));
% Total grayordinates across regions == active grayordinates
total = sum([reg.numVox]);
active = nnz(s.dat(:,1) ~= 0 & ~isnan(s.dat(:,1)));
verifyEqual(t, total, active, 'regions must partition the active grayordinates.');
end


% =========================================================================
function o = local_obj(D)
n = size(D,1); h = n/2;
mL = struct('struct','CORTEX_LEFT','type','surf','start',1,'count',h,'numvert',32492,'vertlist',0:h-1,'voxlist',[]);
mR = struct('struct','CORTEX_RIGHT','type','surf','start',h+1,'count',h,'numvert',32492,'vertlist',0:h-1,'voxlist',[]);
bm = struct('type','dense','length',n,'models',{{mL,mR}},'vol',[]);
bm.grayordinate_type = 'cortex_only'; bm.cluster = [];
o = fmri_surface_data('dat', single(D), 'brain_model', bm, ...
    'surface_space', 'fsLR_32k', 'intent', 'dscalar');
end
