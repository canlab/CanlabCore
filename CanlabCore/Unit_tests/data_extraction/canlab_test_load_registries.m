function tests = canlab_test_load_registries
% Unit tests for the load_image_set / load_atlas keyword registries:
% newly-added surface (CIFTI) and volumetric map sets, the categorized 'list'
% output, and the Tian subcortical scale atlases.
%
% Requires the full Neuroimaging_Pattern_Masks data (large CIFTI/.mat binaries) on
% the path. Individual tests assumeFail (skip) when a specific data file is not
% present, so this file can never *fail* CI -- but it is intentionally SKIPPED
% entirely under GitHub Actions (see setupOnce) so the CI unit tier does not depend
% on those large data files. It runs normally when executed locally.
%
% Run:    runtests('canlab_test_load_registries')
%
% :See also: load_image_set, load_atlas

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
function setupOnce(t)
% Skip the whole file in GitHub Actions CI (kept out of the CI set by request);
% still runs locally. GITHUB_ACTIONS is set to 'true' on all GitHub runners.
if ~isempty(getenv('GITHUB_ACTIONS'))
    t.assumeFail(['Skipped in GitHub Actions: this registry test needs the full ' ...
        'Neuroimaging_Pattern_Masks data files. It runs locally.']);
end
assert(~isempty(which('load_image_set')), 'load_image_set not on path.');
assert(~isempty(which('load_atlas')),     'load_atlas not on path.');
t.TestData.figvis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
end

function teardownOnce(t)
close all force
% figvis is unset if setupOnce bailed early (e.g. the GitHub Actions skip), so
% guard against it -- an error here would flip filtered tests to "failed".
if isfield(t.TestData, 'figvis')
    set(0, 'DefaultFigureVisible', t.TestData.figvis);
end
end


% -------------------------------------------------------------------------
function test_hcp_groupica_surface(t)
% hcp_ica15/25/50 return fmri_surface_data with the right number of maps.
specs = {'hcp_ica15', 15; 'hcp_ica25', 25; 'hcp_ica50', 50};
for i = 1:size(specs, 1)
    if isempty(which(sprintf('hcp_d%d_ICs.dscalar.nii', specs{i,2})))
        t.assumeFail('HCP group-ICA dscalar not on path.');
    end
    o = load_image_set(specs{i,1}, 'noverbose');
    verifyTrue(t, isa(o, 'fmri_surface_data'), 'HCP ICA must return fmri_surface_data.');
    verifyEqual(t, size(o.dat, 2), specs{i,2}, 'Wrong number of ICA components.');
end
end


% -------------------------------------------------------------------------
function test_spectral_bases_surface(t)
if isempty(which('spectral_bases_200.dscalar.nii'))
    t.assumeFail('spectral_bases_200.dscalar.nii not on path.');
end
o = load_image_set('spectral_bases', 'noverbose');
verifyTrue(t, isa(o, 'fmri_surface_data'), 'spectral_bases must return fmri_surface_data.');
verifyEqual(t, size(o.dat, 2), 200, 'Expected 200 spectral basis functions.');
end


% -------------------------------------------------------------------------
function test_mito_maps_volumetric(t)
if isempty(which('mito_maps.mat'))
    t.assumeFail('mito_maps.mat not on path.');
end
[o, names] = load_image_set('mito_maps', 'noverbose');
verifyTrue(t, isa(o, 'fmri_data'), 'mito_maps must return fmri_data.');
verifyEqual(t, size(o.dat, 2), 6, 'Expected 6 mitochondrial maps.');
verifyEqual(t, numel(names), 6, 'Expected 6 map names.');
end


% -------------------------------------------------------------------------
function test_list_returns_struct_of_tables(t)
evalc('tmp = load_image_set(''list'');');  % capture printed tables to suppress them
verifyClass(t, tmp, 'struct');
expected = {'signatures','datasets','networks','surface','gradients','metaanalysis'};
for f = expected
    verifyTrue(t, isfield(tmp, f{1}), sprintf('list struct must have field %s.', f{1}));
    verifyClass(t, tmp.(f{1}), 'table');
end
% Surface table must advertise the new fmri_surface_data sets.
verifyTrue(t, any(contains(tmp.surface.keyword, 'hcp_ica15')), 'Surface table lists hcp_ica15.');
verifyTrue(t, any(contains(tmp.surface.keyword, 'spectral_bases')), 'Surface table lists spectral_bases.');
end


% -------------------------------------------------------------------------
function test_tian_scale_atlases(t)
% Tian S1-S4 load with the expected subcortical region counts.
if isempty(which('tian_3t_s1_fmriprep20_atlas_object.mat'))
    t.assumeFail('Tian scale atlas objects not on path.');
end
expected = [16 32 50 54];
for s = 1:4
    a = load_atlas(sprintf('tian_s%d', s), 'noverbose');
    verifyTrue(t, isa(a, 'atlas'), 'Tian scale must return an atlas.');
    verifyEqual(t, numel(a.labels), expected(s), ...
        sprintf('Tian S%d should have %d regions.', s, expected(s)));
end
% fsl6 variant also resolves.
a = load_atlas('tian_s2_fsl6', 'noverbose');
verifyEqual(t, numel(a.labels), 32, 'Tian S2 fsl6 should have 32 regions.');
end


% -------------------------------------------------------------------------
function test_load_atlas_list_returns_struct(t)
evalc('la = load_atlas(''list'');');   % capture printed tables to suppress them
verifyClass(t, la, 'struct');
for f = {'combined','cortex','subcortical','thalamus_hypothalamus','brainstem_cerebellum','networks_specialized'}
    verifyTrue(t, isfield(la, f{1}), sprintf('atlas list must have field %s.', f{1}));
    verifyClass(t, la.(f{1}), 'table');
end
end
