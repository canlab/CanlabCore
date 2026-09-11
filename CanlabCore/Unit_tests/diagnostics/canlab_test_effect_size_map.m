function tests = canlab_test_effect_size_map
%CANLAB_TEST_EFFECT_SIZE_MAP Tests for canlab_effect_size_map and its helpers.
%
%   Uses the 30-image emotionreg sample as stand-in first-level contrast
%   images, with synthetic ResMS images and a synthetic first-level design.
%   Voxels are subsampled so the tests run quickly.

tests = functiontests(localfunctions);
end


% -------------------------------------------------------------------------
% Fixtures
% -------------------------------------------------------------------------

function [con_obj, resms_obj, X, c] = get_synthetic_inputs(n_subjects)
% Contrast images from the sample data (subsampled voxels), synthetic
% ResMS images and design.
con_obj = canlab_get_sample_fmri_data();
con_obj = get_wh_image(con_obj, 1:n_subjects);

% Keep every 100th voxel so power calculations are fast
wh_keep = false(size(con_obj.dat, 1), 1);
wh_keep(1:100:end) = true;
con_obj.dat(~wh_keep, :) = 0;

resms_obj = con_obj;
resms_obj.dat = zeros(size(con_obj.dat), 'like', con_obj.dat);
resms_obj.dat(wh_keep, :) = 20 + 5 * rand(sum(wh_keep), n_subjects);

X = [randn(200, 1) ones(200, 1)];
c = [1 0]';
end


function common = fast_options()
% Options shared by the tests: no plots, no printing, Bonferroni only (no
% SPM smoothness estimation), a short N range, fast power approximation.
common = {'noplot', 'verbose', false, 'correction', 'bonferroni', ...
    'N_range', 5:5:60, 'power_method', 'shift'};
end


% -------------------------------------------------------------------------
% canlab_scan_time_allocation
% -------------------------------------------------------------------------

function test_scan_time_allocation_basic(tc)
alloc = canlab_scan_time_allocation(60, 'verbose', false);

tc.verifyEqual(alloc.hours, 60);
tc.verifyTrue(all(diff(alloc.N) > 0));
tc.verifyTrue(all(diff(alloc.hours_per_subject) < 0), 'Hours per subject must fall as N rises');
tc.verifyTrue(all(alloc.sessions_per_subject == round(alloc.sessions_per_subject)));
tc.verifyTrue(all(alloc.images_per_subject > 0));
tc.verifyTrue(all(alloc.effective_images_per_subject >= alloc.parameters.min_images));
tc.verifyEqual(numel(alloc.N), numel(alloc.functional_min_per_subject));

% 3 subjects sharing 60 hours: 20 hours each, 14 sessions of 1.5 h
tc.verifyEqual(alloc.N(1), 3);
tc.verifyEqual(alloc.sessions_per_subject(1), 14);
end


function test_scan_time_allocation_infeasible_errors(tc)
% Not enough time to run even one subject's first-level model
tc.verifyError(@() canlab_scan_time_allocation(0.1, 'verbose', false), ...
    'canlab_scan_time_allocation:NoFeasibleN');
end


% -------------------------------------------------------------------------
% canlab_effect_size_map: object mode
% -------------------------------------------------------------------------

function test_object_mode_with_resms_and_design(tc)
[con_obj, resms_obj, X, c] = get_synthetic_inputs(30);

[results, maps] = canlab_effect_size_map(con_obj, resms_obj, 'X', X, 'c', c, ...
    'hours', 60, 'TR', 2, fast_options{:});

N = 30;
nvox = results.n_voxels;

tc.verifyEqual(results.N, N);
tc.verifyTrue(nvox > 100 && nvox <= sum(any(con_obj.dat ~= 0, 2)));
tc.verifyTrue(results.have_within_decomposition);

% Effect size and variance components
tc.verifyEqual(size(results.d), [nvox 1]);
tc.verifyEqual(results.d, results.t ./ sqrt(N), 'AbsTol', 1e-10);
tc.verifyTrue(all(results.sig2_between >= 0));
tc.verifyTrue(all(results.sig2_within > 0));
tc.verifyEqual(results.sig2_within + results.sig2_between, ...
    max(results.sig2_total, results.sig2_within), 'RelTol', 1e-8);

% Design variance from X and c
pX = pinv(X);
tc.verifyEqual(results.design.desvar, c' * (pX * pX') * c, 'RelTol', 1e-10);
tc.verifyEqual(results.design.nscan, 200);
tc.verifyEqual(results.design.k, 2);
tc.verifyEqual(results.design.min_images, 2);

% Allocation and power
alloc = results.allocation;
nN = numel(alloc.N);
tc.verifyEqual(size(alloc.expected_t), [nvox nN]);
tc.verifyEqual(size(alloc.power_corrected), [nvox nN]);
tc.verifyEqual(size(alloc.power_uncorrected), [nvox nN 2]);
tc.verifyTrue(all(alloc.power_corrected(:) >= 0 & alloc.power_corrected(:) <= 1));
tc.verifyTrue(all(alloc.power_uncorrected(:) >= 0 & alloc.power_uncorrected(:) <= 1));
tc.verifyTrue(all(alloc.power_uncorrected(:, :, 1) >= alloc.power_corrected - 1e-12, 'all'), ...
    'Uncorrected p < .05 power must be at least the FWE-corrected power');
tc.verifyTrue(ismember(alloc.best.N, alloc.N));
tc.verifyEqual(results.power_map, alloc.power_corrected(:, alloc.best.index));

% Bonferroni: NISC equals the number of voxels
tc.verifyEqual(results.smoothness.method, 'bonferroni');
tc.verifyEqual(results.smoothness.NISC, nvox, 'RelTol', 1e-4);

% Replication power
tc.verifyEqual(size(results.replication.power_corrected), [nvox 1]);
tc.verifyTrue(all(results.replication.power_corrected >= 0 & results.replication.power_corrected <= 1));

% Maps
tc.verifyClass(maps.power_best, 'fmri_data');
tc.verifyEqual(size(maps.power_best.dat, 1), nvox);
tc.verifyEqual(double(maps.cohens_d.dat), results.d, 'RelTol', 1e-5, 'AbsTol', 1e-6);
tc.verifyEqual(double(maps.between_subjects_std.dat), sqrt(results.sig2_between), 'RelTol', 1e-5, 'AbsTol', 1e-6);
end


function test_object_mode_without_within_information(tc)
% Contrast images only: all variance is treated as between-subject
[con_obj] = get_synthetic_inputs(20);

results = canlab_effect_size_map(con_obj, fast_options{:});

tc.verifyEqual(results.N, 20);
tc.verifyFalse(results.have_within_decomposition);
tc.verifyTrue(all(results.sig2_within == 0));
tc.verifyEqual(results.sig2_between, results.sig2_total, 'RelTol', 1e-10);
tc.verifyTrue(all(diff(results.allocation.mean_power_corrected) >= -1e-12), ...
    'With no within-subject variance, power must not decrease with N');
end


function test_noncentral_power_method(tc)
% Exact noncentral t power runs and is close to the shift approximation
[con_obj, resms_obj, X, c] = get_synthetic_inputs(15);

opts = {'noplot', 'verbose', false, 'correction', 'bonferroni', 'N_range', 10:10:40};
r_nc = canlab_effect_size_map(con_obj, resms_obj, 'X', X, 'c', c, opts{:}, 'power_method', 'noncentral');
r_sh = canlab_effect_size_map(con_obj, resms_obj, 'X', X, 'c', c, opts{:}, 'power_method', 'shift');

tc.verifyTrue(all(r_nc.allocation.power_corrected(:) >= 0 & r_nc.allocation.power_corrected(:) <= 1));
tc.verifyEqual(r_nc.allocation.mean_power_corrected, r_sh.allocation.mean_power_corrected, 'AbsTol', 0.15);
end


function test_group_maps_must_be_paired(tc)
[con_obj] = get_synthetic_inputs(10);
tc.verifyError(@() canlab_effect_size_map(con_obj, 'group_con', con_obj, fast_options{:}), ...
    'canlab_effect_size_map:GroupMaps');
end


function test_design_variance_requires_nscan(tc)
[con_obj, resms_obj] = get_synthetic_inputs(10);
tc.verifyError(@() canlab_effect_size_map(con_obj, resms_obj, 'design_variance', 0.01, fast_options{:}), ...
    'canlab_effect_size_map:NoNscan');
end


% -------------------------------------------------------------------------
% canlab_effect_size_map: SPM folder mode
% -------------------------------------------------------------------------

function test_spm_folder_mode(tc)
% Build a fake set of first-level SPM folders and read them back
n_subjects = 6;
[con_obj, resms_obj, X, c] = get_synthetic_inputs(n_subjects);

parent = tempname;
mkdir(parent);
tc.addTeardown(@() rmdir(parent, 's'));

pX = pinv(X);
for i = 1:n_subjects
    d = fullfile(parent, sprintf('sub-%02d', i));
    mkdir(d);

    obj = get_wh_image(con_obj, i);
    obj.fullpath = fullfile(d, 'con_0001.nii');
    evalc('write(obj, ''overwrite'')');

    obj = get_wh_image(resms_obj, i);
    obj.fullpath = fullfile(d, 'ResMS.nii');
    evalc('write(obj, ''overwrite'')');

    SPM = struct();
    SPM.xX.X = X;
    SPM.xX.xKXs.X = X;
    SPM.xX.Bcov = pX * pX';
    SPM.xX.erdf = 170;
    SPM.nscan = 200;
    SPM.xY.RT = 2;
    SPM.xCon = struct('name', {'baseline', 'task'}, 'STAT', {'T', 'T'}, ...
        'c', {[0 1]', c}, 'Vcon', {struct('fname', 'con_0002.nii'), struct('fname', 'con_0001.nii')});
    SPM.VResMS = struct('fname', 'ResMS.nii');
    save(fullfile(d, 'SPM.mat'), 'SPM');
end

% Select the contrast by name; 'task' is contrast #2 whose image is con_0001.nii
results = canlab_effect_size_map(parent, 'contrast', 'task', fast_options{:});

tc.verifyEqual(results.N, n_subjects);
tc.verifyEqual(results.contrast_name, 'task');
tc.verifyEqual(numel(results.subject_names), n_subjects);
tc.verifyTrue(results.have_within_decomposition);
tc.verifyEqual(results.design.TR, 2);
tc.verifyEqual(numel(results.design.desvar_per_subject), n_subjects);
tc.verifyEqual(results.design.desvar_per_subject(1), c' * (pX * pX') * c, 'RelTol', 1e-10);
tc.verifyEqual(results.design.df_shrink_factor, 170 / 200, 'RelTol', 1e-10);
tc.verifyTrue(all(results.sig2_between >= 0));

% Same analysis from the objects directly should give the same effect sizes
r_obj = canlab_effect_size_map(con_obj, resms_obj, 'X', X, 'c', c, fast_options{:});
tc.verifyEqual(results.n_voxels, r_obj.n_voxels);
tc.verifyEqual(results.summary.mean_d, r_obj.summary.mean_d, 'AbsTol', 1e-3);

% A bad contrast name gives a helpful error
tc.verifyError(@() canlab_effect_size_map(parent, 'contrast', 'nonexistent', fast_options{:}), ...
    'canlab_effect_size_map:BadContrast');
end


% -------------------------------------------------------------------------
% canlab_power_allocation_curves
% -------------------------------------------------------------------------

function test_power_allocation_curves(tc)
out = canlab_power_allocation_curves('within_std', [20 40], 'between_std', [2 4 6], ...
    'N_range', 5:5:60, 'doplot', false, 'verbose', false, 'power_method', 'shift');

nN = numel(out.alloc.N);
tc.verifyEqual(size(out.power), [2 nN 3]);
tc.verifyTrue(all(out.power(:) >= 0 & out.power(:) <= 1));
tc.verifyEqual(size(out.optimal_N), [2 3]);
tc.verifyTrue(all(ismember(out.optimal_N(:), out.alloc.N)));

% Less between-subject variance and less within-subject noise -> more power
tc.verifyTrue(all(out.power(:, :, 1) >= out.power(:, :, 3), 'all'));
tc.verifyTrue(all(out.power(1, :, :) >= out.power(2, :, :), 'all'));
end


% -------------------------------------------------------------------------
% Deprecated wrapper
% -------------------------------------------------------------------------

function test_deprecated_wrapper_warns(tc)
tc.verifyWarning(@() effect_size_map(), 'effect_size_map:Deprecated');
end
