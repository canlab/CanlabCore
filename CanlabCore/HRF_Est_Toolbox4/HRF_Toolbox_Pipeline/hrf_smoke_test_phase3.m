function hrf_smoke_test_phase3(prefix)
%HRF_SMOKE_TEST_PHASE3 Verify Phase 3 v0 scaffold (fmri_hrf, statistic_hrf).
%
% hrf_smoke_test_phase3                - run constructor-only checks (no real data)
% hrf_smoke_test_phase3(prefix)        - also exercise make_fmri_stat_hrf
%                                         against a NIfTI prefix on disk
%
% Phase 3 v0 ships:
%   1. @fmri_hrf      < fmri_data
%   2. @statistic_hrf < statistic_image
%   3. make_fmri_stat_hrf - canonical pairing entry point
%
% This test exercises:
%   * empty construction
%   * lift-from-fmri_data construction with HRF metadata
%   * lift-from-statistic_image construction
%   * cat across two synthetic objects
%   * to_statistic_hrf round-trip
%   * make_fmri_stat_hrf from a real NIfTI prefix (if provided)

fprintf('\n=== Phase 3 v0 smoke test ===\n');

%% Empty construction ---------------------------------------------------
fprintf('\n[1/6] Empty constructors ...\n');
Hb_empty = fmri_hrf();
Ht_empty = statistic_hrf();
assert(isa(Hb_empty, 'fmri_hrf') && isa(Hb_empty, 'fmri_data'), ...
    'fmri_hrf() should produce an fmri_hrf that is-a fmri_data');
assert(isa(Ht_empty, 'statistic_hrf') && isa(Ht_empty, 'statistic_image'), ...
    'statistic_hrf() should produce a statistic_hrf that is-a statistic_image');
fprintf('    ok\n');

%% Lift from synthetic fmri_data ----------------------------------------
fprintf('\n[2/6] Lift from synthetic fmri_data + HRF metadata ...\n');
% Construct a tiny 4-voxel x 6-volume fmri_data fixture (2 conditions x 3 lags).
n_vox = 4; n_vol = 6;
fake_dat = randn(n_vox, n_vol);
fd = local_make_fake_fmri_data(fake_dat);

meta = table( ...
    (1:n_vol)', ...
    string({'pain','pain','pain','neutral','neutral','neutral'}'), ...
    [1;2;3;1;2;3], ...
    [0;0.8;1.6;0;0.8;1.6], ...
    string({'pain_lag1','pain_lag2','pain_lag3','neutral_lag1','neutral_lag2','neutral_lag3'}'), ...
    'VariableNames', {'volume_index','condition','lag_index','lag_seconds','image_label'});

Hb = fmri_hrf(fd, ...
    'MetadataTable', meta, ...
    'Subject', 'sub-99', ...
    'RunLabel', 'task-fake_run-01', ...
    'ModelName', 'sfir', ...
    'TR', 0.8);

assert(isa(Hb, 'fmri_hrf'));
assert(strcmp(Hb.subject, 'sub-99'));
assert(strcmp(Hb.model_name, 'sfir'));
assert(isequal(sort(Hb.conditions), sort({'pain','neutral'})), ...
    'conditions should auto-derive from metadata when not provided');
assert(height(Hb.metadata_table) == size(Hb.dat, 2));
fprintf('    ok\n');

%% Lift from synthetic statistic_image ----------------------------------
fprintf('\n[3/6] Lift from synthetic statistic_image ...\n');
si = statistic_image(fd);
si.type = 'T';
Ht = statistic_hrf(si, ...
    'MetadataTable', meta, ...
    'Subject', 'sub-99', ...
    'RunLabel', 'task-fake_run-01', ...
    'ModelName', 'sfir', ...
    'TR', 0.8);
assert(isa(Ht, 'statistic_hrf'));
assert(strcmp(Ht.subject, 'sub-99'));
fprintf('    ok\n');

%% cat across two synthetic objects -------------------------------------
fprintf('\n[4/6] cat across (subject, run) ...\n');
Hb2 = fmri_hrf(fd, ...
    'MetadataTable', meta, ...
    'Subject', 'sub-98', ...
    'RunLabel', 'task-fake_run-01', ...
    'ModelName', 'sfir', ...
    'TR', 0.8);
Hb_stacked = cat(1, Hb, Hb2);
assert(isa(Hb_stacked, 'fmri_hrf'));
assert(size(Hb_stacked.dat, 2) == 2 * n_vol, ...
    'stacked volume count should be sum of inputs');
assert(height(Hb_stacked.metadata_table) == 2 * n_vol);
assert(any(strcmp('subject', Hb_stacked.metadata_table.Properties.VariableNames)), ...
    'stacked metadata should carry subject column');
fprintf('    ok\n');

%% to_statistic_hrf round-trip -----------------------------------------
fprintf('\n[5/6] to_statistic_hrf from beta + synthetic SE ...\n');
se_si = statistic_image(fd);
se_si.dat = abs(fd.dat) + 0.1;  % strictly positive SE
se_si.type = 'HRF beta standard error';
Ht_derived = to_statistic_hrf(Hb, 'SE', se_si);
assert(isa(Ht_derived, 'statistic_hrf'));
assert(strcmp(Ht_derived.subject, Hb.subject), 'metadata should match parent fmri_hrf');
assert(isequal(size(Ht_derived.dat), size(Hb.dat)));
fprintf('    ok\n');

%% make_fmri_stat_hrf from NIfTI prefix (real data) ---------------------
fprintf('\n[6/6] make_fmri_stat_hrf from a NIfTI prefix ...\n');
if nargin < 1 || isempty(prefix)
    fprintf('    (skipped: pass a prefix to exercise this path)\n');
    fprintf('\n=== Phase 3 v0 smoke test complete (synthetic-only) ===\n');
    return
end

assert(exist([prefix '_beta.nii'], 'file') == 2, ...
    'expected %s_beta.nii to exist', prefix);

[Hb_real, Ht_real] = make_fmri_stat_hrf(prefix, ...
    'Subject', 'sub-real', ...
    'RunLabel', 'task-real_run-01', ...
    'ModelName', 'sfir', ...
    'TR', 0.8);

assert(isa(Hb_real, 'fmri_hrf'));
assert(isa(Ht_real, 'statistic_hrf') || isempty(Ht_real));
assert(~isempty(Hb_real.dat), 'real Hb should have voxel data');
assert(height(Hb_real.metadata_table) == size(Hb_real.dat, 2), ...
    'metadata rows must match volumes');
fprintf('    Hb: %d voxels x %d volumes; metadata %d rows\n', ...
    size(Hb_real.dat, 1), size(Hb_real.dat, 2), height(Hb_real.metadata_table));
if ~isempty(Ht_real.dat)
    fprintf('    Ht: %d voxels x %d volumes; metadata %d rows\n', ...
        size(Ht_real.dat, 1), size(Ht_real.dat, 2), height(Ht_real.metadata_table));
end

fprintf('\n=== Phase 3 v0 smoke test complete ===\n');
end


function fd = local_make_fake_fmri_data(dat)
% Build a minimal fmri_data fixture that the fmri_hrf constructor will accept.
% We use the canlabCore "noverbose" pattern via direct construction from data.
fd = fmri_data();
fd.dat = dat;
fd.removed_voxels = false(size(dat, 1), 1);
fd.removed_images = false(size(dat, 2), 1);
end
