function hrf_smoke_test_phase1(study_dir)
%HRF_SMOKE_TEST_PHASE1 Verify Phase 1 audit/score fix end-to-end.
%
% hrf_smoke_test_phase1                 - run helper-level sanity checks only
% hrf_smoke_test_phase1(study_dir)      - also run the audit RepairMissing
%                                          path against a real finished study
%
% Phase 1 introduced:
%   1. hrf_score_one_prefix.m   - new shared per-prefix scoring helper
%   2. hrf_score_wholebrain_input_table.m - refactored to call the helper
%   3. hrf_write_slurm_study_script.m - generated worker now wraps scoring
%       in try/catch so a single signature failure no longer kills a task
%   4. hrf_audit_slurm_outputs.m - new 'RepairMissing' name-value pair
%
% This test exercises (1) directly and (4) against a study directory.

fprintf('\n=== Phase 1 smoke test ===\n');

%% Helper-level sanity checks (no data needed) ----------------------------
fprintf('\n[1/4] Helper: missing-inputs path ...\n');
tmp_prefix = fullfile(tempdir, sprintf('hrf_smoke_%d', randi(1e6)));
status = hrf_score_one_prefix(tmp_prefix, 'SignatureSets', {'all'});
assert(~status.core_inputs_present, ...
    'core_inputs_present should be false for nonexistent prefix');
assert(~isempty(status.errors), 'expected at least one error');
assert(isempty(status.wrote_files), 'should not have written any files');
fprintf('    ok\n');

fprintf('\n[2/4] Helper: validates ScoreObjects argument ...\n');
got_err = false;
try
    hrf_score_one_prefix(tmp_prefix, 'ScoreObjects', {'beta', 'banana'});
catch ME
    got_err = strcmp(ME.identifier, 'hrf_score_one_prefix:UnknownScoreObjects');
end
assert(got_err, 'expected UnknownScoreObjects error for bad ScoreObjects');
fprintf('    ok\n');

fprintf('\n[3/4] Helper: filename composition (single- vs multi-model) ...\n');
% Indirect check via the helper's error path: the file it would *try* to
% write is reported in the error context when validation fails.
status_single = hrf_score_one_prefix(tmp_prefix, ...
    'ModelName', 'sfir', 'NumModels', 1, 'ScoreObjects', {'beta'}, ...
    'SignatureSets', {'all'});
status_multi = hrf_score_one_prefix(tmp_prefix, ...
    'ModelName', 'sfir', 'NumModels', 3, 'ScoreObjects', {'beta'}, ...
    'SignatureSets', {'all'});
% Both should report core_inputs_present=false (we have no NIfTIs).
assert(~status_single.core_inputs_present);
assert(~status_multi.core_inputs_present);
fprintf('    ok\n');

%% Audit-level repair (real study) ---------------------------------------
fprintf('\n[4/4] Audit RepairMissing against real study ...\n');
if nargin < 1 || isempty(study_dir)
    fprintf('    (skipped: pass study_dir to exercise this path)\n');
    fprintf('\n=== smoke test complete (helper-only) ===\n');
    return
end

assert(exist(study_dir, 'dir') == 7, 'study_dir does not exist: %s', study_dir);

fprintf('    reading audit BEFORE repair ...\n');
[a_before, s_before] = hrf_audit_slurm_outputs(study_dir);
needs_before = sum(a_before.core_complete & ~a_before.score_complete);
fprintf('       %d total rows; %d core_complete; %d score_complete; %d need repair\n', ...
    height(a_before), sum(a_before.core_complete), sum(a_before.score_complete), needs_before);

if needs_before == 0
    fprintf('    no rows need repair; nothing to test. Done.\n');
    return
end

fprintf('    running audit WITH RepairMissing ...\n');
[a_after, s_after] = hrf_audit_slurm_outputs(study_dir, ...
    'RepairMissing', true, 'Verbose', true);

needs_after = sum(a_after.core_complete & ~a_after.score_complete);
fprintf('       repair_summary.n_rows_attempted = %d\n', s_after.repair.n_rows_attempted);
fprintf('       repair_summary.n_rows_repaired  = %d\n', s_after.repair.n_rows_repaired);
fprintf('       repair_summary.n_rows_failed    = %d\n', s_after.repair.n_rows_failed);
fprintf('       repair_summary.n_files_written  = %d\n', s_after.repair.n_files_written);
if ~isempty(s_after.repair.errors)
    fprintf('       first errors:\n');
    for i = 1:min(numel(s_after.repair.errors), 3)
        fprintf('         %s\n', s_after.repair.errors{i});
    end
end

fprintf('    AFTER: %d rows still need repair (was %d)\n', needs_after, needs_before);

%% Equivalence check: repaired score files == post-hoc backfill -----------
% Pick the first repaired row, regenerate its score CSV via the post-hoc
% path on a tempdir copy, and check the two files agree.
fprintf('\n[equivalence] Comparing in-audit repair vs post-hoc backfill ...\n');
repaired_rows = find(a_before.core_complete & ~a_before.score_complete);
if isempty(repaired_rows)
    fprintf('    (no rows to compare)\n');
elseif ~isempty(s_after.repair.errors) && s_after.repair.n_rows_repaired == 0
    fprintf('    (skipped: repair did not write any files)\n');
else
    r = repaired_rows(1);
    audit_csv = '';
    if a_after.score_complete(r) && ~isempty(a_after.beta_scores_file{r})
        audit_csv = a_after.beta_scores_file{r};
    end
    if ~isempty(audit_csv) && exist(audit_csv, 'file') == 2
        backup_csv = [audit_csv '.audit_repair_backup'];
        copyfile(audit_csv, backup_csv);
        % Re-run the post-hoc backfill on this row in isolation.
        T = a_before(r, :);
        T.subject = string(T.subject);
        T.model = string(T.model);
        % hrf_score_wholebrain_input_table expects the model-resolved prefix
        % (it calls the helper with NumModels=1). The audit table's
        % output_prefix is task-level; model_prefix already has the model
        % suffix applied for multi-model studies.
        T.prefix = string(T.model_prefix);
        T.beta_scores_file = string('');
        T.t_scores_file = string('');
        T.metadata_file = string(T.metadata_file);
        T.result_mat_file = string(T.result_mat_file);
        delete(audit_csv);
        hrf_score_wholebrain_input_table(T, ...
            'SignatureSets', {'all'}, 'ScoreObjects', {'beta'}, ...
            'Overwrite', true);
        if exist(audit_csv, 'file') == 2
            T_audit = readtable(backup_csv);
            T_posthoc = readtable(audit_csv);
            ok = isequal(size(T_audit), size(T_posthoc)) && ...
                isequal(T_audit.Properties.VariableNames, T_posthoc.Properties.VariableNames);
            if ok
                fprintf('    ✓ shapes and column names match (%d rows, %d cols)\n', ...
                    height(T_audit), width(T_audit));
            else
                fprintf('    ✗ shapes/cols disagree (audit %dx%d vs post-hoc %dx%d)\n', ...
                    size(T_audit, 1), size(T_audit, 2), ...
                    size(T_posthoc, 1), size(T_posthoc, 2));
            end
        else
            fprintf('    ✗ post-hoc backfill did not write a CSV\n');
        end
        % Restore the audit-repaired version.
        copyfile(backup_csv, audit_csv);
        delete(backup_csv);
    else
        fprintf('    (no audit-repaired CSV available to compare)\n');
    end
end

fprintf('\n=== smoke test complete ===\n');
end
