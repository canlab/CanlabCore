%% RSM reliability (ICC) across sessions
% Reproduces the 01282025 RSM Reliability workflow: per-subject ICC of the
% RSM across that subject's sessions, aggregated across subjects, at several
% grouping levels.

dat = make_synthetic_rsa_data('n_sub', 9, 'n_ses', 5);

%% Per-session RSM (one slice per subject x session)
R_sess = compute_rsm(dat, 'group_by', {'condition','bodysite'}, ...
                          'level', 'session', 'subject_var', 'sub', ...
                          'session_var', 'sesno', 'metric', 'spearman');

%% Whole-RSM reliability: per-subject ICC(3,k), then aggregate
out = R_sess.reliability('icc_type', '3-k');
fprintf('Whole-RSM ICC across %d subjects: mean=%.3f, median=%.3f, std=%.3f\n', ...
    out.summary.n_subjects, out.summary.mean, out.summary.median, out.summary.std);
disp(out.per_subject)

%% Per-grouping reliability (per superordinate condition + per bodysite)
T_rel_g = R_sess.reliability_by_grouping('icc_type', '3-k');
disp(T_rel_g)
% NOTE: condition-level rows (hot/warm/imagine, 28 cells) are reliable;
% bodysite-level rows (3 cells) are flagged as statistically unreliable.

%% Per-individual-condition reliability (row-level)
T_rel_pc = R_sess.reliability_per_condition('icc_type', '3-k');
disp(head(T_rel_pc, 8))

%% Bar plot of subject reliabilities
figure;
barplot_columns(out.per_subject.ICC', 'names', out.per_subject.sub', ...
    'noviolin', 'nostars', 'noind');
ylim([0 1.1]); ylabel('Reliability (ICC 3-k)'); title('RSM reliability by subject');
