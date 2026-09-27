%% RSA quick-start
% Build an RSM, visualize it, run within/between contrasts, and compare to
% model RDMs. Runs end-to-end on synthetic data.
%
% For real data, replace make_synthetic_rsa_data() with your fmri_data
% object and set group_by / subject_var to your metadata column names.

dat = make_synthetic_rsa_data();

%% 1. Build a per-subject RSM (24 conditions = 3 conditions x 8 bodysites)
R = compute_rsm(dat, 'group_by', {'condition','bodysite'}, ...
                     'subject_var', 'sub', 'metric', 'spearman');
disp(R)

%% 2. Group-mean RSM heatmap with condition block borders
figure; plot(mean(R), 'block_borders_by', 'condition');
title('Group-mean RSM');

%% 3. Within vs between condition contrasts (FDR-corrected)
spec = {
    'within_hot',   'hot',                 [];
    'within_warm',  'warm',                [];
    'within_imag',  'imagine',             [];
    'hot_vs_warm',  'hot',                 'warm';
    'HW_vs_HI',     {'hot','warm'},        {'hot','imagine'};
};
T = R.ttest_contrasts(spec, 'tail', 'both', 'correction', 'fdr');
disp(T)

%% 4. Bar plot of within-condition similarities
figure;
plot_rsm_contrast_bars(R.cells_table({'hot','warm','imagine'}), ...
    'colors', {[1 .6 .6],[.6 1 .6],[.6 .6 1]}, ...
    'title', 'Within-condition similarities');

%% 5. Formal comparison to model RDMs
result = R.compare({'condition','bodysite'}, 'correlation_type', 'kendall_taua');
fprintf('\nModel-RDM correlations (Kendall tau-a):\n');
for i = 1:numel(result.candidate_names)
    fprintf('  %-12s r=%.3f, p=%.4g %s\n', result.candidate_names{i}, ...
        result.r_mean(i), result.relatedness_p_corr(i), ...
        ternary(result.relatedness_sig(i), '*', ''));
end
fprintf('  noise ceiling: [%.3f %.3f]\n', result.noise_ceiling(1), result.noise_ceiling(2));

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
