%% AcceptMap RSA pipeline
% Full representational-similarity analysis of the WASABI AcceptMap dataset.
%
% DESIGN
%   8 subjects. The experimental factor is ACCEPTANCE (0 = experience,
%   1 = accept). Shared-anchor bodysite design: every subject has a
%   "leftface" site plus ONE idiosyncratic "other" site. We recode bodysite
%   to {leftface, other}.
%
% IMPORTANT -- why this differs from DistractMap:
%   acceptance is CONFOUNDED with session (each session contains only one
%   acceptance level). So there is NO session-to-session replication axis and
%   classic RSM reliability across sessions is undefined. The natural analysis
%   is the CROSSNOBIS representational dissimilarity matrix over the 4
%   conditions {acceptance x bodysite_type}, with cross-validation folds built
%   from the within-condition image repeats (fold_var='occurrence').
%   This reproduces the generate_RSA_accept_crossnobis design.
%
% PIPELINE
%   1. Load + prepare metadata
%   2. Build the crossnobis RDM (per subject + group)
%   3. Visualize the representational geometry
%   4. Cell-level contrasts on the RDM (acceptance vs bodysite effects)
%   5. Compare the brain RDM to model RDMs (acceptance / bodysite models)
%   6. (Optional) Searchlight RSA of the acceptance model

%% 1. Load and prepare ----------------------------------------------------
datafile = "\\dartfs-hpc\rc\lab\C\CANlab\labdata\projects\WASABI\WASABI_N_of_Few\analysis\04302026 MultiModal Run Maps.mat";
load(datafile, 'acceptmap_run');

T = acceptmap_run.metadata_table;
T.bodysite_type = rsa_recode_reference(T.bodysite, 'leftface', 'other_label', 'other');
acceptmap_run.metadata_table = T;

fprintf('AcceptMap: %d images, %d subjects, acceptance levels: %s\n', ...
    height(T), numel(unique(T.subject_id)), num2str(unique(T.acceptance)'));

%% 2. Build the crossnobis RDM -------------------------------------------
% 4 conditions: {acc0,acc1} x {leftface,other}. Folds = within-condition
% occurrence rank (no metadata column lines the conditions into shared folds,
% so use the 'occurrence' auto-mode, ordered by run_number).
R = compute_rsm(acceptmap_run, 'group_by', {'acceptance','bodysite_type'}, ...
    'subject_var', 'subject_id', 'metric', 'crossnobis', ...
    'fold_var', 'occurrence', 'run_var', 'run_number', 'level', 'subject');

fprintf('Crossnobis RDM: [%s] (4 conditions x %d subjects)\n', num2str(size(R)), size(R,3));
disp('Labels:'); disp(R.labels');

%% 3. Visualize the representational geometry ----------------------------
figure('Name','AcceptMap group RDM');
plot(mean(R));
title('AcceptMap group-mean crossnobis RDM');

figure('Name','AcceptMap RDM MDS');
plot(mean(R), 'mode', 'mds');
title('AcceptMap MDS (4 conditions)');

% Group-mean RDM as a table
m = mean(R.dat, 3, 'omitnan');
disp('Group-mean crossnobis RDM:');
disp(array2table(m, 'RowNames', matlab.lang.makeValidName(R.labels), ...
    'VariableNames', matlab.lang.makeValidName(R.labels)));

%% 4. Cell-level contrasts on the RDM ------------------------------------
% Auto-attached groupings include the acceptance levels (0, 1) and the
% bodysite_types (leftface, other). For a dissimilarity matrix, larger cells
% = more distinct representations.
disp('Available groupings:'); disp(fieldnames(R.groupings)');

% Cross-acceptance distance (same bodysite, different acceptance) vs
% cross-bodysite distance (same acceptance, different bodysite).
% In RDM space these are between-grouping cells.
v_cross_acc = R.cells('x0', 'x1', 'transform', 'none');   % acc0 vs acc1 cells
fprintf('\nMean cross-acceptance distance per subject:\n');
disp(v_cross_acc');
[~, p_acc, ~, st_acc] = ttest(v_cross_acc, 0, 'tail', 'right');
fprintf('Cross-acceptance distance > 0: t(%d)=%.2f, p=%.4g\n', st_acc.df, st_acc.tstat, p_acc);

%% 5. Compare brain RDM to model RDMs ------------------------------------
% Does the brain's geometry match an "acceptance" model (conditions cluster
% by acceptance) or a "bodysite" model (cluster by bodysite_type)?
result = R.compare({'acceptance','bodysite_type'}, ...
    'correlation_type', 'kendall_taua');
fprintf('\nModel-RDM comparison (Kendall tau-a):\n');
for i = 1:numel(result.candidate_names)
    star = ''; if result.relatedness_sig(i), star = '*'; end
    fprintf('  %-14s r=%+.3f  p=%.4g %s\n', result.candidate_names{i}, ...
        result.r_mean(i), result.relatedness_p_corr(i), star);
end
fprintf('  noise ceiling: [%.3f %.3f]\n', result.noise_ceiling(1), result.noise_ceiling(2));
if ~isnan(result.differences_p(1,2))
    fprintf('  acceptance vs bodysite model difference: p=%.4g\n', result.differences_p(1,2));
end

%% 6. Parcelwise brain maps ----------------------------------------------
% Where in the brain is the acceptance distinction represented? Run the
% crossnobis RDM per parcel and map the cross-acceptance vs cross-bodysite
% contrast. (Groupings here are the acceptance levels x0/x1 and bodysite
% types leftface/other; reference those auto-built names.)
atlas = load_atlas('canlab2024');
results = acceptmap_run.rsa_parcelwise('atlas', atlas, ...
    'group_by', {'acceptance','bodysite_type'}, 'subject_var', 'subject_id', ...
    'metric', 'crossnobis', 'fold_var', 'occurrence', 'run_var', 'run_number', ...
    'contrasts', {'crossAcc_vs_crossBodysite', {'x0','x1'}, {'leftface','other'}}, ...
    'correction', 'fdr', 'tail', 'right');
% NOTE: statistic_image has BOTH a `threshold` property and a `threshold`
% method, so the dot form `map.threshold(0.05,'unc')` indexes the property
% (badsubscript error). Always call threshold() with function syntax.
map = threshold(results.maps.crossAcc_vs_crossBodysite, 0.05, 'unc');
montage(map);

%% 7. Searchlight RSA -----------------------------------------------------
% Where does the LOCAL geometry match the acceptance model? Restrict to a
% gray-matter mask for speed.
gray = fmri_data(which('gray_matter_mask.img'), 'noverbose');
sl = acceptmap_run.searchlight_rsa('acceptance', 'radius', 3, ...
    'group_by', {'acceptance','bodysite_type'}, 'subject_var', 'subject_id', ...
    'metric', 'correlation', 'compare', 'spearman', 'mask', gray);

% IMPORTANT: with permutations=0 (default) the map holds RAW correlations and
% its .p values are placeholders (all 1). View it DIRECTLY -- do NOT call
% threshold(map,0.05,'unc') here (you'd get an empty map / "no results").
montage(sl.maps.acceptance);                      % raw r map (correct view)

% For p-thresholded inference, pass permutations>0 (slower); then threshold works:
%   sl = acceptmap_run.searchlight_rsa('acceptance', 'radius', 3, ...
%       'group_by', {'acceptance','bodysite_type'}, 'subject_var', 'subject_id', ...
%       'metric', 'correlation', 'compare', 'spearman', 'mask', gray, ...
%       'permutations', 500);
%   montage(threshold(sl.maps.acceptance, 0.05, 'unc'));

fprintf('\nAcceptMap pipeline complete.\n');
