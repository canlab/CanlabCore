%% DistractMap RSA pipeline
% Full representational-similarity analysis of the WASABI DistractMap dataset
% using the CanlabCore RSA tools.
%
% DESIGN
%   7 subjects, 2 sessions each, 4 conditions
%   (nback_heat, nback_nostimblock, rest_nostim, rest_stim).
%   Shared-anchor bodysite design: every subject has a "Left Face" bodysite
%   plus ONE idiosyncratic "other" site that differs per subject. We recode
%   bodysite to {Left Face, Other Body Site} so conditions are comparable
%   across subjects (see rsa_recode_reference).
%
% PIPELINE
%   1. Load + prepare metadata
%   2. Build subject- and session-level RSMs
%   3. Visualize representational geometry
%   4. Reliability (ICC across sessions)
%   5. Within / between condition contrasts
%   6. Multi-level LME modeling
%   7. (Optional) Parcelwise brain maps
%
% Requires: CanlabCore on the path. Reload class methods if edited:
%   clear classes; rehash path; addpath(genpath('.../CanlabCore'));

%% 1. Load and prepare ----------------------------------------------------
datafile = "\\dartfs-hpc\rc\lab\C\CANlab\labdata\projects\WASABI\WASABI_N_of_Few\analysis\04302026 MultiModal Run Maps.mat";
load(datafile, 'distractmap_run');

T = distractmap_run.metadata_table;

% Two-level session index per subject (0 = earlier, 1 = later)
G = findgroups(T.subject_id);
min_ses = splitapply(@min, T.session_number, G);
T.sesno = double(T.session_number ~= min_ses(G));

% Recode bodysite into the shared anchor vs the idiosyncratic "other"
T.bodysite_type = rsa_recode_reference(T.bodySite, 'Left Face', ...
                                       'other_label', 'Other Body Site');

distractmap_run.metadata_table = T;

fprintf('DistractMap: %d images, %d subjects, conditions: %s\n', ...
    height(T), numel(unique(T.subject_id)), ...
    strjoin(cellstr(string(unique(T.condition)))', ', '));

%% 2. Build RSMs ----------------------------------------------------------
% Subject-collapsed RSM (sessions averaged) -- for contrasts / geometry
R = compute_rsm(distractmap_run, 'group_by', {'condition','bodysite_type'}, ...
    'subject_var', 'subject_id', 'metric', 'spearman', ...
    'nan_policy', 'skip_replicate');

% Per-session RSM -- for reliability and drift
R_sess = compute_rsm(distractmap_run, 'group_by', {'condition','bodysite_type'}, ...
    'subject_var', 'subject_id', 'session_var', 'session_number', ...
    'level', 'session', 'metric', 'spearman', 'nan_policy', 'skip_replicate');

fprintf('Subject RSM: [%s];  Session RSM: [%s]\n', ...
    num2str(size(R)), num2str(size(R_sess)));

%% 3. Visualize representational geometry --------------------------------
figure('Name','DistractMap group RSM');
plot(mean(R), 'block_borders_by', 'condition');
title('DistractMap group-mean RSM');

figure('Name','DistractMap subject RSMs');
plot(R, 'mode', 'grid');
sgtitle('Per-subject RSMs');

% MDS of the 8 conditions -- shows how conditions cluster in representational
% space (nback vs rest; Left Face vs Other within each condition).
figure('Name','DistractMap MDS');
plot(mean(R), 'mode', 'mds');

% Hierarchical clustering of the conditions
figure('Name','DistractMap dendrogram');
plot(mean(R), 'mode', 'dendrogram');

%% 4. Reliability (ICC across sessions) ----------------------------------
% Whole-RSM reliability per subject, aggregated
rel = R_sess.reliability('icc_type', '3-k');
fprintf('\nWhole-RSM ICC across %d subjects: mean=%.3f (median=%.3f)\n', ...
    rel.summary.n_subjects, rel.summary.mean, rel.summary.median);

% Per-condition reliability + bar plot
out = R_sess.reliability_per_condition('icc_type', '3-k');
figure('Name','DistractMap reliability');
barplot_columns(out.Mean_ICC', 'names', format_strings_for_legend(out.Condition'), ...
    'custom_se', out.Std_ICC', 'noviolin', 'nofig');
ylabel('Reliability (ICC 3-k)'); title('DistractMap per-condition reliability');
hline(0, 'k-');
disp(out(:, {'Condition','Mean_ICC','Std_ICC','n_subjects'}));
% NOTE: only 2 sessions/subject -> ICC noisy; rest_nostim has little task
% structure to replicate (expect low/negative ICC).

%% 5. Within / between condition contrasts -------------------------------
% compute_rsm auto-attaches groupings: one per condition (nback_heat, ...)
% spanning that condition's two bodysite_types, and one per bodysite_type
% (Left_Face, Other_Body_Site) spanning the four conditions.
disp('Available groupings:'); disp(fieldnames(R.groupings)');

% Contrast battery over the superordinate conditions. cells(A,A) is the
% within-condition LeftFace<->Other similarity; cells(A,B) is between-condition.
spec = {
    % name                  cells_A         cells_B
    'within_nback_heat',    'nback_heat',    [];
    'within_nback_nostim',  'nback_nostimblock', [];
    'within_rest_stim',     'rest_stim',     [];
    'within_rest_nostim',   'rest_nostim',   [];
    'heat_vs_reststim',     'nback_heat',    'rest_stim';
    'nbackHeat_vs_nostim',  'nback_heat',    'nback_nostimblock';
};
T_contrasts = R.ttest_contrasts(spec, 'tail', 'both', 'correction', 'fdr');
disp(T_contrasts);

% Per-subject cell table for plotting
T_cells = R.cells_table({'nback_heat','nback_nostimblock','rest_nostim','rest_stim'});
figure('Name','DistractMap within-condition similarity');
plot_rsm_contrast_bars(T_cells, 'title', 'Within-condition (LeftFace <-> Other) similarity');

%% 6. Multi-level LME modeling -------------------------------------------
% Model pairwise similarity as a function of same-condition / same-bodysite-
% type / same-session, with subject as a random effect.
mdl = distractmap_run.rsa_lme( ...
    'predictors', {'condition','bodysite_type','sesno'}, ...
    'subject_var', 'subject_id');
disp(mdl.Coefficients);

% Variance components
icc = rsa_lme_icc(mdl);
disp(icc.summary);

% Nested model comparison (subject_var 'subject_id' -> short name 'Subject')
seq = rsa_model_sequence('Y ~ 1 + (1|Subject)');
seq = seq.add_term('SameCondition');
seq = seq.add_term('SameBodysitetype');
seq = seq.add_term('SameSesno');
[T_models, best] = distractmap_run.rsa_compare_models(seq.formulas, ...
    'predictors', {'condition','bodysite_type','sesno'}, ...
    'subject_var', 'subject_id', 'select_by', 'aic', 'verbose', false);
disp(T_models);
fprintf('Best model by AIC: %s\n', seq.formulas{best});

%% 7. Parcelwise brain maps ----------------------------------------------
% Per-parcel RSA inference, FDR-corrected across parcels, projected onto the
% brain as statistic_image maps. Reference groupings by the auto-built names
% (nback_heat, rest_stim, ...) -- NOT the composite labels.
atlas = load_atlas('canlab2024');

% (a) Contrast map: within-condition similarity for nback_heat vs rest_stim
results = distractmap_run.rsa_parcelwise('atlas', atlas, ...
    'group_by', {'condition','bodysite_type'}, 'subject_var', 'subject_id', ...
    'metric', 'spearman', ...
    'contrasts', {'nbackheat_vs_reststim', 'nback_heat', 'rest_stim'}, ...
    'correction', 'fdr', 'tail', 'right');
% NOTE: statistic_image has BOTH a `threshold` property and a `threshold`
% method, so the dot form `map.threshold(0.05,'unc')` indexes the property
% (badsubscript error). Always call threshold() with function syntax.
map = threshold(results.maps.nbackheat_vs_reststim, 0.05, 'unc');
montage(map);

% (b) LME map: where is the SameCondition effect strongest?
results_lme = distractmap_run.rsa_parcelwise('atlas', atlas, ...
    'group_by', {'condition','bodysite_type'}, 'subject_var', 'subject_id', ...
    'lme', 'Y ~ SameCondition + SameBodysitetype + (1|Subject)', ...
    'predictors', {'condition','bodysite_type'});
montage(results_lme.maps.SameCondition);

%% 8. Searchlight RSA -----------------------------------------------------
% Where does the LOCAL representational geometry match the condition model?
% A spherical searchlight builds an RSM in each sphere and correlates it with
% the model RDM (positive = local geometry matches the model). Restrict to a
% gray-matter mask for speed.
gray = fmri_data(which('gray_matter_mask.img'), 'noverbose');
sl = distractmap_run.searchlight_rsa('condition', 'radius', 3, ...
    'group_by', {'condition','bodysite_type'}, 'subject_var', 'subject_id', ...
    'metric', 'correlation', 'compare', 'spearman', 'mask', gray);

% IMPORTANT: with permutations=0 (default) the map holds RAW correlations and
% its .p values are placeholders (all 1). View it DIRECTLY -- do NOT call
% threshold(map,0.05,'unc') here (p<0.05 is never true, so you'd get an empty
% map / "no results"). To emphasize strong matches, threshold by magnitude.
montage(sl.maps.condition);                       % raw r map (correct view)
% slmap = sl.maps.condition; slmap.dat(slmap.dat < 0.2) = 0; montage(slmap);

% For p-thresholded inference, pass permutations>0 (slower). Then .p is a real
% permutation p-value and the usual threshold() call works:
%   sl = distractmap_run.searchlight_rsa('condition', 'radius', 3, ...
%       'group_by', {'condition','bodysite_type'}, 'subject_var', 'subject_id', ...
%       'metric', 'correlation', 'compare', 'spearman', 'mask', gray, ...
%       'permutations', 500);
%   montage(threshold(sl.maps.condition, 0.05, 'unc'));

fprintf('\nDistractMap pipeline complete.\n');
