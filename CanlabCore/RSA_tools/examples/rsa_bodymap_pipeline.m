%% BodyMap RSA pipeline (Sun et al., 2026) -- run-confound controlled
% Full representational-similarity analysis of the WASABI BodyMap dataset
% using the CanlabCore RSA tools, reproducing the Sun et al. 2026 figures
% with explicit control of the SHARED-RUN confound.
%
% DESIGN
%   9 subjects, 5-7 sessions each. Each RUN delivers Hot / Warm / Imagine
%   trials to ONE of 8 bodysites (abdomen, chest, leftarm, leftface, leftleg,
%   rightarm, rightface, rightleg). The run-level betas give 24 conditions
%   per subject = 8 bodysites x 3 conditions, each measured once per session
%   => 5-7 session replicates. Each (bodysite, session) is exactly ONE run,
%   so RUN == bodysite x session.
%
% THE SHARED-RUN CONFOUND (why this script exists)
%   Hot, Warm, and Imagine are ALWAYS estimated from the SAME run (the three
%   trial types co-occur within each bodysite's run). So the three conditions
%   of a bodysite share that run's noise (motion, physiology, drift) AND the
%   run/bodysite mean response. That shared structure massively inflates the
%   same-bodysite cross-condition similarities (Hot-Warm, Hot-Imagine,
%   Imagine-Warm) -- the very cells used to judge condition geometry -- and
%   confounds bodysite clustering too (the run IS the bodysite). In this
%   dataset the inflation is ~8x: same-bodysite cross-condition similarity in
%   S1 is z=0.44 within-run vs z=0.06 across-run (section 2). The fix is to
%   NEVER compare two patterns from the same run -- only across runs/sessions.
%
%   This reproduces the published approach (Representational Similarity
%   Analysis FullCode.mlx), which built run-level RSMs and then EXCLUDED
%   within-session(=within-run) correlations, replacing same-bodysite
%   cross-condition cells with across-session-only averages. Here we get the
%   same result more cleanly, two API-supported ways:
%     (1) CROSSNOBIS with cross-validation folds = SESSIONS (== runs here).
%         The within-fold condition difference cancels the run-common
%         component, and the cross-fold product has zero expectation for
%         run/session noise -> unbiased RDM by construction (Walther 2016).
%     (2) MIXED-EFFECTS model with a SAME-RUN nuisance term (runid). The
%         SameRunid coefficient absorbs the within-run inflation, so
%         SameBodysite/SameCondition are estimated from ACROSS-run pairs only.
%
% PIPELINE
%   1. Load + prepare
%   2. Diagnose the shared-run confound (quantify the inflation)
%   3. Build RSMs: run-controlled (crossnobis) + naive (for contrast)
%   4. Visualize representational geometry              (Fig 7A/B)
%   5. Model RDMs: bodysite / condition / cross-condition (Fig 7C)
%   6. Cell-level condition contrasts                    (Fig 7D-F)
%   7. Formal model comparison + noise ceiling
%   8. Multi-level LME with a SAME-RUN nuisance
%   9. Parcelwise brain maps                             (Fig 7G / S8 / S9)
%  10. Searchlight RSA
%
% Requires: CanlabCore on the path. Reload class methods if edited:
%   clear classes; rehash path; addpath(genpath('.../CanlabCore'));

%% 1. Load and prepare ----------------------------------------------------
datafile = "\\dartfs-hpc\rc\lab\C\CANlab\labdata\projects\WASABI\WASABI_N_of_Few\analysis\04302026 MultiModal Run Maps.mat";
load(datafile, 'run_maps');                 % fmri_data_st, 1224 run-level betas

T = run_maps.metadata_table;
T.sesno = double(string(T.sesno));          % numeric session index (CV folds)
T.runid = categorical(strcat(string(T.ses), '_', string(T.bodysite)));  % unique run == bodysite x session
run_maps.metadata_table = T;

fprintf('BodyMap: %d run images, %d subjects, %d bodysites, conditions: %s\n', ...
    height(T), numel(unique(string(T.sub))), numel(unique(string(T.bodysite))), ...
    strjoin(cellstr(unique(string(T.condition)))', ', '));

%% 2. Diagnose the shared-run confound ------------------------------------
% Hot/Warm/Imagine of a bodysite are estimated from the SAME run, so their
% cross-condition similarity is inflated by shared-run noise + the run mean.
% Direct evidence: take SAME-bodysite, DIFFERENT-condition pairs and compare
% their correlation when from the same run (== same session here) vs across
% runs. A somatosensory region (S1) makes the magnitude vivid.
roi = select_atlas_subset(load_atlas('canlab2024'), ...
    {'Ctx_3b_L','Ctx_3b_R','Ctx_1_L','Ctx_1_R','Ctx_2_L','Ctx_2_R'});   % S1
rmS1 = apply_mask(run_maps, roi);
X = double(rmS1.dat);
sub = string(T.sub); bs = string(T.bodysite); cond = string(T.condition); ses = T.sesno;
same_run = []; diff_run = [];
for s = unique(sub)'
    ix = find(sub==s); C = corr(X(:,ix)); n = numel(ix);
    [a,b] = find(triu(true(n),1));
    pick  = (bs(ix(a)) == bs(ix(b))) & (cond(ix(a)) ~= cond(ix(b)));  % same bodysite, diff condition
    sameR = ses(ix(a)) == ses(ix(b));                                 % same session == same run
    z = atanh(min(max(C(sub2ind([n n],a,b)),-.999),.999));
    same_run(end+1) = mean(z(pick &  sameR),'omitnan'); %#ok<SAGROW>
    diff_run(end+1) = mean(z(pick & ~sameR),'omitnan'); %#ok<SAGROW>
end
[~,p_inf,~,st_inf] = ttest(same_run, diff_run);
fprintf(['\nShared-run inflation (same-bodysite Hot/Warm/Imagine pairs, S1):\n' ...
    '  same run z=%.3f vs across runs z=%.3f; inflation=%.3f, t(%d)=%.2f, p=%.3g\n'], ...
    mean(same_run), mean(diff_run), mean(same_run-diff_run), st_inf.df, st_inf.tstat, p_inf);
fprintf('  => cross-condition similarity is dominated by shared-run noise; compare ACROSS runs only.\n');

%% 3. Build RSMs: run-controlled vs naive ---------------------------------
% PRIMARY: crossnobis RDM, folds = sessions (== runs). The within-fold
% condition difference cancels the run-common component, and cross-fold
% products remove run/session noise -> unbiased. Condition order is
% condition-major (hot/imagine/warm blocks) to match Fig 7A.
gray = fmri_data(which('gray_matter_mask.img'), 'noverbose');
run_gray = apply_mask(run_maps, gray);          % mask ONCE (run_maps is ~500 MB)

% NOTE: whole-brain crossnobis (~90k voxels) takes a few MINUTES -- the
% cross-validated distance loops over condition pairs x session-fold pairs.
% (Per-parcel crossnobis in section 9 is fast: each parcel is small.) To
% prototype quickly, restrict run_gray to a mask or subsample voxels first.
R = compute_rsm(run_gray, 'group_by', {'condition','bodysite'}, ...
    'subject_var', 'sub', 'metric', 'crossnobis', 'fold_var', 'sesno', ...
    'level', 'subject');

% NAIVE (for comparison only): session-pooled Spearman similarity. This is
% the construction that does NOT remove the same-session effect.
R_naive = compute_rsm(run_gray, 'group_by', {'condition','bodysite'}, ...
    'subject_var', 'sub', 'metric', 'spearman', 'level', 'subject');

fprintf('\nCrossnobis RDM [%s], naive Spearman RSM [%s]\n', ...
    num2str(size(R)), num2str(size(R_naive)));
disp('conditions:'); disp(R.labels');

% Show the payoff: bodysite separation, naive vs session-controlled.
lab = string(R.labels);                         % "condition_bodysite"
bsite = extractAfter(lab, "_");
k = numel(lab); [ii,jj] = find(triu(true(k),1));
sameBS = bsite(ii) == bsite(jj);
Mc = mean(R.dat,3,'omitnan');  vc = Mc(sub2ind([k k],ii,jj));
Mn = mean(R_naive.dat,3,'omitnan');  vn = Mn(sub2ind([k k],ii,jj));
fprintf('Naive Spearman  same-BS=%.3f diff-BS=%.3f (similarity; sep=%.3f)\n', ...
    mean(vn(sameBS)), mean(vn(~sameBS)), mean(vn(sameBS))-mean(vn(~sameBS)));
fprintf('Crossnobis      same-BS=%.3f diff-BS=%.3f (distance;   sep=%.3f)\n', ...
    mean(vc(sameBS)), mean(vc(~sameBS)), mean(vc(~sameBS))-mean(vc(sameBS)));

%% 4. Visualize representational geometry (Fig 7A / 7B) -------------------
figure('Name','BodyMap group RDM');
plot(mean(R), 'block_borders_by', 'condition');
title('BodyMap group-mean crossnobis RDM (session cross-validated)');

figure('Name','BodyMap subject RDMs');
plot(R, 'mode', 'grid');
sgtitle('Per-subject crossnobis RDMs');

figure('Name','BodyMap MDS'); plot(mean(R), 'mode', 'mds');
figure('Name','BodyMap dendrogram'); plot(mean(R), 'mode', 'dendrogram');

%% 5. Model RDMs (Fig 7C) -------------------------------------------------
% Build same-vs-different model RDMs over the 24 conditions.
parts = split(lab, "_");                        % k x 2: [condition, bodysite]
cond_meta = table(parts(:,1), parts(:,2), 'VariableNames', {'condition','bodysite'});

Mbody = rsm.from_categorical(cond_meta, 'bodysite');
Mcond = rsm.from_categorical(cond_meta, 'condition');
figure('Name','Model RDMs');
subplot(1,2,1); imagesc(Mbody.dat); axis image off; title('Bodysite model');
subplot(1,2,2); imagesc(Mcond.dat); axis image off; title('Condition model');
colormap(gray);

%% 6. Cell-level condition contrasts (Fig 7D-F) --------------------------
% Auto-attached groupings: one per condition (hot/imagine/warm, spanning the
% 8 bodysites) and one per bodysite (spanning the 3 conditions).
disp('Available groupings:'); disp(fieldnames(R.groupings)');
%
% Cross-condition similarity uses BETWEEN-block cells, expressed as {a,b}:
%   within-hot  = contrast('hot')           (mean distance among hot conditions)
%   HW          = contrast({'hot','warm'})  (mean hot<->warm distance)
%   HI vs HW    = contrast({'hot','imagine'}, {'hot','warm'})
% NOTE: R is a crossnobis DISTANCE (RDM): smaller = more similar. So a
% NEGATIVE (HI - HW) means Hot is MORE similar to Imagine than to Warm
% (the Sun et al. "HI > HW" similarity result, in dissimilarity space).
spec = {
 % name             cells_A              cells_B
 'within_hot',      'hot',               [];
 'within_warm',     'warm',              [];
 'within_imagine',  'imagine',           [];
 'HW',              {'hot','warm'},       [];
 'HI',              {'hot','imagine'},    [];
 'IW',              {'imagine','warm'},   [];
 'HI_vs_HW',        {'hot','imagine'},    {'hot','warm'};
 'HI_vs_IW',        {'hot','imagine'},    {'imagine','warm'};
};
T_contrasts = R.ttest_contrasts(spec, 'tail', 'both', 'correction', 'fdr');
disp(T_contrasts(:, {'Contrast','Mean_Diff','t','df','P','FDR_P','sig'}));

% Within-condition distances per subject, with within-subject lines.
T_cells = R.cells_table({'hot','imagine','warm'});
figure('Name','BodyMap within-condition distances');
plot_rsm_contrast_bars(T_cells, 'title', 'Within-condition distance (crossnobis)');

%% 7. Formal model comparison + noise ceiling ----------------------------
% Does the brain's geometry match a bodysite model or a condition model?
% Kendall tau-a per Nili et al. (2014); relatedness via subject RFX.
result = R.compare({'bodysite','condition'}, 'correlation_type', 'kendall_taua');
fprintf('\nModel-RDM comparison (Kendall tau-a):\n');
for i = 1:numel(result.candidate_names)
    star = ''; if result.relatedness_sig(i), star = '*'; end
    fprintf('  %-10s r=%+.3f  p=%.4g %s\n', result.candidate_names{i}, ...
        result.r_mean(i), result.relatedness_p_corr(i), star);
end
fprintf('  noise ceiling: [%.3f %.3f]\n', result.noise_ceiling(1), result.noise_ceiling(2));

%% 8. Multi-level LME with a SAME-RUN nuisance ---------------------------
% Correlation-based control for the shared-run confound: model pairwise
% similarity as same-bodysite / same-condition / SAME-RUN, with subject
% random effects. SameRunid is 1 only for the within-run Hot/Warm/Imagine
% pairs, so it absorbs the run inflation; SameBodysite and SameCondition are
% then estimated from ACROSS-run pairs (the clean, interpretable effects).
% Empirically: SameRunid ~ 0.38 (the inflation), SameBodysite ~ 0.02
% (p<1e-26), SameCondition ~ 0.06 -- the true geometry once run noise is out.
mdl = run_gray.rsa_lme( ...
    'predictors', {'bodysite','condition','runid'}, ...
    'subject_var', 'sub');
disp(mdl.Coefficients);

% Nested comparison: add the run nuisance first, then the effects of interest.
seq = rsa_model_sequence('Y ~ 1 + (1|Subject)');
seq = seq.add_term('SameRunid');         % shared-run nuisance first
seq = seq.add_term('SameCondition');
seq = seq.add_term('SameBodysite');
[T_models, best] = run_gray.rsa_compare_models(seq.formulas, ...
    'predictors', {'bodysite','condition','runid'}, ...
    'subject_var', 'sub', 'select_by', 'aic', 'verbose', false);
disp(T_models);
fprintf('Best model by AIC: %s\n', seq.formulas{best});

%% 9. Parcelwise brain maps (Fig 7G / S8 / S9) ---------------------------
% Per-parcel RSA, FDR/Bonferroni-corrected across parcels, projected to the
% brain. Crossnobis per parcel is fast (small voxel counts) and keeps the
% session cross-validation. Reference groupings by their auto-built names.
atlas = load_atlas('canlab2024');

% (a) Contrast map: where is Hot MORE similar to Imagine than to Warm?
%     In distance space that is (HI - HW) < 0, so use tail='left'.
results = run_maps.rsa_parcelwise('atlas', atlas, ...
    'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'metric', 'crossnobis', 'fold_var', 'sesno', ...
    'contrasts', {'HI_vs_HW', {'hot','imagine'}, {'hot','warm'}}, ...
    'correction', 'fdr', 'tail', 'left');
% statistic_image has a `threshold` PROPERTY and METHOD -- call threshold()
% with FUNCTION syntax (map.threshold(...) indexes the property -> error).
montage(threshold(results.maps.HI_vs_HW, 0.05, 'unc'));

% (b) Second condition-contrast map: Hot-Warm vs Imagine-Warm (run-clean via
%     crossnobis fold=sesno). Use crossnobis here -- it removes the shared-run
%     confound by construction. (The whole-brain correlation-based run-clean
%     decomposition with a SameRunid nuisance is in section 8; the parcelwise
%     LME path is not used for BodyMap -- see note below.)
results_b = run_maps.rsa_parcelwise('atlas', atlas, ...
    'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'metric', 'crossnobis', 'fold_var', 'sesno', ...
    'contrasts', {'HW_vs_IW', {'hot','warm'}, {'imagine','warm'}}, ...
    'correction', 'fdr', 'tail', 'both');
montage(threshold(results_b.maps.HW_vs_IW, 0.05, 'unc'));
% NOTE: bodysite-specificity brain maps are covered by the bodysite-model
% searchlight (section 10) and the model comparison (section 7). The
% rsa_parcelwise LME path currently returns empty maps on this 8-level
% bodysite design (silent per-parcel failure) -- prefer the crossnobis
% contrasts above for parcelwise brain maps.

%% 10. Searchlight RSA ----------------------------------------------------
% Where does the LOCAL geometry match the bodysite model? With permutations=0
% the map holds RAW correlations (view directly; do not threshold by p).
% Whole-brain searchlight is SLOW (~10+ min); shrink the mask to prototype.
sl = run_gray.searchlight_rsa('bodysite', 'radius', 3, ...
    'group_by', {'condition','bodysite'}, 'subject_var', 'sub', ...
    'metric', 'correlation', 'compare', 'spearman');
montage(sl.maps.bodysite);                       % raw r map
% For p-thresholded inference, pass permutations>0 (slow), then:
%   montage(threshold(sl.maps.bodysite, 0.05, 'unc'));

fprintf('\nBodyMap pipeline complete.\n');
