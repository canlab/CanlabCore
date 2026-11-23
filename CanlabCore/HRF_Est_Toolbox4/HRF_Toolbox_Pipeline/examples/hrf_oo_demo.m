%% HRF Toolbox -- an object-oriented tour
% How the HRF estimation pipeline plugs into CanlabCore's object model, using
% the same operations as a real study (the WASABI DistractMap "HRF Estimation
% Animations and Atlas Regions" analysis): fit HRF models, carry the whole-brain
% condition x lag estimates as first-class objects, animate them over the
% peristimulus window, summarize them by atlas region, and take them to the
% group -- all off the standard fmri_data / statistic_image surface.
%
%   image_vector (abstract superclass)
%     |
%     +-- fmri_data        --[ subclass ]--> fmri_hrf       (HRF beta / amplitude maps)
%     +-- statistic_image  --[ subclass ]--> statistic_hrf  (HRF t / p / SE maps)
%
%   timeseries --[ hrf_fit_all_models ]-------> fitted HRF curves (+ model SE)
%   wholebrain --[ make_fmri_stat_hrf ]-------> paired (fmri_hrf, statistic_hrf)
%   fmri_hrf   --[ hrf_make_montage_animation ]--> condition x lag brain movie
%   fmri_hrf   --[ apply_parcellation ]--------> per-region HRF curves
%   subject stack --[ ttest ]-----------------> group statistic_hrf (t maps)
%   study      --[ hrf_time_unfolding_stats ]--> corrected per-lag group inference
%
% Both classes ARE fmri_data / statistic_image subclasses, so every inherited
% method (montage, threshold, region, apply_parcellation, extract_roi_averages,
% ...) works on HRF estimates unchanged. Sections 1-7 run in-memory on a
% real brain geometry; Section 8 is the on-disk study API with the actual
% DistractMap calls.
%
% Requires: CanlabCore + HRF_Est_Toolbox4 + SPM12 (+ Neuroimaging_Pattern_Masks).
% Run top-to-bottom; brain movies are gated behind SHOW (default off).

SHOW = false;      % set true to render the brain montages / animations (slow)
TR   = 1;


%% 1. Fit HRF models to a timeseries: hrf_fit_all_models
% The functional entry point: given a 1D timeseries and per-condition stick
% functions, fit a family of HRF models (fir, sfir, canonical, spline, ...) with
% one consistent output API -- each carrying the fitted curve and, for linear
% models, model-based standard errors.
win = 30; T = 300; onsets = 10:30:T-40;
stick = zeros(T, 1); stick(onsets) = 1;
hrf_true = spm_hrf(TR);
bold = conv(stick, hrf_true); bold = bold(1:T);
rng(0); tc = bold + 0.4 * std(bold) * randn(T, 1);

fits = hrf_fit_all_models(tc, TR, {stick}, win, {'fir', 'sfir', 'canonical'});
fprintf('fitted models: %s ; sfir peak at lag %d, model SE available: %d\n', ...
    strjoin(fieldnames(fits), ', '), ...
    find(fits.sfir.hrf == max(fits.sfir.hrf), 1), isfield(fits.sfir, 'se'));

figure('Color', 'w'); hold on;
plot((0:numel(hrf_true)-1) * TR, hrf_true / max(hrf_true), 'k--', 'LineWidth', 1.5);
plot((0:numel(fits.sfir.hrf)-1) * TR, fits.sfir.hrf / max(fits.sfir.hrf), 'LineWidth', 1.5);
plot((0:numel(fits.fir.hrf)-1) * TR, fits.fir.hrf / max(fits.fir.hrf), 'LineWidth', 1);
legend({'true (canonical)', 'sfir fit', 'fir fit'}); xlabel('Seconds'); ylabel('Normalized HRF');
title('hrf\_fit\_all\_models -- recovered HRF shapes'); box off;


%% 2. fmri_hrf: whole-brain condition x lag betas as an fmri_data subclass
% A whole-brain fit produces one beta per (condition x lag) volume. fmri_hrf
% wraps those maps with the metadata to interpret them (condition, lag_index,
% lag_seconds) while remaining a full fmri_data underneath. We build one here in
% a real brain space (borrowing a sample dataset's geometry), planting a
% canonical HRF time-course into the 'pain' condition so the maps are realistic.
base = load_image_set('emotionreg'); nvox = size(base.dat, 1);
conds = {'pain', 'neutral'}; nLag = 20; nVol = numel(conds) * nLag;
hshape = spm_hrf(TR); hshape = hshape(1:nLag); hshape = hshape / max(hshape);
rng(1); B = 0.2 * randn(nvox, nVol);
cond_col = {}; lag_idx = []; lag_sec = []; img_lab = {}; col = 0;
for c = 1:numel(conds)
    for L = 1:nLag
        col = col + 1;
        if strcmp(conds{c}, 'pain'), B(:, col) = B(:, col) + hshape(L); end
        cond_col{end+1} = conds{c};  lag_idx(end+1) = L;                 %#ok<SAGROW>
        lag_sec(end+1) = (L-1) * TR; img_lab{end+1} = sprintf('%s_lag%02d', conds{c}, L); %#ok<SAGROW>
    end
end
mt = table(cond_col(:), lag_idx(:), lag_sec(:), img_lab(:), ...
    'VariableNames', {'condition', 'lag_index', 'lag_seconds', 'image_label'});
fd = base; fd.dat = B; fd.removed_images = 0; fd.image_names = ''; fd.fullpath = '';

Hb = fmri_hrf(fd, 'MetadataTable', mt, 'Subject', 'sub-01', 'RunLabel', 'run-01', ...
    'ModelName', 'sfir', 'TR', TR, 'Conditions', conds);
disp(Hb)
fprintf('isa(Hb, ''fmri_data'') = %d\n', isa(Hb, 'fmri_data'));


%% 3. statistic_hrf: the paired inferential object (statistic_image subclass)
% make_fmri_stat_hrf builds the beta side (fmri_hrf) and the t side
% (statistic_hrf) together, so their HRF metadata stays aligned. The t side is a
% statistic_image subclass -> threshold / region / montage apply.
bimg = statistic_image; bimg.dat = B; bimg.volInfo = base.volInfo;
bimg.removed_voxels = base.removed_voxels; bimg.removed_images = 0;
se = 0.3 + abs(randn(nvox, nVol));
timg = statistic_image; timg.dat = B ./ se; timg.volInfo = base.volInfo;
timg.removed_voxels = base.removed_voxels; timg.removed_images = 0; timg.type = 'T';

[Hb2, Ht] = make_fmri_stat_hrf(struct('b', bimg, 't', timg), ...
    'Subject', 'sub-01', 'ModelName', 'sfir', 'TR', TR, 'MetadataTable', mt);  %#ok<ASGLU>
disp(Ht)

% Slice the pain peak-lag t-map and push it through the standard stat chain.
% (The pipeline stores .p/.dfe; we set them here for this synthetic object.)
peakLag = find(hshape == max(hshape), 1);
idx = find(strcmp(mt.condition, 'pain') & mt.lag_index == peakLag);
painT = get_wh_image(Ht, idx);
painT.dfe = 100; painT.p = 2 * (1 - tcdf(abs(painT.dat), painT.dfe));
painT = threshold(painT, 0.01, 'unc');
fprintf('pain peak-lag t-map: threshold -> %d regions\n', numel(region(painT)));
if SHOW, montage(painT); end   %#ok<UNRCH>  (SHOW is a demo toggle)


%% 4. HRF Estimation ANIMATION -- straight off the object
% hrf_make_montage_animation accepts any image_vector, so a fmri_hrf slice
% animates directly. Slice the 'pain' condition (its 20 lag volumes) and sweep
% the montage across the peristimulus window -- the DistractMap "HRF animation",
% now driven by the object instead of a NIfTI path.
painCond = get_wh_image(Hb, find(strcmp(mt.condition, 'pain')));   % 20-volume fmri_hrf
fprintf('pain condition slice: %d lag volumes to animate\n', size(painCond.dat, 2));
if SHOW    %#ok<UNRCH>
    hrf_make_montage_animation(painCond, 'pain_hrf.mp4', ...
        'FrameStep', 1, 'FPS', 8, 'TitlePrefix', 'pain lag', 'MontageType', 'compact2');
    % thresholded t-map animation is the same call on the statistic_hrf:
    % hrf_make_montage_animation(get_wh_image(Ht, find(strcmp(mt.condition,'pain'))), ...
    %     'pain_t.mp4', 'Threshold', 2.0, 'FPS', 12);
end


%% 5. Atlas Regions -- per-region HRF curves via inherited apply_parcellation
% Because fmri_hrf IS an fmri_data, apply_parcellation reduces the whole-brain
% condition x lag maps to region means: [volumes x regions]. Slice the 'pain'
% rows and you have one HRF curve per atlas region -- the object-level version
% of the score-CSV atlas curves that plot_hrf_atlas_curves renders, and a basis
% for the same RankBy summaries (peak, AUC, ...).
at = select_atlas_subset(load_atlas('canlab2024'), 1:40);
pm = apply_parcellation(Hb, at);                       % [nVol x nRegion]
painRows = pm(strcmp(mt.condition, 'pain'), :);        % [nLag x nRegion], one curve/region
region_peak = max(abs(painRows), [], 1);               % rank regions by |peak| (as RankBy='peak_abs')
[~, ord] = sort(region_peak, 'descend');
fprintf('mean region pain curve corr with true HRF = %.2f ; top region: %s\n', ...
    corr(mean(painRows, 2), hshape(:)), at.labels{ord(1)});

figure('Color', 'w'); tt = (0:nLag-1) * TR;
plot(tt, painRows(:, ord(1:min(6, end))), 'LineWidth', 1.3);
xlabel('Peristimulus time (s)'); ylabel('HRF beta'); box off;
legend(at.labels(ord(1:min(6, end))), 'Interpreter', 'none', 'Location', 'best');
title('Per-region pain HRF curves (top 6 by |peak|)');


%% 6. To the GROUP: a group statistic_hrf of one-sample t maps
% Stack subjects' condition x lag betas and one-sample-t each volume across
% subjects -> a group statistic_hrf (t/p per condition x lag). This is the
% object core of hrf_make_average_montage_animations, whose .group_t is exactly
% this kind of statistic_image. Vectorized here across all 40 volumes at once.
Nsub = 8; rng(5);
sigvox = false(nvox, 1); sigvox(1:6000) = true;        % a focal set of "responding" voxels
stack = zeros(nvox, nVol, Nsub);
for s = 1:Nsub
    signal = zeros(nvox, nVol); signal(sigvox, :) = B(sigvox, :);
    stack(:, :, s) = signal + 0.3 * randn(nvox, nVol);
end
mu = mean(stack, 3); sd = std(stack, 0, 3);
tg = mu ./ (sd ./ sqrt(Nsub));
gsi = statistic_image; gsi.dat = tg; gsi.volInfo = base.volInfo;
gsi.removed_voxels = base.removed_voxels; gsi.removed_images = 0; gsi.type = 'T';
gsi.dfe = Nsub - 1; gsi.p = 2 * (1 - tcdf(abs(tg), Nsub - 1));
Hg = statistic_hrf(gsi, 'MetadataTable', mt, 'Subject', 'group', ...
    'ModelName', 'sfir', 'TR', TR, 'Conditions', conds);
disp(Hg)

painGT = get_wh_image(Hg, idx);                         % group t at pain peak lag
painGT = threshold(painGT, 0.001, 'unc');
fprintf('group pain peak-lag t-map: %d suprathreshold voxels; the whole Hg animates like section 4\n', ...
    sum(painGT.sig(:, 1)));
if SHOW, montage(painGT); end   %#ok<UNRCH>


%% 7. Corrected per-lag group inference over the HRF timecourse
% hrf_time_unfolding_stats builds a per-subject [subject x lag] contrast and
% tests each lag. Small-n per-lag inference needs multiple-comparison control:
% 'Correction' routes the lags through the shared hrf_group_stats engine
% (sign-flip / label max-|t| FWER, temporal cluster-mass, or FDR).
nsubj = 9; nt = win; rng(7);
study = struct(); study.subject_ids = arrayfun(@(i) sprintf('sub-%02d', i), 1:nsubj, 'uni', 0);
study.results = cell(1, nsubj);
for s = 1:nsubj
    h = 0.02 * randn(nt, numel(conds)); h(8:14, 1) = h(8:14, 1) + 0.7;
    r0 = struct('conditions', {conds}); r0.fits = struct('sfir', struct('hrf', h));
    study.results{s} = r0;
end
stats = hrf_time_unfolding_stats(study, 'Model', 'sfir', ...
    'ConditionA', 'pain', 'ConditionB', 'neutral', ...
    'Correction', 'cluster', 'MissingPolicy', 'silent');
fprintf('per-lag significance: %d uncorrected vs %d cluster-corrected (lags %s)\n', ...
    sum(stats.significant), sum(stats.significant_corrected), ...
    num2str(find(stats.significant_corrected')));


%% 8. The full study API on disk (the actual DistractMap workflow)
% Sections 1-7 run in-memory; a real study fits every subject/run to whole-brain
% NIfTIs, scores them by atlas/signature/imageset, and indexes them in an
% input_table. Every feature below consumes that table or its output dirs. These
% are the exact calls from the DistractMap "Animations and Atlas Regions"
% analysis; run them against your own hrf_outputs directories.
%
%   at  = load_atlas('canlab2024');
%   src = struct('label', {'distractmap','distractmap'}, ...
%                'table', {input_table_lf, input_table_obs});   % pool 2 body sites
%
%   % --- Atlas-region HRF curves, ranked (RankBy: peak_abs|auc_abs|hrf_match|
%   %     peak_t|n_sig|snr|shape_r2|sd), with condition contrasts and non-atlas
%   %     sources (signatures, bucknerlab/marg/hansen22 imagesets):
%   plot_hrf_atlas_curves(src, 'AtlasObj', at, 'Model', 'sfir', 'Object', 'beta', ...
%       'Conditions', {'nback-stimblock_ttl_1','rest_stim_ttl_1'}, 'TopN', 20, 'RankBy', 'hrf_match');
%   plot_hrf_curves(src, 'Source', 'signature', 'Set', 'all', 'RankBy', 'hrf_match', 'TopN', 44);
%
%   % --- Per-curve shape metrics + group one-sample t (FDR across rows):
%   T = hrf_curve_summaries(input_table_lf, 'Sources', {'atlas_CANLab2024'}, ...
%           'Objects', {'beta'}, 'Conditions', {'nback-stimblock_ttl_1','rest_stim_ttl_1'});
%   G = hrf_curve_summary_groupstats(T, 'GroupBy', {'condition','model','source_name'}, ...
%           'Across', 'subject', 'Metrics', {'peak_value','peak_lag_seconds','auc'}, 'Correction', 'fdr');
%   G(G.sig, :)
%
%   % --- Group whole-brain t-map animations (returns statistic_image .group_t):
%   avg = hrf_make_average_montage_animations(input_table_lf, 'Model', 'sfir', 'Object', 'beta', ...
%       'Conditions', {'nback-stimblock_ttl_1','rest_stim_ttl_1'}, ...
%       'GroupStatistic', 't', 'GroupCorrection', 'unc', 'GroupAlpha', 0.001);
%   avg.group_t          % statistic_image objects -> threshold/montage/region
%
%   % --- Term / network / region wordcloud animations over the HRF, FWE-corrected:
%   hrf_animate_wordcloud({lf,obs}, 'Set','neurosynth_topics_ri', 'Model','sfir', ...
%       'Condition','nback-stimblock_ttl_1', 'Correction','permutation', 'Threshold',0.05);
%
%   % --- Directed connectivity (HRF-deconvolved Granger) + trial-level mediation:
%   R  = hrf_causality({lf,obs}, 'Unit','signature');   R.remove.net_group
%   Rr = hrf_causality({lf,obs}, 'Condition','rest_stim');
%   Rn = hrf_causality({lf,obs}, 'Condition','nback-stimblock');
%   hrf_plot_causality(hrf_causality_contrast(Rr, Rn, 'Label1','rest', 'Label2','nback'));
%   Rmed = hrf_causality_mediation({lf,obs}, 'X','temp', 'M','NPS', 'Y','rating', 'TrialType','*stimblock*');


%% 9. The takeaway
% One object model, one method vocabulary, extended to HRFs:
%
%   timeseries --hrf_fit_all_models-->  fitted curves (+ model SE)
%   wholebrain --make_fmri_stat_hrf-->  fmri_hrf (IS-A fmri_data)
%                                       statistic_hrf (IS-A statistic_image)
%   fmri_hrf --hrf_make_montage_animation--> condition x lag brain movie
%   fmri_hrf --apply_parcellation--> per-region HRF curves (RankBy summaries)
%   subject stack --ttest--> group statistic_hrf --threshold/montage/region
%   study --hrf_time_unfolding_stats(Correction)--> corrected per-lag inference
%
% fmri_hrf and statistic_hrf inherit the entire fmri_data / statistic_image
% surface, add condition x lag HRF metadata and an HRF-aware display, and every
% HRF analysis -- animation, atlas-region curves, group maps, wordclouds,
% causality -- lands back in a standard CANlab container.
disp('HRF object-oriented tour complete.');
