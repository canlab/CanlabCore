%% HRF Toolbox -- an object-oriented tour
% How the HRF estimation pipeline plugs into CanlabCore's object model. Two new
% classes carry HRF estimates as first-class citizens of the CANlab hierarchy,
% so everything you already know about fmri_data / statistic_image -- masking,
% plotting, thresholding, region extraction, montages, atlases -- works on HRF
% results unchanged, and the HRF-specific layer (condition x lag metadata, model
% fitting, deconvolution, causality, group inference) is added on top.
%
%   image_vector (abstract superclass)
%     |
%     +-- fmri_data        --[ subclass ]--> fmri_hrf       (HRF beta / amplitude maps)
%     +-- statistic_image  --[ subclass ]--> statistic_hrf  (HRF t / p / SE maps)
%
%   timeseries --[ hrf_fit_all_models ]--> fitted HRF curves (+ model SE)
%   wholebrain --[ make_fmri_stat_hrf ]--> paired (fmri_hrf, statistic_hrf)
%   fmri_hrf   --[ to_statistic_hrf ]-----> statistic_hrf
%   study      --[ hrf_time_unfolding_stats ]--> per-lag group inference (corrected)
%   timeseries --[ hrf_deconv_timeseries -> hrf_granger_causality ]--> directed flow
%
% Both classes are fmri_data / statistic_image SUBCLASSES: a fmri_hrf IS an
% fmri_data with condition x lag metadata, so any inherited method just works.
%
% Requires: CanlabCore + HRF_Est_Toolbox4 (legacy fitters) + SPM12 on the path
% (+ Neuroimaging_Pattern_Masks for load_atlas). Run top-to-bottom; heavy
% montages are gated behind SHOW (default off).

SHOW = false;      % set true to pop the brain montages / stat plots
TR   = 1;          % seconds


%% 1. Fit HRF models to a timeseries: hrf_fit_all_models
% The functional entry point. Given a 1D timeseries and per-condition stick
% functions, it fits a family of HRF models with one consistent output API
% (fir, sfir, canonical, spline, ...), each carrying the fitted curve and,
% for linear models, model-based standard errors.
win = 30;                                        % HRF window (seconds)
T = 300; onsets = 10:30:T-40;
stick = zeros(T, 1); stick(onsets) = 1;          % one condition's event sticks
hrf_true = spm_hrf(TR);                          % ground-truth HRF (SPM canonical)
bold = conv(stick, hrf_true); bold = bold(1:T);
rng(0); tc = bold + 0.4 * std(bold) * randn(T, 1);   % noisy BOLD timeseries

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


%% 2. fmri_hrf: HRF betas as a first-class fmri_data subclass
% A whole-brain fit produces one beta per (condition x lag) volume. fmri_hrf
% wraps those maps with the metadata needed to interpret them (condition,
% lag_index, lag_seconds) while remaining a full fmri_data underneath. We build
% one here in a real brain space (borrowing a sample dataset's geometry).
base = load_image_set('emotionreg');             % real volInfo / MNI space
nvox = size(base.dat, 1);
conds = {'pain', 'neutral'}; nLag = 20; nVol = numel(conds) * nLag;
rng(1); B = randn(nvox, nVol);                   % stand-in condition x lag betas
cond_col = {}; lag_idx = []; lag_sec = []; img_lab = {};
for c = 1:numel(conds)
    for L = 1:nLag
        cond_col{end+1} = conds{c};  lag_idx(end+1) = L;                 %#ok<SAGROW>
        lag_sec(end+1) = (L-1) * TR; img_lab{end+1} = sprintf('%s_lag%02d', conds{c}, L); %#ok<SAGROW>
    end
end
mt = table(cond_col(:), lag_idx(:), lag_sec(:), img_lab(:), ...
    'VariableNames', {'condition', 'lag_index', 'lag_seconds', 'image_label'});
fd = base; fd.dat = B; fd.removed_images = 0; fd.image_names = ''; fd.fullpath = '';

Hb = fmri_hrf(fd, 'MetadataTable', mt, 'Subject', 'sub-01', 'RunLabel', 'run-01', ...
    'ModelName', 'sfir', 'TR', TR, 'Conditions', conds);
disp(Hb)                                          % HRF-aware summary
fprintf('isa(Hb, ''fmri_data'') = %d\n', isa(Hb, 'fmri_data'));

% Because it IS an fmri_data, inherited methods just work. Slice one condition x
% lag volume by its metadata, and run an atlas parcellation on it:
idx = find(strcmp(mt.condition, 'pain') & mt.lag_index == 6);   % pain, peak lag
painBeta = get_wh_image(Hb, idx);                 % -> a 1-volume fmri_hrf
atl = select_atlas_subset(load_atlas('canlab2024'), 1:20);
parcel_means = apply_parcellation(painBeta, atl); % inherited fmri_data method
fprintf('apply_parcellation(fmri_hrf, atlas) -> %d parcel means\n', numel(parcel_means));


%% 3. statistic_hrf: the paired inferential object (subclass of statistic_image)
% make_fmri_stat_hrf builds the beta side (fmri_hrf) and the t side
% (statistic_hrf) together, so their HRF metadata is guaranteed aligned. The t
% side is a statistic_image subclass, so threshold / region / montage apply.
bimg = statistic_image; bimg.dat = B; bimg.volInfo = base.volInfo;
bimg.removed_voxels = base.removed_voxels; bimg.removed_images = 0;
se = 0.3 + abs(randn(nvox, nVol));
timg = statistic_image; timg.dat = B ./ se; timg.volInfo = base.volInfo;
timg.removed_voxels = base.removed_voxels; timg.removed_images = 0; timg.type = 'T';

[Hb2, Ht] = make_fmri_stat_hrf(struct('b', bimg, 't', timg), ...
    'Subject', 'sub-01', 'RunLabel', 'run-01', 'ModelName', 'sfir', ...
    'TR', TR, 'MetadataTable', mt);               %#ok<ASGLU>
disp(Ht)

% Slice the pain peak-lag t-map and push it through the standard stat chain.
% (The pipeline stores .p/.dfe; we set them here for this synthetic object.)
painT = get_wh_image(Ht, idx);
painT.dfe = 100; painT.p = 2 * (1 - tcdf(abs(painT.dat), painT.dfe));
painT = threshold(painT, 0.01, 'unc');            % same threshold() as any stat map
r = region(painT);                                % -> region objects
fprintf('pain peak-lag t-map: threshold -> %d regions (region/table/montage all apply)\n', numel(r));
if SHOW, montage(painT); end   %#ok<UNRCH>  (SHOW is a demo toggle)


%% 4. Group inference over the HRF timecourse, with correction
% hrf_time_unfolding_stats builds a per-subject [subject x lag] contrast and
% tests each lag. Small-n per-lag inference needs multiple-comparison control:
% 'Correction' routes the tested lags through the shared hrf_group_stats engine
% (sign-flip / label max-|t| FWER, temporal cluster-mass, or FDR).
nsub = 9; nt = win;
study = struct(); study.subject_ids = arrayfun(@(i) sprintf('sub-%02d', i), 1:nsub, 'uni', 0);
study.results = cell(1, nsub);
for s = 1:nsub
    h = 0.02 * randn(nt, numel(conds));           % [lag x condition] HRF betas
    h(8:14, 1) = h(8:14, 1) + 0.7;                % plant a pain response 8-14 s
    r0 = struct('conditions', {conds}); r0.fits = struct('sfir', struct('hrf', h));
    study.results{s} = r0;
end
stats = hrf_time_unfolding_stats(study, 'Model', 'sfir', ...
    'ConditionA', 'pain', 'ConditionB', 'neutral', ...
    'Correction', 'cluster', 'MissingPolicy', 'silent');
fprintf('per-lag significance: %d uncorrected vs %d cluster-corrected (lags %s)\n', ...
    sum(stats.significant), sum(stats.significant_corrected), ...
    num2str(find(stats.significant_corrected')));
if SHOW, plot_hrf_time_unfolding_stats(stats); end   %#ok<UNRCH>  (SHOW is a demo toggle)


%% 5. Beyond: deconvolution -> causality, misspecification, animation
% The same objects/timeseries feed the rest of the toolbox. These are pointers
% (they need real timeseries or the whole-brain score CSVs); see the runnable
% causality walkthrough in examples/hrf_causality_synth_demo.m.
%
%   % HRF-informed deconvolution then directed connectivity (removes the
%   % regional-HRF-latency confound of BOLD Granger causality; David et al. 2008):
%   xhat = hrf_deconv_timeseries(y, fits.sfir.hrf, 'TR', TR);
%   G    = hrf_granger_causality([xhat_region1, xhat_region2, ...]);
%
%   % Curve-shape misspecification vs a reference HRF (needs whole-brain score
%   % CSVs / an input table from the study pipeline):
%   Tmis = hrf_misspec_metrics(input_table, 'Reference', 'canonical');
%
%   % Whole-brain condition x lag t-map animation over the HRF:
%   hrf_animate_wholebrain_stats(results.wholebrain_by_model.sfir, ...
%       'Object', 't', 'Condition', 'pain', 'OutputFile', 'pain_tmaps.mp4');


%% 6. The takeaway
% One object model, one method vocabulary, extended to HRFs:
%
%   timeseries --hrf_fit_all_models-->  fitted curves (+ model SE)
%   wholebrain --make_fmri_stat_hrf-->  fmri_hrf  (IS-A fmri_data)
%                                       statistic_hrf (IS-A statistic_image)
%   statistic_hrf --threshold--> region --table/montage
%   study --hrf_time_unfolding_stats(Correction)--> corrected per-lag inference
%
% fmri_hrf and statistic_hrf are first-class citizens of the CANlab hierarchy:
% they inherit the entire fmri_data / statistic_image surface, add condition x
% lag HRF metadata and an HRF-aware display, and every HRF analysis lands back
% in a standard CANlab container -- so nothing you already know has to change.
disp('HRF object-oriented tour complete.');
