function out = hrf_fit_wholebrain_stats(fmri_nii, events_tsv, varargin)
%HRF_FIT_WHOLEBRAIN_STATS Vectorized whole-brain HRF beta and T maps.
%
% out = hrf_fit_wholebrain_stats(fmri_nii, events_tsv, ...)
%
% This fits a linear HRF model to every voxel in a 4D fMRI
% image and returns two CANlab statistic_image objects:
%   out.b  - beta/HRF amplitude maps, one 3D volume per condition x lag
%   out.t  - T maps for the same condition x lag volumes
%
% Both objects include .ste, .p, .dfe, .N, .image_labels, and .volInfo.
% If OutputPrefix is provided, 4D NIfTI files are written to disk.
%
% SPM GKWY compatibility
% ----------------------
% By default this is OLS on a constant-only baseline -- it is NOT identical
% to an SPM first-level GLM, which fits the grand-mean-scaled, high-pass
% filtered, prewhitened data (SPM's "gKWY"; see Misc_utilities/spmify.m).
% Two tiers make the fit SPM-comparable:
%
%   Tier B (exact): pass 'SPM', '/path/to/SPM.mat' (or an estimated SPM
%     struct). The data are transformed to gKWY (global scale g, high-pass K,
%     ReML whitening W) and the design is filtered to KWX -- reproducing
%     spm_spm. Use 'SPMRun' to pick the run for a multi-run SPM.
%
%   Tier A (no SPM.mat): 'HighpassSeconds' (default 128, SPM's default) adds
%     DCT high-pass confounds to the design, replicating K and removing the
%     low-frequency drift that otherwise leaks into long FIR/sFIR lags as a
%     spurious sustained baseline. ScaleMode 'grandmean' replicates g.
%     'Whiten' ('ar1' or 'ar2') estimates the noise autocorrelation from the
%     residuals of a high-variance voxel subsample (one global AR model,
%     SPM-style) and prewhitens data + design -- giving GLS-valid SE/t
%     without an SPM.mat. Set HighpassSeconds [] / 0 / Inf and Whiten 'none'
%     to recover the legacy raw-OLS behavior.
%
% ScaleMode: 'none' (default), 'zscore' (per-voxel), or 'grandmean' (single
% global scalar to mean 100, SPM-style). Ignored when 'SPM' is supplied.

p = inputParser;
p.addRequired('fmri_nii', @(x) ischar(x) || isstring(x) || isa(x, 'fmri_data'));
p.addRequired('events_tsv', @(x) ischar(x) || isstring(x) || istable(x));
p.addParameter('TR', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('MaskNii', '', @(x) ischar(x) || isstring(x) || isa(x, 'fmri_mask_image'));
p.addParameter('Conditions', {}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('WindowSeconds', 30, @(x) isnumeric(x) && all(x(:) > 0));
p.addParameter('Mode', 'FIR', @(x) ischar(x) || isstring(x));
p.addParameter('Nuisance', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('OutputPrefix', '', @(x) ischar(x) || isstring(x));
p.addParameter('Overwrite', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('PThresh', [], @(x) isempty(x) || (isscalar(x) && x > 0 && x < 1));
p.addParameter('ThreshType', 'unc', @(x) ischar(x) || isstring(x));
p.addParameter('WriteThresholdedT', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('ChunkSize', 50000, @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x, 1) == 0);
p.addParameter('ScaleMode', 'none', @(x) ischar(x) || isstring(x));
% --- SPM GKWY-compatibility controls -------------------------------------
% The fit is OLS on a constant-only baseline unless you ask for SPM-style
% conditioning. Two tiers (see help block):
%   Tier B  'SPM'              - path to / struct of an ESTIMATED SPM.mat.
%                                Applies exact g (global scale), K (high-pass)
%                                and W (ReML whitening), matching spm_spm.
%           'SPMRun'           - which run in a multi-run SPM the image is (1).
%   Tier A  'HighpassSeconds'  - DCT high-pass cutoff in seconds (default 128,
%                                SPM's default). Replicates K when no SPM.mat
%                                is available. Set [] / 0 / Inf to disable.
% Whitening (W) is only available via Tier B. With an SPM.mat, ScaleMode and
% HighpassSeconds are ignored (g and K come from the SPM).
p.addParameter('SPM', [], @(x) isempty(x) || isstruct(x) || ischar(x) || isstring(x));
p.addParameter('SPMRun', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x, 1) == 0);
p.addParameter('HighpassSeconds', 128, @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
% Prewhitening for the no-SPM path: 'none' (default), 'ar1' or 'ar2'. The
% autocorrelation is estimated from the OLS residuals of a high-variance
% voxel subsample (one global model, SPM-style) and applied to data + design,
% giving GLS-valid SE/t without an SPM.mat. SPM-exact whitening uses 'SPM'.
p.addParameter('Whiten', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(fmri_nii, events_tsv, varargin{:});
opts = p.Results;

if isa(fmri_nii, 'fmri_data')
    data_obj = fmri_nii;
    fmri_name = '';
    if isempty(opts.TR)
        if isprop(data_obj, 'image_metadata') && isfield(data_obj.image_metadata, 'TR_in_sec') && ...
                ~isnan(data_obj.image_metadata.TR_in_sec)
            TR = data_obj.image_metadata.TR_in_sec;
        else
            error('TR is required when fmri_nii is an fmri_data object without image_metadata.TR_in_sec.');
        end
    else
        TR = opts.TR;
    end
    n_tp = size(data_obj.dat, 2);
else
    fmri_name = char(fmri_nii);
    if ~exist(fmri_name, 'file'), error('fMRI file not found: %s', fmri_name); end
    [TR, n_tp] = local_get_tr_and_ntp(fmri_name, opts.TR);

    mask_input = opts.MaskNii;
    if isstring(mask_input), mask_input = char(mask_input); end

    if isempty(mask_input)
        data_obj = fmri_data(fmri_name, 'noverbose');
    else
        data_obj = fmri_data(fmri_name, mask_input, 'noverbose');
    end
end

if istable(events_tsv)
    E = events_tsv;
else
    events_name = char(events_tsv);
    if ~exist(events_name, 'file'), error('Events file not found: %s', events_name); end
    E = hrf_load_events_tsv(events_name);
end

if isempty(opts.Conditions)
    cond_names = unique(E.trial_type, 'stable');
else
    cond_names = cellstr(string(opts.Conditions));
end

[Runc, condition_groups] = hrf_build_stick_functions(E, cond_names, TR, n_tp);
cond_names = {condition_groups.label};
[X, design_info] = local_build_wholebrain_design(Runc, cond_names, TR, opts.WindowSeconds, opts.Mode, opts.Nuisance);

if size(X, 1) ~= size(data_obj.dat, 2)
    error('Design rows (%d) must match fMRI time points (%d).', size(X, 1), size(data_obj.dat, 2));
end

% --- SPM GKWY conditioning (Tier B exact, or Tier A high-pass) -----------
% Transforms data_obj.dat and X up front so the existing vectorized OLS loop
% below becomes the consistent (filtered / whitened) fit. dof_correction
% accounts for high-pass confounds that Tier B projects out of the data.
dof_correction = 0;
gkwy_info = struct('mode', 'none', 'highpass_seconds', [], 'n_highpass', 0, ...
    'whitened', false, 'grandmean_scaled', false);

if ~isempty(opts.SPM)
    SPM = local_load_spm(opts.SPM);
    [data_obj, X, dof_correction, n_hp] = local_apply_spm_gkwy(data_obj, X, SPM, opts.SPMRun);
    gkwy_info.mode = 'spm_gkwy';
    gkwy_info.n_highpass = n_hp;
    gkwy_info.whitened = true;
    gkwy_info.grandmean_scaled = true;
    if ~strcmpi(char(opts.ScaleMode), 'none')
        warning('hrf_fit_wholebrain_stats:ScaleModeIgnored', ...
            'ScaleMode=''%s'' ignored: global scaling comes from SPM.xGX.gSF.', char(opts.ScaleMode));
        opts.ScaleMode = 'none';
    end
    if ~strcmpi(char(opts.Whiten), 'none')
        warning('hrf_fit_wholebrain_stats:WhitenIgnored', ...
            'Whiten=''%s'' ignored: whitening comes from SPM.xX.W.', char(opts.Whiten));
    end
    if opts.Verbose
        fprintf('SPM GKWY applied (run %d): global scale + high-pass (%d confounds) + ReML whitening.\n', ...
            opts.SPMRun, n_hp);
    end
else
    hp = opts.HighpassSeconds;
    if ~isempty(hp) && isfinite(hp) && hp > 0
        [X, design_info, n_hp] = local_add_highpass(X, design_info, TR, hp);
        gkwy_info.mode = 'highpass';
        gkwy_info.highpass_seconds = hp;
        gkwy_info.n_highpass = n_hp;
        if opts.Verbose
            fprintf('High-pass filter: DCT, cutoff %.0f s, %d confound regressors added.\n', hp, n_hp);
        end
    end
end

% Grand-mean scaling (Tier A replicate of SPM's g) needs one global scalar
% over all in-mask voxels x time, computed before the chunked loop.
gm_scale = 1;
if strcmpi(char(opts.ScaleMode), 'grandmean')
    gm = mean(double(data_obj.dat(:)), 'omitnan');
    if gm == 0 || ~isfinite(gm)
        warning('hrf_fit_wholebrain_stats:GrandMeanZero', ...
            'Grand mean is zero/non-finite; skipping grand-mean scaling.');
    else
        gm_scale = 100 / gm;
        gkwy_info.grandmean_scaled = true;
    end
end

% Data-estimated prewhitening (Tier A+). Only when no SPM.mat supplied; with
% an SPM the whitening already came from W above. Estimated on a high-variance
% voxel subsample, then applied to data + design so the loop below is GLS.
if isempty(opts.SPM) && ~strcmpi(strtrim(char(opts.Whiten)), 'none')
    [Wmat, ar_coef] = local_estimate_ar_whitening(data_obj.dat, X, char(opts.Whiten));
    data_obj.dat = single((Wmat * double(data_obj.dat'))');   % whiten data (vox x time)
    X = Wmat * X;                                             % whiten design
    gkwy_info.whitened = true;
    gkwy_info.whiten_mode = lower(strtrim(char(opts.Whiten)));
    gkwy_info.ar_coef = ar_coef;
    if opts.Verbose
        fprintf('Prewhitening: %s, AR coef = [%s].\n', upper(strtrim(char(opts.Whiten))), ...
            strtrim(sprintf('%.3f ', ar_coef(:)')));
    end
end

[PX, pen] = local_get_pseudoinverse(X, design_info, opts.Mode);
hat_trace = trace(X * PX);
dfe = max(size(X, 1) - hat_trace - dof_correction, 1);
coef_var_scale = max(diag(design_info.output_lift * (PX * PX') * design_info.output_lift'), 0);

n_keep = size(design_info.output_lift, 1);
n_vox = size(data_obj.dat, 1);

beta_dat = zeros(n_vox, n_keep, 'single');
ste_dat = zeros(n_vox, n_keep, 'single');
t_dat = zeros(n_vox, n_keep, 'single');
p_dat = zeros(n_vox, n_keep, 'single');

chunk_size = min(opts.ChunkSize, n_vox);
if opts.Verbose
    fprintf('Fitting whole-brain %s model: %d voxels, %d time points, %d output maps\n', ...
        upper(char(opts.Mode)), n_vox, size(X, 1), n_keep);
end

for first_vox = 1:chunk_size:n_vox
    last_vox = min(first_vox + chunk_size - 1, n_vox);
    wh = first_vox:last_vox;
    Y = double(data_obj.dat(wh, :)');

    switch lower(char(opts.ScaleMode))
        case 'none'
            % leave data in native units
        case {'zscore', 'z'}
            Y = zscore(Y, 0, 1);
            Y(isnan(Y)) = 0;
        case 'grandmean'
            Y = Y * gm_scale;   % single global scalar (SPM-style grand-mean to 100)
        otherwise
            error('Unknown ScaleMode: %s. Use ''none'', ''zscore'', or ''grandmean''.', char(opts.ScaleMode));
    end

    B = PX * Y;
    R = Y - X * B;
    mse = sum(R .^ 2, 1) ./ dfe;
    SE = sqrt(coef_var_scale * mse);

    Bkeep = design_info.output_lift * B;
    SEkeep = SE;
    Tkeep = Bkeep ./ SEkeep;
    Pkeep = 2 * (1 - tcdf(abs(Tkeep), dfe));
    Pkeep(Pkeep == 0) = eps;

    beta_dat(wh, :) = single(Bkeep');
    ste_dat(wh, :) = single(SEkeep');
    t_dat(wh, :) = single(Tkeep');
    p_dat(wh, :) = single(Pkeep');

    if opts.Verbose
        fprintf('  voxels %d-%d / %d\n', first_vox, last_vox, n_vox);
    end
end

labels = design_info.labels;
meta_table = design_info.metadata_table;
meta_table.N = repmat(size(X, 1), height(meta_table), 1);
meta_table.dfe = repmat(dfe, height(meta_table), 1);
meta_table.TR = repmat(TR, height(meta_table), 1);
meta_table.mode = repmat(string(upper(char(opts.Mode))), height(meta_table), 1);

b_obj = statistic_image;
b_obj.type = sprintf('%s HRF beta', upper(char(opts.Mode)));
b_obj.dat = beta_dat;
b_obj.p = p_dat;
b_obj.p_type = sprintf('Two-tailed P-values from beta/SE, dfe = %.3f', dfe);
b_obj.ste = ste_dat;
b_obj.sig = true(size(beta_dat));
b_obj.N = size(X, 1);
b_obj.dfe = dfe;
b_obj.volInfo = data_obj.volInfo;
b_obj.removed_voxels = data_obj.removed_voxels;
b_obj.removed_images = false(1, n_keep);
b_obj.image_labels = labels;
b_obj.dat_descrip = sprintf('Whole-brain %s HRF beta values: condition x lag maps', upper(char(opts.Mode)));

t_obj = statistic_image;
t_obj.type = 'T';
t_obj.dat = t_dat;
t_obj.p = p_dat;
t_obj.p_type = sprintf('Two-tailed P-values from beta/SE, dfe = %.3f', dfe);
t_obj.ste = ste_dat;
t_obj.sig = true(size(t_dat));
t_obj.N = size(X, 1);
t_obj.dfe = dfe;
t_obj.volInfo = data_obj.volInfo;
t_obj.removed_voxels = data_obj.removed_voxels;
t_obj.removed_images = false(1, n_keep);
t_obj.image_labels = labels;
t_obj.dat_descrip = sprintf('Whole-brain %s HRF T values: condition x lag maps', upper(char(opts.Mode)));

if ~isempty(opts.PThresh)
    t_obj = threshold(t_obj, opts.PThresh, char(opts.ThreshType), 'noverbose');
end

out = struct();
out.b = b_obj;
out.t = t_obj;
out.design_matrix = X;
out.design_info = design_info;
out.design_info.penalty = pen;
out.metadata_table = meta_table;
out.conditions = cond_names;
out.condition_groups = condition_groups;
out.TR = TR;
out.dfe = dfe;
out.N = size(X, 1);
out.gkwy = gkwy_info;
out.input_fmri = fmri_name;
out.paths = struct();

if ~isempty(opts.OutputPrefix)
    out.paths = local_write_outputs(out, char(opts.OutputPrefix), logical(opts.Overwrite), logical(opts.WriteThresholdedT));
end
end

function SPM = local_load_spm(spm_in)
% Resolve the SPM arg to a struct (accepts a struct or a path to SPM.mat).
if isstruct(spm_in)
    SPM = spm_in;
elseif ischar(spm_in) || isstring(spm_in)
    f = char(spm_in);
    if exist(f, 'file') ~= 2
        error('hrf_fit_wholebrain_stats:SPMNotFound', 'SPM.mat not found: %s', f);
    end
    S = load(f, 'SPM');
    if ~isfield(S, 'SPM')
        error('hrf_fit_wholebrain_stats:NoSPMVar', '%s does not contain an SPM variable.', f);
    end
    SPM = S.SPM;
else
    error('hrf_fit_wholebrain_stats:BadSPM', 'SPM must be a struct or a path to SPM.mat.');
end
end

function [data_obj, X, dof_corr, n_hp] = local_apply_spm_gkwy(data_obj, X, SPM, run)
% Replicate spmify()/spm_spm: scale data by g (SPM.xGX.gSF), high-pass (K) and
% whiten (W) both data and design. Handles two shapes of fMRI input against a
% (possibly multi-run, e.g. per-session-concatenated) SPM:
%   (a) the WHOLE concatenated session -- n_tp == total SPM scans: apply every
%       run's K block, the full (block-diagonal) W, and per-scan gSF.
%   (b) a SINGLE run of a multi-run SPM -- n_tp == that run's scan count: select
%       SPMRun and remap K(run).row (concatenated indices) to local 1:n_tp.
% Column layout of X is preserved, so design_info.output_lift / sFIR penalty
% stay valid.
if ~isfield(SPM, 'xX') || ~isfield(SPM.xX, 'K') || ~isfield(SPM.xX, 'W')
    error('hrf_fit_wholebrain_stats:SPMNotEstimated', ...
        'SPM.mat must be ESTIMATED (needs xX.K, xX.W, xGX.gSF). Run spm_spm first.');
end
K = SPM.xX.K;
nrun = numel(K);
n_tp = size(data_obj.dat, 2);
total_scans = size(SPM.xX.W, 1);

if n_tp == total_scans
    % (a) Whole concatenated session.
    Kapply = K;                 % spm_filter loops over each K(s).row block
    rows = 1:total_scans;
    n_hp = local_sum_x0_cols(K);
else
    % (b) One run from the SPM; remap its rows to local frames.
    if run < 1 || run > nrun
        error('hrf_fit_wholebrain_stats:BadSPMRun', ...
            'SPMRun=%d out of range (SPM has %d run(s)).', run, nrun);
    end
    rows = K(run).row;
    if numel(rows) ~= n_tp
        error('hrf_fit_wholebrain_stats:SPMRunMismatch', ...
            ['fMRI image has %d time points, but neither the full SPM (%d scans) ' ...
             'nor SPM run %d (%d scans) matches. For a per-session SPM with ' ...
             'concatenated runs, pass the concatenated-session fMRI (matches the ' ...
             'full scan count), or set SPMRun to this run''s index within the session.'], ...
            n_tp, total_scans, run, numel(rows));
    end
    Kapply = K(run);
    Kapply.row = 1:n_tp;        % concatenated indices -> local frames
    n_hp = local_sum_x0_cols(Kapply);
end

W = SPM.xX.W(rows, rows);

% Whiten + high-pass the data. spm_filter works on [time x columns].
Y = double(data_obj.dat');      % time x vox
KWY = spm_filter(Kapply, W * Y);

% Global (grand-mean) scaling per scan, exactly as spmify does.
if isfield(SPM, 'xGX') && isfield(SPM.xGX, 'gSF') && ~isempty(SPM.xGX.gSF)
    g = SPM.xGX.gSF(rows);
    KWY = KWY .* g(:);
else
    warning('hrf_fit_wholebrain_stats:NoGSF', ...
        'SPM.xGX.gSF missing; skipped global scaling (K and W still applied).');
end
data_obj.dat = single(KWY');    % vox x time

% Apply the SAME K and W to the design.
X = full(spm_filter(Kapply, W * X));

dof_corr = n_hp;   % high-pass confounds are projected out, not columns
end

function n = local_sum_x0_cols(K)
% Total high-pass confound regressors across all K blocks.
n = 0;
for s = 1:numel(K)
    if isfield(K(s), 'X0') && ~isempty(K(s).X0)
        n = n + size(K(s).X0, 2);
    end
end
end

function [X, info, n_added] = local_add_highpass(X, info, TR, cutoff)
% Tier A: append SPM-style DCT high-pass confounds as nuisance columns.
% This removes low-frequency drift from the task betas (equivalent to K) and
% the added columns are counted in the OLS dof automatically.
D = local_dct_highpass(size(X, 1), TR, cutoff);
n_added = size(D, 2);
if n_added == 0, return; end
X = [X, D];
nkeep = size(info.output_lift, 1);
info.output_lift = [info.output_lift, zeros(nkeep, n_added)];
info.highpass_columns = (size(X, 2) - n_added + 1):size(X, 2);
info.highpass_seconds = cutoff;
end

function [W, ar] = local_estimate_ar_whitening(dat, X, wmode)
% Estimate a single global AR(p) whitening matrix from OLS residuals of a
% high-variance voxel subsample, then return W such that W*Y is ~white.
switch lower(strtrim(wmode))
    case 'ar1', p = 1;
    case 'ar2', p = 2;
    otherwise
        error('hrf_fit_wholebrain_stats:UnknownWhiten', ...
            'Whiten must be ''none'', ''ar1'', or ''ar2''. Got ''%s''.', wmode);
end
n_tp = size(dat, 2);

% Subsample voxels (highest variance, finite) to estimate AR cheaply.
v = var(double(dat), 0, 2);
v(~isfinite(v)) = 0;
n_pos = sum(v > 0);
nsamp = min(2000, max(n_pos, 1));
[~, ord] = sort(v, 'descend');
samp = ord(1:nsamp);
Ys = double(dat(samp, :)');            % time x nsamp
R = Ys - X * (pinv(X) * Ys);           % OLS residuals on the (high-passed) design

% Pooled normalized autocovariance across sampled voxels (lags 0..p).
acf = zeros(p + 1, 1);
for k = 0:p
    rr = R(1:end - k, :) .* R(1 + k:end, :);
    acf(k + 1) = mean(rr(:));
end
if acf(1) <= 0
    W = speye(n_tp); ar = zeros(p, 1); return;
end
acf = acf / acf(1);

% Yule-Walker for AR(p), then enforce stationarity.
ar = toeplitz(acf(1:p)) \ acf(2:p + 1);
ar = local_make_stationary(ar);

% Theoretical AR(p) autocorrelation over all lags, then whiten via inv(chol).
full_acf = local_ar_acf(ar, n_tp);
V = toeplitz(full_acf);
L = chol(V + 1e-8 * eye(n_tp), 'lower');
W = inv(L);
end

function ar = local_make_stationary(ar)
% Shrink AR coefficients toward zero until the process is stationary.
for it = 1:50
    if isscalar(ar)
        ok = abs(ar(1)) < 0.999;
    else
        ok = abs(ar(2)) < 0.999 && (ar(1) + ar(2)) < 0.999 && (ar(2) - ar(1)) < 0.999;
    end
    if ok, return; end
    ar = ar * 0.95;
end
ar = zeros(size(ar));
end

function a = local_ar_acf(ar, n)
% Theoretical autocorrelation of an AR(p) process, lags 0..n-1 (a(1)=lag 0).
a = zeros(n, 1); a(1) = 1;
if isscalar(ar)
    rho = max(min(ar(1), 0.999), -0.999);
    a = (rho .^ (0:n - 1))';
    return;
end
ar1 = ar(1); ar2 = ar(2);
if n >= 2, a(2) = ar1 / (1 - ar2); end
for k = 3:n
    a(k) = ar1 * a(k - 1) + ar2 * a(k - 2);
end
end

function D = local_dct_highpass(n_tp, TR, cutoff)
% SPM's high-pass basis: spm_dctmtx(N, k) with k = fix(2*N*TR/cutoff + 1),
% dropping the DC (constant) column since X already carries an intercept.
k = fix(2 * (n_tp * TR) / cutoff + 1);
if k <= 1
    D = zeros(n_tp, 0);
    return;
end
D = spm_dctmtx(n_tp, k);
D = D(:, 2:end);
end

function [TR, n_tp] = local_get_tr_and_ntp(fmri_name, requested_tr)
info = niftiinfo(fmri_name);
sz = info.ImageSize;
if numel(sz) < 4
    error('Expected a 4D fMRI image, got %d dimensions.', numel(sz));
end
n_tp = sz(4);

TR = requested_tr;
if isempty(TR) && isfield(info, 'PixelDimensions') && numel(info.PixelDimensions) >= 4
    TR = info.PixelDimensions(4);
end
if isempty(TR) || TR <= 0
    error('Could not infer TR from NIfTI header. Pass ''TR'' explicitly.');
end
end

function [X, info] = local_build_wholebrain_design(Runc, cond_names, TR, window_seconds, mode, nuisance)
switch lower(char(mode))
    case {'fir', 'sfir'}
        [X, info] = local_build_fir_design(Runc, cond_names, TR, window_seconds, mode, nuisance);
    case 'canonical'
        [X, info] = local_build_canonical_design(Runc, cond_names, TR, window_seconds, mode, nuisance);
    case 'spline'
        [X, info] = local_build_spline_design(Runc, cond_names, TR, window_seconds, mode, nuisance);
    otherwise
        error('Unknown Mode: %s. Use ''FIR'', ''sFIR'', ''canonical'', or ''spline''.', char(mode));
end
end

function [X, info] = local_build_fir_design(Runc, cond_names, TR, window_seconds, mode, nuisance)
numstim = numel(Runc);
len = numel(Runc{1});
if isscalar(window_seconds)
    window_seconds = repmat(window_seconds, 1, numstim);
end
if numel(window_seconds) ~= numstim
    error('WindowSeconds must be scalar or one value per condition.');
end

Runs = zeros(len, numstim);
for i = 1:numstim
    Runs(:, i) = Runc{i}(:);
end

DX_all = cell(1, numstim);
tlen_all = zeros(1, numstim);
fir_columns_by_condition = cell(1, numstim);
labels = {};
condition_col = {};
condition_index = [];
lag_index = [];
lag_seconds = [];

start_col = 1;
for i = 1:numstim
    t = 1:TR:window_seconds(i);
    tlen_all(i) = numel(t);
    DX_i = tor_make_deconv_mtx3(Runs(:, i), tlen_all(i), 1);
    DX_all{i} = DX_i(:, 1:tlen_all(i));

    fir_columns_by_condition{i} = start_col:(start_col + tlen_all(i) - 1);
    for j = 1:tlen_all(i)
        label = sprintf('%s_lag%03d_%0.3gs', local_safe_label(cond_names{i}), j, t(j));
        labels{end + 1, 1} = label; %#ok<AGROW>
        condition_col{end + 1, 1} = cond_names{i}; %#ok<AGROW>
        condition_index(end + 1, 1) = i; %#ok<AGROW>
        lag_index(end + 1, 1) = j; %#ok<AGROW>
        lag_seconds(end + 1, 1) = t(j); %#ok<AGROW>
    end
    start_col = start_col + tlen_all(i);
end

DX = horzcat(DX_all{:});
X = [DX ones(len, 1)];

if ~isempty(nuisance)
    if size(nuisance, 1) ~= len
        error('Nuisance matrix must have %d rows to match fMRI time points.', len);
    end
    X = [X nuisance];
end

info = struct();
info.mode = char(mode);
info.TR = TR;
info.tlen = tlen_all;
info.fir_columns = 1:sum(tlen_all);
info.fir_columns_by_condition = fir_columns_by_condition;
info.intercept_column = sum(tlen_all) + 1;
info.nuisance_columns = (info.intercept_column + 1):size(X, 2);
info.labels = labels;
info.metadata_table = table((1:numel(labels))', condition_col, condition_index, lag_index, lag_seconds, labels, ...
    'VariableNames', {'volume_index', 'condition', 'condition_index', 'lag_index', 'lag_seconds', 'image_label'});
info.output_lift = zeros(numel(labels), size(X, 2));
for i = 1:numel(labels)
    info.output_lift(i, i) = 1;
end
end

function [X, info] = local_build_canonical_design(Runc, cond_names, TR, window_seconds, mode, nuisance)
[~, t, tlen] = local_scalar_window(window_seconds, TR, mode);
numstim = numel(Runc);
len = numel(Runc{1});
h = local_canonical_basis(TR, tlen);

Xtask = zeros(len, numstim);
for i = 1:numstim
    v = conv(Runc{i}(:), h(:));
    Xtask(:, i) = v(1:len);
end
X = [ones(len, 1), Xtask];
if ~isempty(nuisance)
    if size(nuisance, 1) ~= len
        error('Nuisance matrix must have %d rows to match fMRI time points.', len);
    end
    X = [X nuisance];
end

[labels, condition_col, condition_index, lag_index, lag_seconds] = ...
    local_condition_lag_metadata(cond_names, t, mode);

info = struct();
info.mode = char(mode);
info.TR = TR;
info.tlen = repmat(tlen, 1, numstim);
info.intercept_column = 1;
info.nuisance_columns = (numstim + 2):size(X, 2);
info.labels = labels;
info.metadata_table = table((1:numel(labels))', condition_col, condition_index, lag_index, lag_seconds, labels, ...
    'VariableNames', {'volume_index', 'condition', 'condition_index', 'lag_index', 'lag_seconds', 'image_label'});
info.output_lift = zeros(numel(labels), size(X, 2));
for c = 1:numstim
    col = c + 1;
    rows = ((c - 1) * tlen + 1):(c * tlen);
    info.output_lift(rows, col) = h(:);
end
end

function [X, info] = local_build_spline_design(Runc, cond_names, TR, window_seconds, mode, nuisance)
[window_seconds, t, tlen] = local_scalar_window(window_seconds, TR, mode); %#ok<ASGLU>
numstim = numel(Runc);
len = numel(Runc{1});
K = 8;
norder = 4;

try
    basis = create_bspline_basis([0, tlen], K + 3, norder);
    Bbasis = eval_basis((1:tlen), basis);
catch
    error('Mode=''spline'' requires the FDA package on the MATLAB path (missing create_bspline_basis/eval_basis; see https://github.com/markgewhite/fda).');
end
Bbasis = Bbasis(:, 3:end-1);

Wi = zeros(len, numstim * K);
for c = 1:numstim
    Wc = tor_make_deconv_mtx3(Runc{c}(:), tlen, 1);
    Wi(:, (c - 1) * K + 1:c * K) = Wc(:, 1:tlen) * Bbasis;
end
X = [ones(len, 1), Wi];
if ~isempty(nuisance)
    if size(nuisance, 1) ~= len
        error('Nuisance matrix must have %d rows to match fMRI time points.', len);
    end
    X = [X nuisance];
end

[labels, condition_col, condition_index, lag_index, lag_seconds] = ...
    local_condition_lag_metadata(cond_names, t, mode);

info = struct();
info.mode = char(mode);
info.TR = TR;
info.tlen = repmat(tlen, 1, numstim);
info.intercept_column = 1;
info.nuisance_columns = (numstim * K + 2):size(X, 2);
info.labels = labels;
info.metadata_table = table((1:numel(labels))', condition_col, condition_index, lag_index, lag_seconds, labels, ...
    'VariableNames', {'volume_index', 'condition', 'condition_index', 'lag_index', 'lag_seconds', 'image_label'});
info.output_lift = zeros(numel(labels), size(X, 2));
for c = 1:numstim
    rows = ((c - 1) * tlen + 1):(c * tlen);
    cols = (c - 1) * K + 2:c * K + 1;
    info.output_lift(rows, cols) = Bbasis;
end
end

function [window_seconds, t, tlen] = local_scalar_window(window_seconds, TR, mode)
if ~isscalar(window_seconds)
    if all(window_seconds == window_seconds(1))
        window_seconds = window_seconds(1);
    else
        error('WindowSeconds must be scalar for whole-brain Mode=''%s''.', char(mode));
    end
end
t = 1:TR:window_seconds;
tlen = numel(t);
end

function h = local_canonical_basis(TR, tlen)
len = max(round(30 / TR), tlen);
xBF.dt = TR;
xBF.length = len;
xBF.name = 'hrf (with time and dispersion derivatives)';
xBF = spm_get_bf(xBF);
h = xBF.bf(1:tlen, 1);
h = h ./ max(h);
end

function [labels, condition_col, condition_index, lag_index, lag_seconds] = local_condition_lag_metadata(cond_names, t, mode)
labels = {};
condition_col = {};
condition_index = [];
lag_index = [];
lag_seconds = [];
for c = 1:numel(cond_names)
    for j = 1:numel(t)
        label = sprintf('%s_%s_lag%03d_%0.3gs', local_safe_label(mode), local_safe_label(cond_names{c}), j, t(j));
        labels{end + 1, 1} = label; %#ok<AGROW>
        condition_col{end + 1, 1} = cond_names{c}; %#ok<AGROW>
        condition_index(end + 1, 1) = c; %#ok<AGROW>
        lag_index(end + 1, 1) = j; %#ok<AGROW>
        lag_seconds(end + 1, 1) = t(j); %#ok<AGROW>
    end
end
end

function [PX, pen] = local_get_pseudoinverse(X, info, mode)
mode = lower(char(mode));
pen = zeros(size(X, 2));

switch mode
    case {'fir', 'canonical', 'spline'}
        PX = pinv(X);
    case 'sfir'
        start_idx = 1;
        for i = 1:numel(info.tlen)
            tlen = info.tlen(i);
            C = (1:tlen)' * ones(1, tlen);
            h = sqrt(1 / (7 / info.TR));
            v = 0.1;
            sig = 1;
            R = v * exp(-h / 2 * (C - C') .^ 2);
            RI = R \ eye(size(R));
            end_idx = start_idx + tlen - 1;
            pen(start_idx:end_idx, start_idx:end_idx) = sig ^ 2 * RI;
            start_idx = end_idx + 1;
        end
        PX = (X' * X + pen) \ X';
    otherwise
        error('Unknown Mode: %s. Use ''FIR'', ''sFIR'', ''canonical'', or ''spline''.', mode);
end
end

function paths = local_write_outputs(out, prefix, overwrite, write_thresholded_t)
paths = struct();
paths.beta = [prefix '_beta.nii'];
paths.t = [prefix '_t.nii'];
paths.se = [prefix '_se.nii'];
paths.p = [prefix '_p.nii'];
paths.metadata = [prefix '_metadata.csv'];

local_write_image(out.b, paths.beta, overwrite);

local_write_image(out.t, paths.t, overwrite);

se_obj = out.b;
se_obj.type = 'HRF beta standard error';
se_obj.dat = out.b.ste;
se_obj.dat_descrip = 'Whole-brain HRF beta standard error maps';
local_write_image(se_obj, paths.se, overwrite);

p_obj = out.t;
p_obj.type = 'p';
p_obj.dat = out.t.p;
p_obj.p = out.t.p;
p_obj.dat_descrip = 'Whole-brain HRF two-tailed p-value maps';
local_write_image(p_obj, paths.p, overwrite);

local_delete_if_overwrite(paths.metadata, overwrite);
writetable(out.metadata_table, paths.metadata);

if write_thresholded_t
    paths.t_thresholded = [prefix '_t_thresh.nii'];
    local_write_image(out.t, paths.t_thresholded, overwrite, 'thresh');
end
end

function local_write_image(obj, fname, overwrite, varargin)
local_delete_if_overwrite(fname, overwrite);
write_args = [{'fname', fname}, varargin(:)'];
if overwrite
    write_args{end + 1} = 'overwrite';
end
obj = local_prepare_image_for_write(obj);
write(obj, write_args{:});
end

function obj = local_prepare_image_for_write(obj)
if ~isprop(obj, 'removed_images') || isempty(obj.removed_images) || isempty(obj.dat)
    return
end

removed = logical(obj.removed_images(:)');
n_images = size(obj.dat, 2);
if numel(removed) ~= n_images
    return
end

if any(removed)
    % threshold() can mark 4D maps as removed even when .dat still contains
    % the full volume series.  image_vector/write() reinserts removed images,
    % so reset this flag when the image count is already complete.
    obj.removed_images = false(size(obj.removed_images));
end
end

function local_delete_if_overwrite(fname, overwrite)
if exist(fname, 'file') ~= 2
    return
end
if ~overwrite
    error('Output file already exists: %s. Use Overwrite=true to replace it.', fname);
end

try
    fileattrib(fname, '+w');
catch
end
try
    delete(fname);
catch err
    error('Could not overwrite existing output file %s: %s', fname, err.message);
end
end

function s = local_safe_label(s)
s = char(s);
s = matlab.lang.makeValidName(s);
end
