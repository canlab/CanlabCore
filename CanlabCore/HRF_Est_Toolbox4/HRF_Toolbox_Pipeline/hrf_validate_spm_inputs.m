function report = hrf_validate_spm_inputs(fmri_files, spm_files, varargin)
% Pre-flight check that per-subject SPM.mat files line up with fMRI inputs
% for exact Tier B GKWY conditioning (see hrf_fit_wholebrain_stats).
%
% For each fMRI file / SPM.mat pair it confirms the SPM is ESTIMATED (carries
% the whitening xX.W, high-pass xX.K, and global scaling xGX.gSF) and that the
% fMRI time-point count matches the SPM in one of two ways:
%   - 'full-session'  : n_tp == total SPM scans (sum over concatenated runs)
%   - 'single-run'    : n_tp == the SPMRun-th run's scan count
% Anything else is flagged so it fails here, on the login node, rather than
% mid-array on the cluster.
%
% :Usage:
% ::
%     report = hrf_validate_spm_inputs(fmri_files, spm_files, ['SPMRuns', runs])
%
% :Inputs:
%
%   **fmri_files:**
%        cellstr / string array of 4D fMRI paths (one per task), .nii or .nii.gz.
%
%   **spm_files:**
%        cellstr / string array of estimated SPM.mat paths, one per fMRI file.
%        An empty entry ('' ) means "no SPM for this file" (Tier A fallback) and
%        is reported as 'no-spm' (not an error).
%
% :Optional Inputs:
%
%   **'SPMRuns':**
%        numeric vector, one run index per fMRI file (default all 1). Only used
%        when the fMRI is a single run of a multi-run (e.g. per-session) SPM.
%
%   **'Verbose' / 'doverbose':**
%        logical, print a per-row table and summary (default true).
%
% :Outputs:
%
%   **report:**
%        table with one row per input: index, fmri_file, spm_file, spm_run,
%        n_tp, total_scans, run_scans, mode, status, ok, message. `ok` is true
%        only for 'full-session' or 'single-run'.
%
% :Examples:
% ::
%     % Build per-subject SPM paths, then validate before submitting SLURM.
%     spm_files = cellfun(@(s) fullfile(spm_root, s, 'SPM.mat'), subject_ids, 'unif', 0);
%     report = hrf_validate_spm_inputs(fmri_files, spm_files);
%     assert(all(report.ok), 'Fix the flagged SPM/fMRI mismatches before launching.');
%
% See also: hrf_fit_wholebrain_stats, hrf_write_slurm_study_script, spmify.

% ----------------------------- parse inputs ------------------------------
p = inputParser;
p.addRequired('fmri_files', @(x) ischar(x) || iscell(x) || isstring(x));
p.addRequired('spm_files',  @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('SPMRuns', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(fmri_files, spm_files, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);
if ~isempty(opts.doverbose), verbose = logical(opts.doverbose); end

fmri_files = local_to_cellstr(fmri_files);
spm_files  = local_to_cellstr(spm_files);
n = numel(fmri_files);
if numel(spm_files) ~= n
    error('hrf_validate_spm_inputs:LengthMismatch', ...
        'fmri_files (%d) and spm_files (%d) must have the same length.', n, numel(spm_files));
end
spm_runs = opts.SPMRuns;
if isempty(spm_runs), spm_runs = ones(1, n); elseif isscalar(spm_runs), spm_runs = repmat(spm_runs, 1, n); end
if numel(spm_runs) ~= n
    error('hrf_validate_spm_inputs:SPMRunsLength', 'SPMRuns must be empty, scalar, or length %d.', n);
end

% ----------------------------- per-row check -----------------------------
idx = (1:n)';
fmri_col = fmri_files(:);
spm_col = spm_files(:);
spm_run_col = spm_runs(:);
n_tp = nan(n, 1); total_scans = nan(n, 1); run_scans = nan(n, 1);
mode = strings(n, 1); status = strings(n, 1); ok = false(n, 1); message = strings(n, 1);

for i = 1:n
    [n_tp(i), total_scans(i), run_scans(i), mode(i), status(i), ok(i), message(i)] = ...
        local_check_one(fmri_files{i}, spm_files{i}, spm_runs(i));
end

report = table(idx, fmri_col, spm_col, spm_run_col, n_tp, total_scans, run_scans, mode, status, ok, message, ...
    'VariableNames', {'index', 'fmri_file', 'spm_file', 'spm_run', 'n_tp', ...
    'total_scans', 'run_scans', 'mode', 'status', 'ok', 'message'});

if verbose
    disp(report(:, {'index', 'spm_run', 'n_tp', 'total_scans', 'run_scans', 'status'}));
    n_ok = sum(ok); n_nospm = sum(status == "no-spm"); n_bad = n - n_ok - n_nospm;
    fprintf('hrf_validate_spm_inputs: %d/%d ok (%d Tier B), %d no-spm (Tier A), %d PROBLEM.\n', ...
        n_ok, n, n_ok, n_nospm, n_bad);
    if n_bad > 0
        bad = report(~ok & status ~= "no-spm", :);
        for b = 1:height(bad)
            fprintf('  [%d] %s -> %s: %s\n', bad.index(b), local_short(bad.fmri_file{b}), ...
                bad.status(b), bad.message(b));
        end
    end
end
end


% =========================================================================
function [n_tp, total_scans, run_scans, mode, status, ok, msg] = local_check_one(fmri_file, spm_file, run)
n_tp = NaN; total_scans = NaN; run_scans = NaN;
mode = ""; status = ""; ok = false; msg = "";

% --- fMRI time points ---
try
    n_tp = local_ntp(fmri_file);
catch err
    status = "fmri-unreadable"; msg = string(err.message); return
end

% --- SPM present? ---
if isempty(strtrim(char(spm_file)))
    status = "no-spm"; msg = "no SPM supplied (Tier A high-pass fallback)"; return
end
if exist(char(spm_file), 'file') ~= 2
    status = "spm-missing"; msg = "SPM.mat path does not exist"; return
end

% --- SPM estimated? ---
try
    S = load(char(spm_file), 'SPM');
catch err
    status = "spm-unloadable"; msg = string(err.message); return
end
if ~isfield(S, 'SPM')
    status = "spm-no-variable"; msg = "file has no SPM variable"; return
end
SPM = S.SPM;
if ~isfield(SPM, 'xX') || ~isfield(SPM.xX, 'K') || ~isfield(SPM.xX, 'W')
    status = "spm-not-estimated"; msg = "missing xX.K / xX.W (run spm_spm first)"; return
end
if ~isfield(SPM, 'xGX') || ~isfield(SPM.xGX, 'gSF') || isempty(SPM.xGX.gSF)
    msg = "xGX.gSF missing -- K and W will apply, global scaling skipped";
end

% --- scan-count match ---
K = SPM.xX.K;
total_scans = size(SPM.xX.W, 1);
if run >= 1 && run <= numel(K)
    run_scans = numel(K(run).row);
end

if n_tp == total_scans
    mode = "full-session"; ok = true;
    if msg == "", msg = sprintf("matches full SPM (%d scans across %d run(s))", total_scans, numel(K)); end
elseif ~isnan(run_scans) && n_tp == run_scans
    mode = "single-run"; ok = true;
    if msg == "", msg = sprintf("matches SPM run %d (%d scans)", run, run_scans); end
else
    status = "scan-count-mismatch";
    msg = sprintf("fMRI=%d tp; full SPM=%d; run %d=%s. Set SPMRuns or pass the concatenated-session fMRI.", ...
        n_tp, total_scans, run, local_num2str(run_scans));
    return
end
status = mode;
end

function n_tp = local_ntp(fmri_file)
info = niftiinfo(char(fmri_file));
sz = info.ImageSize;
if numel(sz) < 4
    error('expected a 4D image, got %d dims', numel(sz));
end
n_tp = sz(4);
end

function c = local_to_cellstr(x)
if ischar(x)
    c = {x};
elseif isstring(x)
    c = cellstr(x);
elseif iscell(x)
    c = cellfun(@(v) char(string(v)), x, 'UniformOutput', false);
else
    error('hrf_validate_spm_inputs:BadType', 'Expected char, string, or cell of paths.');
end
c = c(:)';
end

function s = local_short(p)
[~, nm, ext] = fileparts(char(p));
s = [nm ext];
end

function s = local_num2str(x)
if isnan(x), s = 'n/a'; else, s = num2str(x); end
end
