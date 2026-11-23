function [spm_runs, report] = hrf_resolve_spm_runs(fmri_files, spm_files, varargin)
% Resolve, for each single-run fMRI file, its run index within a per-session
% (concatenated-runs) SPM.mat -- the value to pass as 'SPMRuns' to
% hrf_fit_wholebrain_stats / hrf_write_slurm_study_script for Tier B GKWY.
%
% Getting this index right matters: if a session's runs share the same scan
% count, a guessed index would pass the scan-count check yet apply the WRONG
% run's high-pass/whitening block. This reads each SPM's actual per-run source
% images and matches each fMRI file to its block, preferring the BIDS run token.
%
% Match order (most to least specific), per file:
%   1. 'run-XX' token shared between the fMRI file and one SPM run's source
%      images (robust to preprocessing prefixes/suffixes).
%   2. fMRI basename contained in one SPM run's source image names.
%   3. unique scan-count match (only one SPM run has numel(rows) == n_tp).
%   Otherwise -> NaN, flagged 'ambiguous' / 'nomatch' in the report.
%
% :Usage:
% ::
%     [spm_runs, report] = hrf_resolve_spm_runs(fmri_files, spm_files)
%
% :Inputs:
%
%   **fmri_files:**
%        cellstr/string of single-run 4D fMRI paths (.nii/.nii.gz).
%
%   **spm_files:**
%        cellstr/string of the per-session SPM.mat path for each fMRI file
%        (the SPM that concatenated that file's session). '' entries -> run 1.
%
% :Optional Inputs:
%
%   **'Verbose' / 'doverbose':**
%        logical, print the resolution table (default true).
%
% :Outputs:
%
%   **spm_runs:**
%        1 x N numeric, the resolved run index per fMRI file (NaN where it
%        could not be resolved -- inspect `report` and resolve those by hand).
%
%   **report:**
%        table: index, fmri_file, spm_file, n_tp, n_runs, spm_run, method,
%        message.
%
% :Examples:
% ::
%     [spm_runs, report] = hrf_resolve_spm_runs(fmri_files, spm_files);
%     assert(~any(isnan(spm_runs)), 'resolve the NaN rows before launching');
%     % then: hrf_write_slurm_study_script(..., 'SPMFiles', spm_files, 'SPMRuns', spm_runs)
%
% See also: hrf_validate_spm_inputs, hrf_fit_wholebrain_stats, hrf_write_slurm_study_script.

p = inputParser;
p.addRequired('fmri_files', @(x) ischar(x) || iscell(x) || isstring(x));
p.addRequired('spm_files',  @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(fmri_files, spm_files, varargin{:});
verbose = logical(p.Results.Verbose);
if ~isempty(p.Results.doverbose), verbose = logical(p.Results.doverbose); end

fmri_files = local_to_cellstr(fmri_files);
spm_files  = local_to_cellstr(spm_files);
n = numel(fmri_files);
if numel(spm_files) ~= n
    error('hrf_resolve_spm_runs:LengthMismatch', ...
        'fmri_files (%d) and spm_files (%d) must match.', n, numel(spm_files));
end

% Cache loaded SPMs (sessions repeat across runs).
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

spm_runs = nan(1, n);
n_tp = nan(n, 1); n_runs = nan(n, 1);
method = strings(n, 1); message = strings(n, 1);

for i = 1:n
    [spm_runs(i), n_tp(i), n_runs(i), method(i), message(i)] = ...
        local_resolve_one(fmri_files{i}, spm_files{i}, cache);
end

report = table((1:n)', fmri_files(:), spm_files(:), n_tp, n_runs, spm_runs(:), method, message, ...
    'VariableNames', {'index', 'fmri_file', 'spm_file', 'n_tp', 'n_runs', 'spm_run', 'method', 'message'});

if verbose
    disp(report(:, {'index', 'n_tp', 'n_runs', 'spm_run', 'method'}));
    n_bad = sum(isnan(spm_runs));
    fprintf('hrf_resolve_spm_runs: %d/%d resolved, %d unresolved.\n', n - n_bad, n, n_bad);
    for i = find(isnan(spm_runs))
        fprintf('  [%d] %s: %s\n', i, local_short(fmri_files{i}), message(i));
    end
end
end


% =========================================================================
function [run, n_tp, n_runs, method, msg] = local_resolve_one(fmri_file, spm_file, cache)
run = NaN; n_tp = NaN; n_runs = NaN; method = ""; msg = "";

if isempty(strtrim(char(spm_file)))
    run = 1; method = "no-spm"; msg = "no SPM -> run 1 (Tier A fallback)"; return
end

try
    SPM = local_load_spm(spm_file, cache);
catch err
    method = "spm-error"; msg = string(err.message); return
end

[run_sources, n_runs] = local_run_source_names(SPM);
if n_runs == 0
    method = "no-sessions"; msg = "SPM has no Sess/nscan structure"; return
end

fbase = local_basename(fmri_file);
ftok = local_run_token(fbase);

% 1. run-token match
if ~isempty(ftok)
    hit = [];
    for r = 1:n_runs
        if any(strcmp(ftok, run_sources{r}.run_tokens))
            hit(end + 1) = r; %#ok<AGROW>
        end
    end
    if isscalar(hit)
        run = hit; method = "run-token"; msg = sprintf("matched run-%s", ftok); return
    elseif numel(hit) > 1
        method = "ambiguous"; msg = sprintf("run-%s matches %d SPM runs", ftok, numel(hit)); return
    end
end

% 2. basename containment
hit = [];
for r = 1:n_runs
    if any(contains(run_sources{r}.bases, fbase)) || any(strcmp(run_sources{r}.bases, fbase))
        hit(end + 1) = r; %#ok<AGROW>
    end
end
if isscalar(hit)
    run = hit; method = "basename"; msg = "matched by source filename"; return
end

% 3. unique scan-count match
try
    n_tp = local_ntp(fmri_file);
catch
    n_tp = NaN;
end
if ~isnan(n_tp)
    counts = cellfun(@(s) s.nscan, run_sources);
    hit = find(counts == n_tp);
    if isscalar(hit)
        run = hit; method = "scan-count"; msg = sprintf("only run %d has %d scans", hit, n_tp); return
    elseif numel(hit) > 1
        method = "ambiguous"; msg = sprintf("%d runs have %d scans; set SPMRuns by hand", numel(hit), n_tp); return
    end
end

method = "nomatch"; msg = "no run-token, filename, or unique scan-count match";
end

function SPM = local_load_spm(spm_file, cache)
key = char(spm_file);
if isKey(cache, key), SPM = cache(key); return; end
if exist(key, 'file') ~= 2, error('SPM.mat not found: %s', key); end
S = load(key, 'SPM');
if ~isfield(S, 'SPM'), error('no SPM variable in %s', key); end
SPM = S.SPM;
cache(key) = SPM; %#ok<NASGU>
end

function [run_sources, n_runs] = local_run_source_names(SPM)
% Per run: the source image basenames, their run tokens, and scan count.
run_sources = {};
n_runs = 0;
if ~isfield(SPM, 'Sess') || isempty(SPM.Sess), return; end
n_runs = numel(SPM.Sess);

% Per-scan source filenames, if available.
vy_names = {};
if isfield(SPM, 'xY') && isfield(SPM.xY, 'VY') && ~isempty(SPM.xY.VY)
    try
        vy_names = {SPM.xY.VY.fname};
    catch
        vy_names = {};
    end
elseif isfield(SPM, 'xY') && isfield(SPM.xY, 'P') && ~isempty(SPM.xY.P)
    vy_names = cellstr(SPM.xY.P);
end

run_sources = cell(1, n_runs);
for r = 1:n_runs
    rows = SPM.Sess(r).row;
    s = struct('bases', {{}}, 'run_tokens', {{}}, 'nscan', numel(rows));
    if ~isempty(vy_names) && max(rows) <= numel(vy_names)
        names = unique(cellfun(@local_basename, vy_names(rows), 'UniformOutput', false), 'stable');
        s.bases = names;
        toks = cellfun(@local_run_token, names, 'UniformOutput', false);
        s.run_tokens = unique(toks(~cellfun(@isempty, toks)), 'stable');
    end
    run_sources{r} = s;
end
end

function n_tp = local_ntp(fmri_file)
info = niftiinfo(char(fmri_file));
sz = info.ImageSize;
if numel(sz) < 4, error('not 4D'); end
n_tp = sz(4);
end

function tok = local_run_token(name)
t = regexp(char(name), 'run[-_]?([0-9A-Za-z]+)', 'tokens', 'once');
if isempty(t), tok = ''; else, tok = t{1}; end
end

function b = local_basename(p)
[~, nm, ext] = fileparts(char(p));
b = [nm ext];
b = regexprep(b, ',\d+$', '');          % strip SPM volume index ',N'
b = regexprep(b, '\.nii(\.gz)?$', '');  % strip nifti extension
end

function c = local_to_cellstr(x)
if ischar(x), c = {x};
elseif isstring(x), c = cellstr(x);
elseif iscell(x), c = cellfun(@(v) char(string(v)), x, 'UniformOutput', false);
else, error('hrf_resolve_spm_runs:BadType', 'Expected char, string, or cell of paths.');
end
c = c(:)';
end

function s = local_short(p)
[~, nm, ext] = fileparts(char(p));
s = [nm ext];
end
