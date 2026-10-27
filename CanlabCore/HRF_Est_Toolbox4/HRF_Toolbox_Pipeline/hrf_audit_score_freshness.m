function report = hrf_audit_score_freshness(output_dir, varargin)
% Post-run completeness/freshness audit of HRF score CSVs against their maps.
%
% Catches the failure modes a "the job finished" check misses: a score CSV
% that is STALE relative to the map it was scored from (older than the
% beta/t NIfTI -> scored from a since-overwritten fit), score CSVs with an
% inconsistent atlas-column count across subjects (a partial / wrong-config
% scoring), and tasks whose logs report a scoring failure. Run this after
% every SLURM batch before trusting the group plots.
%
% :Usage:
% ::
%     report = hrf_audit_score_freshness(output_dir)
%     report = hrf_audit_score_freshness(output_dir, 'ScoreObjects', {'beta','t'})
%
% :Inputs:
%
%   **output_dir:**
%        folder holding the *_<object>_map_scores.csv and *_<object>.nii files
%        (the SLURM OutputDir). Multi-model files are *_<model>_<object>_*.
%
% :Optional Inputs:
%
%   **'ScoreObjects':**  cellstr, default {'beta','t'}.
%   **'LogDir':**        SLURM log dir, default <output_dir>/logs. Scanned for
%                        'Score failure' / 'Cannot score' lines.
%   **'ExpectedAtlasCols':** scalar; if given, CSVs whose atlas-column count
%                        differs are flagged. If empty, the modal count across
%                        CSVs is used as the reference.
%   **'Verbose' / 'doverbose':** print the summary (default true).
%
% :Outputs:
%
%   **report:**
%        table, one row per score CSV: object, csv_file, status
%        ('fresh' | 'STALE' | 'missing_nii' | 'unreadable'), csv_time,
%        nii_time, atlas_cols, atlas_col_flag, and (if logs found)
%        a log_scoring_failure flag.
%
% :Examples:
% ::
%     report = hrf_audit_score_freshness('/path/to/hrf_outputs');
%     bad = report(report.status~="fresh" | report.atlas_col_flag | report.log_scoring_failure, :);
%     assert(isempty(bad), 'Re-score the %d flagged outputs before plotting.', height(bad));
%
% See also: hrf_audit_slurm_outputs, hrf_score_one_prefix, hrf_apply_maps_to_wholebrain.

p = inputParser;
p.addRequired('output_dir', @(x) ischar(x) || isstring(x));
p.addParameter('ScoreObjects', {'beta', 't'}, @(x) ischar(x) || iscell(x) || isstring(x));
p.addParameter('LogDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('ExpectedAtlasCols', [], @(x) isempty(x) || isscalar(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(output_dir, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);
if ~isempty(opts.doverbose), verbose = logical(opts.doverbose); end

od = char(output_dir);
objs = cellstr(string(opts.ScoreObjects));
log_dir = char(opts.LogDir);
if isempty(log_dir), log_dir = fullfile(od, 'logs'); end

% --- gather scoring-failure task hints from logs (once) ---
log_failures = local_scan_logs(log_dir);

obj_col = {}; csv_col = {}; status = strings(0, 1);
csv_time = strings(0, 1); nii_time = strings(0, 1); atlas_cols = []; log_fail = false(0, 1);

for o = 1:numel(objs)
    obj = objs{o};
    csvs = dir(fullfile(od, sprintf('*_%s_map_scores.csv', obj)));
    for k = 1:numel(csvs)
        csv_name = csvs(k).name;
        csv_path = fullfile(csvs(k).folder, csv_name);
        % matching map: ..._<obj>_map_scores.csv -> ..._<obj>.nii
        nii_name = regexprep(csv_name, sprintf('_%s_map_scores\\.csv$', obj), sprintf('_%s.nii', obj));
        nii_path = fullfile(od, nii_name);
        nii_d = dir(nii_path);

        ncol = local_atlas_col_count(csv_path);
        if isnan(ncol)
            st = "unreadable"; ntime = "";
        elseif isempty(nii_d)
            st = "missing_nii"; ntime = "";
        else
            ntime = local_tstr(nii_d.datenum);
            if csvs(k).datenum < nii_d.datenum - 1/1440   % CSV >1 min older than its map
                st = "STALE";
            else
                st = "fresh";
            end
        end

        obj_col{end + 1, 1} = obj; %#ok<AGROW>
        csv_col{end + 1, 1} = csv_name; %#ok<AGROW>
        status(end + 1, 1) = st; %#ok<AGROW>
        csv_time(end + 1, 1) = local_tstr(csvs(k).datenum); %#ok<AGROW>
        nii_time(end + 1, 1) = ntime; %#ok<AGROW>
        atlas_cols(end + 1, 1) = ncol; %#ok<AGROW>
        log_fail(end + 1, 1) = local_csv_failed_in_logs(csv_name, log_failures); %#ok<AGROW>
    end
end

% --- atlas-column consistency ---
ref = opts.ExpectedAtlasCols;
if isempty(ref)
    finite_counts = atlas_cols(~isnan(atlas_cols));
    if isempty(finite_counts), ref = NaN; else, ref = mode(finite_counts); end
end
atlas_col_flag = ~isnan(atlas_cols) & atlas_cols ~= ref;

report = table(obj_col, csv_col, status, csv_time, nii_time, atlas_cols, atlas_col_flag, log_fail, ...
    'VariableNames', {'object', 'csv_file', 'status', 'csv_time', 'nii_time', ...
    'atlas_cols', 'atlas_col_flag', 'log_scoring_failure'});

if verbose
    n = height(report);
    n_stale = sum(report.status == "STALE");
    n_missing = sum(report.status == "missing_nii" | report.status == "unreadable");
    n_atlas = sum(report.atlas_col_flag);
    n_log = sum(report.log_scoring_failure);
    fprintf('hrf_audit_score_freshness: %d score CSVs | reference atlas cols = %g\n', n, ref);
    fprintf('  STALE (older than their map): %d\n', n_stale);
    fprintf('  missing/unreadable map or CSV: %d\n', n_missing);
    fprintf('  atlas-col count != %g: %d\n', ref, n_atlas);
    fprintf('  subject flagged with a scoring failure in logs: %d\n', n_log);
    bad = report(report.status ~= "fresh" | report.atlas_col_flag | report.log_scoring_failure, :);
    if isempty(bad)
        fprintf('  ALL CLEAR -- every score CSV is fresh, consistent, and log-clean.\n');
    else
        fprintf('  ---- %d flagged outputs ----\n', height(bad));
        for b = 1:height(bad)
            fprintf('   [%s] %s  (atlas_cols=%g%s%s)\n', bad.status(b), bad.csv_file{b}, bad.atlas_cols(b), ...
                local_tag(bad.atlas_col_flag(b), ' ATLAS-COUNT'), local_tag(bad.log_scoring_failure(b), ' LOG-FAIL'));
        end
    end
end
end


% =========================================================================
function n = local_atlas_col_count(csv_path)
n = NaN;
fid = fopen(csv_path, 'r');
if fid < 0, return; end
c = onCleanup(@() fclose(fid));
hdr = fgetl(fid);
if ischar(hdr)
    cols = strsplit(hdr, ',');
    n = sum(startsWith(strtrim(cols), 'atlas_'));
end
end

function failmap = local_scan_logs(log_dir)
% Map run-token -> set of 'model/object' that the worker reported a scoring
% failure for, from the MOST RECENT SLURM job only (highest job id in the
% *_<jobid>_<task>.{out,err} names) so superseded runs don't taint the result.
% Precise (subject,model,object) matching lets us flag only the exact CSVs
% that failed, not every CSV of a subject that had any failure.
failmap = containers.Map('KeyType', 'char', 'ValueType', 'any');
if exist(log_dir, 'dir') ~= 7, return; end
files = [dir(fullfile(log_dir, '*.out')); dir(fullfile(log_dir, '*.err'))];
if isempty(files), return; end
jobids = nan(numel(files), 1);
for i = 1:numel(files)
    tok = regexp(files(i).name, '_(\d+)_\d+\.(out|err)$', 'tokens', 'once');
    if ~isempty(tok), jobids(i) = str2double(tok{1}); end
end
if any(~isnan(jobids))
    files = files(jobids == max(jobids));
end
for i = 1:numel(files)
    try
        txt = fileread(fullfile(files(i).folder, files(i).name));
    catch
        continue
    end
    banner = regexp(txt, 'Running HRF task \d+/\d+:\s*(\S+)\s*\|\s*(\S+)', 'tokens', 'once');
    if isempty(banner), continue; end
    runtoken = [banner{1} '_' banner{2}];
    pairs = regexp(txt, 'Score failure for model=(\w+), object=(\w+)', 'tokens');
    if isempty(pairs), continue; end
    if isKey(failmap, runtoken), s = failmap(runtoken); else, s = {}; end
    for j = 1:numel(pairs)
        s{end + 1} = [lower(pairs{j}{1}) '/' lower(pairs{j}{2})]; %#ok<AGROW>
    end
    failmap(runtoken) = unique(s);
end
end

function tf = local_csv_failed_in_logs(csv_name, failmap)
% Flag a CSV only if its exact (run-token, model, object) reported a failure.
tf = false;
if isempty(failmap) || ~isa(failmap, 'containers.Map'), return; end
t = regexp(csv_name, '^(.*)_hrf_(\w+)_(beta|t)_map_scores\.csv$', 'tokens', 'once');
if isempty(t)
    t = regexp(csv_name, '^(.*)_hrf_(beta|t)_map_scores\.csv$', 'tokens', 'once');
    if isempty(t), return; end
    runtoken = t{1}; key = lower(t{2});            % single-model: object only
else
    runtoken = t{1}; key = [lower(t{2}) '/' lower(t{3})];
end
if ~isKey(failmap, runtoken), return; end
s = failmap(runtoken);
tf = any(strcmp(key, s)) || any(endsWith(s, ['/' key]));
end

function s = local_tag(flag, label)
if flag, s = label; else, s = ''; end
end

function s = local_tstr(dn)
s = string(datetime(dn, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm'));
end
