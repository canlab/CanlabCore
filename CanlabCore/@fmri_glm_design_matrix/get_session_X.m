function [Xs, delta, C, B, names] = get_session_X(obj, s)
% Get design matrix (predictors) for one session of an fmri_glm_design_matrix
% object, using the basis functions defined in the object and the onsets for
% one session (s).
%
% :Usage:
% ::
%
%     [Xs, delta, C, B, names] = get_session_X(obj, session number)
%
% Method: a high-resolution (microtime) pipeline. Onsets and event durations
% are sampled into a high-resolution indicator at RES samples per second
% (RES = 16, TR-independent), convolved with the per-condition basis set
% (resampled to the same resolution), and downsampled to the TR with
% getPredictors. This matches onsets2fmridesign and is robust to fractional
% TRs. Each condition may have its own basis set (one or more basis
% functions), and parametric modulators are supported.
%
% :Outputs:
%   **Xs:**    [nscan x k] of-interest predictors (conditions x basis fns, then PM regressors)
%   **delta:** [nscan x nconds] TR-resolution onset indicator
%   **C:**     user-specified covariates for this session ([] if none)
%   **B:**     [nscan x 1] baseline (intercept) for this session
%   **names:** 1 x k cell of regressor names
%
% :See also: getPredictors, onsets2fmridesign, fmri_glm_design_matrix.build

% ..
%    High-resolution microtime sampling rate (samples per second). TR-
%    independent, matching onsets2fmridesign (res = 16).
% ..
res = 16;

nsess = length(obj.Sess);
if s > nsess, error('Session %3.0f does not exist', s); end

TR = obj.xY.RT;
nscan = obj.nscan(s);

if isempty(nscan) || nscan < 1
    error('fmri_glm_design_matrix:get_session_X', ...
        'obj.nscan(%d) is empty or invalid; set the number of scans for each session.', s);
end

nconds = length(obj.Sess(s).U);

% Units of onsets/durations: 'secs' (default) or 'scans'/'TR'
units = '';
if ~isempty(obj.xBF) && isfield(obj.xBF(1), 'UNITS') && ~isempty(obj.xBF(1).UNITS)
    units = obj.xBF(1).UNITS;
end
if any(strcmpi(units, {'scans', 'tr', 'trs'})), to_sec = TR; else, to_sec = 1; end

% Session length in high-resolution samples. nscan*TR*res is an integer
% multiple of res*TR, so getPredictors downsamples cleanly to nscan rows.
len_hires = round(nscan * TR * res);

% ----------------------------------------------
% Of-interest predictors: one block per condition (HRF-convolved)
% ----------------------------------------------
Xs       = cell(1, nconds);
delta    = zeros(nscan, nconds);
condname = cell(1, nconds);

for i = 1:nconds

    U = obj.Sess(s).U(i);

    condname{i} = local_getfield(U, 'name', sprintf('R%d', i));
    if iscell(condname{i})
        if isempty(condname{i}), condname{i} = sprintf('R%d', i); else, condname{i} = condname{i}{1}; end
    end

    ons_sec = U.ons(:) * to_sec;

    % Durations (sec): empty/scalar -> applied to all events; else per-event
    dur = local_getfield(U, 'dur', []);
    if isempty(dur), dur = 0; end
    dur = dur(:);
    if isscalar(dur), dur = repmat(dur, numel(ons_sec), 1); end
    dur_sec = dur * to_sec;

    % High-resolution onset/epoch indicator for this condition
    dhr = local_hires_indicator(ons_sec, dur_sec, res, len_hires);

    % Basis set for this condition, resampled to res samples/second
    bf = local_bf_at_res(obj.xBF(min(i, numel(obj.xBF))), res);

    Xi = getPredictors(dhr, bf, 'dsrate', res * TR, 'force_delta');
    Xs{i} = local_fit_rows(Xi, nscan);

    % TR-resolution indicator (for reference/output)
    tr_idx = round(ons_sec / TR) + 1;
    tr_idx = tr_idx(tr_idx >= 1 & tr_idx <= nscan);
    delta(tr_idx, i) = 1;
end

Xs = cat(2, Xs{:});

% ----------------------------------------------
% Parametric modulators (appended after the of-interest block)
% ----------------------------------------------
Xs_pm = cell(1, nconds);

for i = 1:nconds

    U = obj.Sess(s).U(i);
    if ~local_has_pm(U), continue, end

    ons_sec = U.ons(:) * to_sec;
    dur = local_getfield(U, 'dur', []); if isempty(dur), dur = 0; end
    dur = dur(:); if isscalar(dur), dur = repmat(dur, numel(ons_sec), 1); end
    dhr = local_hires_indicator(ons_sec, dur * to_sec, res, len_hires);

    bf = local_bf_at_res(obj.xBF(min(i, numel(obj.xBF))), res);

    pm_vals = U.P.P(:);
    Xpm = getPredictors(dhr, bf, 'dsrate', res * TR, 'force_delta', ...
        'parametric_singleregressor', {pm_vals});
    Xs_pm{i} = local_fit_rows(Xpm, nscan);
end

Xs = [Xs cat(2, Xs_pm{:})];

% ----------------------------------------------
% Covariates and baseline
% ----------------------------------------------
C = [];
if ~isempty(obj.Sess(s).C) && isfield(obj.Sess(s).C, 'C')
    C = obj.Sess(s).C.C;
end

B = ones(nscan, 1);

% ----------------------------------------------
% Names: one per basis function per condition, then PM names
% ----------------------------------------------
names = {};
for i = 1:nconds
    nbf = size(local_bf_at_res(obj.xBF(min(i, numel(obj.xBF))), res), 2);
    for j = 1:nbf
        names{end + 1} = replaceblanks([condname{i} ' BF' num2str(j)]); %#ok<AGROW>
    end
end

for i = 1:nconds
    U = obj.Sess(s).U(i);
    if ~local_has_pm(U), continue, end
    nbf = size(local_bf_at_res(obj.xBF(min(i, numel(obj.xBF))), res), 2);
    for j = 1:nbf
        names{end + 1} = replaceblanks([U.P.name ' BF' num2str(j)]); %#ok<AGROW>
    end
end

names = names(:)';

end % get_session_X


% =========================================================================
% Local helpers
% =========================================================================
function dhr = local_hires_indicator(ons_sec, dur_sec, res, len_hires)
% Build a high-resolution onset/epoch indicator (len_hires x 1). Events with
% positive duration become epochs; zero-duration events are single impulses.
dhr = zeros(len_hires, 1);
for j = 1:numel(ons_sec)
    k0 = round(ons_sec(j) * res) + 1;          % first scan is time 0 -> element 1
    if k0 < 1 || k0 > len_hires, continue, end
    if numel(dur_sec) >= j && dur_sec(j) > 0
        k1 = min(len_hires, round(k0 + dur_sec(j) * res));
        dhr(k0:k1) = 1;
    else
        dhr(k0) = 1;
    end
end
end


function bf = local_bf_at_res(xBF, res)
% Return the basis-set matrix resampled to res samples/second, using the
% stored sampling interval xBF.dt (seconds/sample). No-op when already at res.
bf = xBF.bf;
dt0 = [];
if isfield(xBF, 'dt') && ~isempty(xBF.dt), dt0 = xBF.dt; end
if isempty(dt0) || dt0 <= 0, dt0 = 1 / res; end
if abs(dt0 - 1 / res) < 1e-9, return, end

n0 = size(bf, 1);
t0 = (0:n0 - 1)' * dt0;                          % times of stored samples (sec)
tnew = (0:1 / res:t0(end))';                     % resample at res Hz
bf = interp1(t0, bf, tnew, 'linear', 'extrap');
end


function X = local_fit_rows(X, nscan)
% Trim or zero-pad X to exactly nscan rows.
if size(X, 1) >= nscan
    X = X(1:nscan, :);
else
    X(end + 1:nscan, :) = 0;
end
end


function tf = local_has_pm(U)
tf = isfield(U, 'P') && ~isempty(U.P) && isfield(U.P, 'P') && ~isempty(U.P.P);
end


function v = local_getfield(s, f, default)
% Return s.(f) if the field exists and is non-empty, else default.
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = default;
end
end


function myname = replaceblanks(myname)
myname(myname == ' ') = '-';
tmp = diff(double(myname)) == 0;
myname([false tmp] & myname == '-') = [];
end
