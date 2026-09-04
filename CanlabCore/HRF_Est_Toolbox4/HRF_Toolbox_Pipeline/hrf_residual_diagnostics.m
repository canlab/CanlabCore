function T = hrf_residual_diagnostics(source, varargin)
%HRF_RESIDUAL_DIAGNOSTICS Residual-based misspecification metrics from HRF fits.
%
% The residual-based half of the Phase 4 misspecification pipeline. Unlike
% hrf_misspec_metrics (which compares the estimated HRF *curve* to a
% reference and works from the whole-brain score CSVs), these metrics need
% the residual TIME COURSE of each fit, which only exists in the 1D
% extracted-signal fits stored in a results / study struct -- not in the
% whole-brain NIfTIs or their score CSVs.
%
% For each (subject, signal, model) fit it reports:
%   * mis_modeling_p   - ResidScan p-value: probability the smoothed
%                        residual's max is that large by chance. SMALL =
%                        task-locked structure remains in the residual =>
%                        the model is misspecified. (Lindquist & Loh 2007.)
%   * power_loss       - PowerLoss: efficiency lost relative to a flexible
%                        sFIR baseline (Loh/Lindquist/Wager). Larger = worse.
%   * resid_task_corr  - |corr(residual, combined stick function)|: direct
%                        check for residual variance tracking task timing.
%   * resid_acf1       - lag-1 autocorrelation of the residual.
%   * durbin_watson    - Durbin-Watson statistic (~2 = white; <2 = positive
%                        autocorrelation, i.e. structured residual).
%   * mse, dfe, n      - fit error variance, error df, n time points.
%
% mis_modeling_p and power_loss are harvested from the stored fit fields
% (computed at fit time with ResidScan FWHM = 4) unless 'Recompute' is
% true, in which case ResidScan is re-run with 'ResidScanFWHM'.
%
% Usage
% -----
%   T = hrf_residual_diagnostics(results)                 % one subject
%   T = hrf_residual_diagnostics(study)                   % study.results{:}
%   T = hrf_residual_diagnostics(study, 'Signal', 'signatures')
%   T = hrf_residual_diagnostics(study, 'Models', {'sfir','canonical'})
%   T = hrf_residual_diagnostics(study, 'Recompute', true, 'ResidScanFWHM', 6)
%
% Input dispatch (first arg)
% --------------------------
%   results struct (has .fits)            -> one subject
%   study struct  (has .results cell)     -> iterate, subject from
%                                            .subject_ids / results{i}.subject_id
%   cell array of results structs         -> iterate
%   struct array .label + .study/.results -> multi-study, rows tagged
%                                            study_label
%
% Name-value parameters
% ---------------------
%   'Signal'        - which 1D signals to diagnose:
%                       'mean' (default) -> results.fits (the mean/selected
%                                           signal)
%                       'signatures'     -> every results.fits_by_signature.*
%                       cellstr          -> named signature fields
%   'Models'        - cellstr filter, e.g. {'sfir','canonical'}. Default all
%                     models present in the fit struct.
%   'Recompute'     - false (default): harvest stored .mis_modeling_p /
%                     .power_loss. true: re-run ResidScan (FWHM below) and,
%                     when the inputs are available, PowerLoss.
%   'ResidScanFWHM' - FWHM (time units) for a recomputed ResidScan. Default 4.
%
% Output
% ------
% Long table; one row per (subject, signal, model) fit:
%   subject, study_label (when multi-study), signal, model,
%   n, dfe, mse, resid_task_corr, resid_acf1, durbin_watson,
%   mis_modeling_p, power_loss
% Consistent with the other Phase 4 / curve tables, so
% hrf_curve_summary_groupstats and friends compose over it.
%
% See also: hrf_misspec_metrics, ResidScan, PowerLoss, run_hrf_pipeline.

p = inputParser;
p.addRequired('input_source');
p.addParameter('Signal', 'mean', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('Models', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Recompute', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('ResidScanFWHM', 4, @(x) isscalar(x) && x > 0);
p.parse(source, varargin{:});
opts = p.Results;

T = local_empty_resid_table();

% Multi-study struct array: a non-scalar struct with a .label field, where
% each element carries either a .study sub-struct or is itself study-like.
if isstruct(source) && ~isscalar(source) && isfield(source, 'label')
    chunks = cell(numel(source), 1);
    for i = 1:numel(source)
        if isfield(source, 'study')
            sub = source(i).study;
        else
            sub = source(i);
        end
        Ti = hrf_residual_diagnostics(sub, varargin{:});
        if ~isempty(Ti) && height(Ti) > 0
            Ti.study_label = repmat(string(source(i).label), height(Ti), 1);
        end
        chunks{i} = Ti;
    end
    chunks = chunks(~cellfun(@isempty, chunks));
    if ~isempty(chunks), T = vertcat(chunks{:}); end
    return
end

% Normalize to a list of (results, subject_id) pairs.
[results_list, subject_ids] = local_results_list(source);

model_filter = local_to_cell(opts.Models);
rows = {};
for i = 1:numel(results_list)
    R = results_list{i};
    subj = subject_ids{i};
    signal_fits = local_collect_signal_fits(R, opts.Signal);
    for sf = 1:numel(signal_fits)
        signal_name = signal_fits(sf).name;
        fits = signal_fits(sf).fits;
        sticks = local_combined_stick(R);
        model_names = fieldnames(fits);
        for mi = 1:numel(model_names)
            mname = model_names{mi};
            if ~isempty(model_filter) && ~any(strcmpi(mname, model_filter)), continue; end
            fit = fits.(mname);
            if ~isstruct(fit) || ~isfield(fit, 'residual') || isempty(fit.residual), continue; end
            m = local_resid_metrics(fit, sticks, R, opts);
            rows{end + 1} = local_resid_row(subj, signal_name, mname, m); %#ok<AGROW>
        end
    end
end
if ~isempty(rows)
    T = vertcat(rows{:});
end
end


% =========================================================================
% Metric computation
% =========================================================================
function m = local_resid_metrics(fit, sticks, R, opts)
m = struct('n', NaN, 'dfe', NaN, 'mse', NaN, ...
    'resid_task_corr', NaN, 'resid_acf1', NaN, 'durbin_watson', NaN, ...
    'mis_modeling_p', NaN, 'power_loss', NaN);

e = double(fit.residual(:));
m.n = numel(e);
if isfield(fit, 'dfe') && ~isempty(fit.dfe), m.dfe = double(fit.dfe); end
if isfield(fit, 'mse') && ~isempty(fit.mse), m.mse = double(fit.mse); end

% Residual autocorrelation (lag-1) and Durbin-Watson.
if numel(e) >= 3 && std(e) > 0
    e0 = e(1:end-1); e1 = e(2:end);
    if std(e0) > 0 && std(e1) > 0
        m.resid_acf1 = corr(e0, e1);
    end
    m.durbin_watson = sum(diff(e) .^ 2) / sum(e .^ 2);
end

% Residual-task correlation (does residual variance track task timing?).
if ~isempty(sticks) && numel(sticks) == numel(e) && std(sticks) > 0 && std(e) > 0
    m.resid_task_corr = abs(corr(e, sticks(:)));
end

% ResidScan / PowerLoss: harvest stored, or recompute.
if logical(opts.Recompute)
    try
        m.mis_modeling_p = ResidScan(e, opts.ResidScanFWHM);
    catch
        m.mis_modeling_p = NaN;
    end
    m.power_loss = local_recompute_powerloss(fit, R);
else
    if isfield(fit, 'mis_modeling_p'), m.mis_modeling_p = double(fit.mis_modeling_p); end
    if isfield(fit, 'power_loss'),     m.power_loss     = double(fit.power_loss); end
end
end


function pl = local_recompute_powerloss(fit, R)
% Re-run PowerLoss when the necessary inputs are reachable from the
% results struct. Falls back to the stored value (or NaN) otherwise.
pl = NaN;
if isfield(fit, 'power_loss') && ~isempty(fit.power_loss), pl = double(fit.power_loss); end
have = isfield(fit, 'residual') && isfield(fit, 'fit') && ...
    isfield(R, 'timeseries') && isfield(R, 'stick_functions') && ...
    isfield(R, 'settings') && isfield(R.settings, 'TR');
if ~have, return; end
try
    e = double(fit.residual(:));
    f = double(fit.fit(:));
    moddf = numel(e) - size(fit.hrf, 1);
    pl = PowerLoss(e, f, max(moddf, 1), double(R.timeseries(:)), ...
        R.settings.TR, R.stick_functions, 0.001);
catch
    % keep the stored value already in pl
end
end


% =========================================================================
% Input normalization
% =========================================================================
function [results_list, subject_ids] = local_results_list(source)
results_list = {};
subject_ids = {};
if iscell(source)
    for i = 1:numel(source)
        if isempty(source{i}), continue; end
        results_list{end + 1} = source{i}; %#ok<AGROW>
        subject_ids{end + 1} = local_subject_of(source{i}, i); %#ok<AGROW>
    end
elseif isstruct(source) && isfield(source, 'results') && iscell(source.results)
    % A study struct.
    ids = {};
    if isfield(source, 'subject_ids'), ids = source.subject_ids; end
    for i = 1:numel(source.results)
        if isempty(source.results{i}), continue; end
        results_list{end + 1} = source.results{i}; %#ok<AGROW>
        if numel(ids) >= i && ~isempty(ids{i})
            subject_ids{end + 1} = char(string(ids{i})); %#ok<AGROW>
        else
            subject_ids{end + 1} = local_subject_of(source.results{i}, i); %#ok<AGROW>
        end
    end
elseif isstruct(source) && isfield(source, 'fits')
    % A single results struct.
    results_list = {source};
    subject_ids = {local_subject_of(source, 1)};
else
    error('hrf_residual_diagnostics:UnknownSource', ...
        ['First arg must be a results struct (with .fits), a study struct ' ...
         '(with .results), or a cell array of results structs.']);
end
end


function s = local_subject_of(R, idx)
s = sprintf('sub-%03d', idx);
if isstruct(R)
    if isfield(R, 'subject_id') && ~isempty(R.subject_id)
        s = char(string(R.subject_id));
    elseif isfield(R, 'subject') && ~isempty(R.subject)
        s = char(string(R.subject));
    end
end
end


function signal_fits = local_collect_signal_fits(R, signal_opt)
% Returns a struct array of (name, fits) to diagnose.
signal_fits = struct('name', {}, 'fits', {});
if ischar(signal_opt) || (isstring(signal_opt) && isscalar(signal_opt))
    so = lower(char(signal_opt));
    if strcmp(so, 'mean')
        if isfield(R, 'fits') && isstruct(R.fits)
            signal_fits(end + 1) = struct('name', 'mean', 'fits', R.fits);
        end
        return
    elseif strcmp(so, 'signatures')
        signal_fits = local_all_signature_fits(R);
        return
    end
    % a single named signature
    signal_fits = local_named_signature_fits(R, {char(signal_opt)});
    return
end
% cellstr / string array of names
signal_fits = local_named_signature_fits(R, cellstr(string(signal_opt)));
end


function sf = local_all_signature_fits(R)
sf = struct('name', {}, 'fits', {});
if isfield(R, 'fits_by_signature') && isstruct(R.fits_by_signature)
    fns = fieldnames(R.fits_by_signature);
    for i = 1:numel(fns)
        sf(end + 1) = struct('name', fns{i}, 'fits', R.fits_by_signature.(fns{i})); %#ok<AGROW>
    end
end
end


function sf = local_named_signature_fits(R, names)
sf = struct('name', {}, 'fits', {});
if ~isfield(R, 'fits_by_signature') || ~isstruct(R.fits_by_signature), return; end
for i = 1:numel(names)
    f = matlab.lang.makeValidName(names{i});
    if isfield(R.fits_by_signature, f)
        sf(end + 1) = struct('name', names{i}, 'fits', R.fits_by_signature.(f)); %#ok<AGROW>
    end
end
end


function sticks = local_combined_stick(R)
% Sum the per-condition stick functions into one task-timing regressor.
sticks = [];
if isfield(R, 'stick_functions') && iscell(R.stick_functions) && ~isempty(R.stick_functions)
    Runc = R.stick_functions;
    len = numel(Runc{1});
    sticks = zeros(len, 1);
    for c = 1:numel(Runc)
        v = double(Runc{c}(:));
        if numel(v) == len, sticks = sticks + v; end
    end
end
end


% =========================================================================
% Table construction / utilities
% =========================================================================
function T = local_empty_resid_table()
T = table('Size', [0 12], ...
    'VariableTypes', {'string','string','string','string', ...
        'double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'subject','study_label','signal','model', ...
        'n','dfe','mse','resid_task_corr','resid_acf1','durbin_watson', ...
        'mis_modeling_p','power_loss'});
end


function row = local_resid_row(subj, signal_name, mname, m)
row = table(string(subj), string(""), string(signal_name), string(mname), ...
    m.n, m.dfe, m.mse, m.resid_task_corr, m.resid_acf1, m.durbin_watson, ...
    m.mis_modeling_p, m.power_loss, ...
    'VariableNames', {'subject','study_label','signal','model', ...
        'n','dfe','mse','resid_task_corr','resid_acf1','durbin_watson', ...
        'mis_modeling_p','power_loss'});
end


function c = local_to_cell(x)
if isempty(x), c = {};
elseif ischar(x), c = {x};
elseif isstring(x), c = cellstr(x);
else, c = x;
end
end
