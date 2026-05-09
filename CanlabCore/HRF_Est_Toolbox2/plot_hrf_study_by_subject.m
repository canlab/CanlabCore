function ax = plot_hrf_study_by_subject(study, varargin)
%PLOT_HRF_STUDY_BY_SUBJECT Plot HRFs by subject, averaging duplicate runs.

p = inputParser;
p.addRequired('study', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x));
p.addParameter('Signature', '', @(x) ischar(x) || isstring(x));
p.addParameter('Unit', 'subject', @(x) ischar(x) || isstring(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.parse(study, varargin{:});
opts = p.Results;
opts.Condition = char(opts.Condition);
opts.Signature = char(opts.Signature);
opts.Model = char(opts.Model);
opts.Unit = char(opts.Unit);

[Y_runs, run_subject_ids, condition_name, skipped] = local_collect_curves(study, opts);
if isempty(Y_runs)
    error('No valid HRF curves found. Skipped %d result(s).', numel(skipped));
end

[Y, subject_ids] = local_aggregate_unit(Y_runs, run_subject_ids, lower(char(opts.Unit)));

figure; ax = axes; hold(ax, 'on');
colors = lines(size(Y, 1));
for s = 1:size(Y, 1)
    plot(ax, Y(s, :)', 'Color', colors(s, :), 'DisplayName', subject_ids{s});
end
legend(ax, 'Interpreter', 'none');
title(ax, local_title(opts, condition_name), 'Interpreter', 'none');
xlabel(ax, 'HRF time bins');
ylabel(ax, 'Response (a.u.)');
hline(0, 'k-');
end

function [Y, subject_ids, condition_name, skipped] = local_collect_curves(study, opts)
Y = [];
subject_ids = {};
condition_name = char(opts.Condition);
skipped = struct('index', {}, 'subject', {}, 'reason', {});
model_name = opts.Model;
missing_policy = lower(char(opts.MissingPolicy));

for s = 1:numel(study.results)
    subject_id = local_subject_id(study, s);
    r = study.results{s};
    if isempty(r)
        skipped = local_skip(skipped, s, subject_id, 'empty result', missing_policy);
        continue
    end

    [fit_struct, ok, reason] = local_fit_struct(r, opts.Signature);
    if ~ok
        skipped = local_skip(skipped, s, subject_id, reason, missing_policy);
        continue
    end
    if ~isfield(fit_struct, model_name)
        skipped = local_skip(skipped, s, subject_id, sprintf('missing model %s', model_name), missing_policy);
        continue
    end

    if isempty(condition_name)
        c = 1;
        if isfield(r, 'conditions') && ~isempty(r.conditions), condition_name = r.conditions{1}; end
    else
        c = find(strcmp(r.conditions, condition_name), 1);
    end
    if isempty(c)
        skipped = local_skip(skipped, s, subject_id, sprintf('missing condition %s', condition_name), missing_policy);
        continue
    end

    y = fit_struct.(model_name).hrf(:, c)';
    if isempty(Y)
        Y = y;
    elseif size(Y, 2) == numel(y)
        Y(end + 1, :) = y; %#ok<AGROW>
    else
        skipped = local_skip(skipped, s, subject_id, 'HRF length mismatch', missing_policy);
        continue
    end
    subject_ids{end + 1, 1} = subject_id; %#ok<AGROW>
end
end

function [Y, subject_ids] = local_aggregate_unit(Y_runs, run_subject_ids, unit)
switch unit
    case 'run'
        Y = Y_runs;
        subject_ids = run_subject_ids(:);
    case 'subject'
        subject_ids = unique(run_subject_ids, 'stable');
        Y = nan(numel(subject_ids), size(Y_runs, 2));
        for i = 1:numel(subject_ids)
            Y(i, :) = local_mean_omitnan(Y_runs(strcmp(run_subject_ids, subject_ids{i}), :), 1);
        end
    otherwise
        error('Unknown Unit: %s. Use ''subject'' or ''run''.', unit);
end
end

function ttl = local_title(opts, condition_name)
sig = char(opts.Signature);
if isempty(sig)
    sig = 'selected signal';
end
ttl = sprintf('%s HRF by %s: %s (%s)', upper(char(opts.Model)), lower(char(opts.Unit)), condition_name, sig);
end

function [fit_struct, ok, reason] = local_fit_struct(r, signature_name)
fit_struct = struct();
ok = false;
reason = '';
if isempty(signature_name)
    if isfield(r, 'fits')
        fit_struct = r.fits;
        ok = true;
    else
        reason = 'missing fits';
    end
    return
end

if ~isfield(r, 'fits_by_signature')
    reason = 'missing fits_by_signature';
    return
end

sig_field = local_signature_field(r, signature_name);
if isempty(sig_field)
    reason = sprintf('missing signature %s', signature_name);
    return
end
fit_struct = r.fits_by_signature.(sig_field);
ok = true;
end

function sig_field = local_signature_field(r, sig)
sig_field = '';
if isfield(r.fits_by_signature, sig)
    sig_field = sig;
    return
end
candidate = matlab.lang.makeValidName(sig);
if isfield(r.fits_by_signature, candidate)
    sig_field = candidate;
    return
end
if isfield(r, 'signature_meta') && isfield(r.signature_meta, 'selected_signatures') && ...
        isfield(r.signature_meta, 'selected_signature_fields')
    names = cellstr(string(r.signature_meta.selected_signatures));
    fields = cellstr(string(r.signature_meta.selected_signature_fields));
    idx = find(strcmp(names, sig), 1);
    if ~isempty(idx) && idx <= numel(fields) && isfield(r.fits_by_signature, fields{idx})
        sig_field = fields{idx};
    end
end
end

function m = local_mean_omitnan(X, dim)
if nargin < 2, dim = 1; end
valid = ~isnan(X);
den = sum(valid, dim);
X(~valid) = 0;
m = sum(X, dim) ./ den;
m(den == 0) = NaN;
end

function skipped = local_skip(skipped, idx, subject_id, reason, missing_policy)
if strcmp(missing_policy, 'error')
    error('Skipping is disabled. Result %d (%s): %s', idx, subject_id, reason);
elseif ~strcmp(missing_policy, 'warn') && ~strcmp(missing_policy, 'silent')
    error('Unknown MissingPolicy: %s. Use ''warn'', ''silent'', or ''error''.', missing_policy);
end

skipped(end + 1) = struct('index', idx, 'subject', subject_id, 'reason', reason);
if strcmp(missing_policy, 'warn')
    warning('plot_hrf_study_by_subject:SkippingResult', 'Skipping result %d (%s): %s', idx, subject_id, reason);
end
end

function subject_id = local_subject_id(study, idx)
if isfield(study, 'subject_ids') && numel(study.subject_ids) >= idx
    subject_id = char(study.subject_ids{idx});
else
    subject_id = sprintf('sub-%03d', idx);
end
end
