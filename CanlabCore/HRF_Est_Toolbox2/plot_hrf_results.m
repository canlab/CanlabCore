function ax = plot_hrf_results(results, varargin)
%PLOT_HRF_RESULTS Quick plotting for run_hrf_pipeline output.
% ax = plot_hrf_results(results, 'Model', 'sfir', 'Conditions', [1 2 3])
% ax = plot_hrf_results(results, 'Model', 'sfir', 'Signature', 'NPS')

p = inputParser;
p.addRequired('results', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('Conditions', [], @(x) isempty(x) || isnumeric(x) || iscell(x));
p.addParameter('Signature', '', @(x) ischar(x) || isstring(x));
p.parse(results, varargin{:});
opts = p.Results;

model_name = char(opts.Model);

if isempty(opts.Conditions)
    cond_idx = 1:numel(results.conditions);
elseif isnumeric(opts.Conditions)
    cond_idx = opts.Conditions;
else
    cond_idx = find(ismember(results.conditions, opts.Conditions));
end

if ~isempty(opts.Signature) && isfield(results, 'fits_by_signature')
    sig = char(opts.Signature);
    if ~isfield(results.fits_by_signature, sig)
        error('Signature %s not found in results.fits_by_signature', sig);
    end
    fit_struct = results.fits_by_signature.(sig);
    ttl = sprintf('%s model HRF by condition (%s signature)', upper(model_name), sig);
else
    fit_struct = results.fits;
    ttl = sprintf('%s model HRF by condition', upper(model_name));
end

if ~isfield(fit_struct, model_name)
    error('Model %s not available in results.', model_name);
end

hrf = fit_struct.(model_name).hrf;
figure; ax = axes; hold(ax, 'on');
plot(ax, hrf(:, cond_idx), 'LineWidth', 1.5);
legend(ax, format_strings_for_legend(results.conditions(cond_idx)), 'Interpreter', 'none');
title(ax, ttl, 'Interpreter', 'none');
xlabel(ax, 'HRF time bins'); ylabel(ax, 'Response (a.u.)');
hline(0, 'k-');
end
