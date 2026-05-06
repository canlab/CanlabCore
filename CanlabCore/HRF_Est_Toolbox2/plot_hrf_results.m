function ax = plot_hrf_results(results, varargin)
%PLOT_HRF_RESULTS Quick plotting for run_hrf_pipeline output.
% ax = plot_hrf_results(results, 'Model', 'sfir', 'Conditions', [1 2 3])
% ax = plot_hrf_results(results, 'Model', 'sfir', 'Signature', 'NPS')

p = inputParser;
p.addRequired('results', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('Conditions', [], @(x) isempty(x) || isnumeric(x) || iscell(x) || isstring(x));
p.addParameter('Signature', '', @(x) ischar(x) || isstring(x));
p.parse(results, varargin{:});
opts = p.Results;

model_name = char(opts.Model);

if isempty(opts.Conditions)
    cond_idx = 1:numel(results.conditions);
elseif isnumeric(opts.Conditions)
    cond_idx = opts.Conditions;
else
    cond_names = cellstr(string(opts.Conditions));
    cond_idx = find(ismember(results.conditions, cond_names));
end

if ~isempty(opts.Signature) && isfield(results, 'fits_by_signature')
    sig = char(opts.Signature);
    sig_field = local_signature_field(results, sig);
    if isempty(sig_field)
        error('Signature %s not found in results.fits_by_signature', sig);
    end
    fit_struct = results.fits_by_signature.(sig_field);
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

% SE shading: use across-signature variability when available and Signature not explicitly set
se_mat = [];
if isempty(opts.Signature) && isfield(results, 'fits_by_signature') && numel(fieldnames(results.fits_by_signature)) > 1
    fn = fieldnames(results.fits_by_signature);
    stack = [];
    for i = 1:numel(fn)
        if isfield(results.fits_by_signature.(fn{i}), model_name)
            stack(:, :, i) = results.fits_by_signature.(fn{i}).(model_name).hrf; %#ok<AGROW>
        end
    end
    if ~isempty(stack)
        se_mat = std(stack, 0, 3) ./ sqrt(size(stack, 3));
    end
end

colors = lines(numel(cond_idx));
for k = 1:numel(cond_idx)
    c = cond_idx(k);
    y = hrf(:, c);
    if ~isempty(se_mat)
        se = se_mat(:, c);
        x = (1:numel(y))';
        fill([x; flipud(x)], [y+se; flipud(y-se)], colors(k,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    plot(ax, y, 'LineWidth', 1.8, 'Color', colors(k,:));
end
legend(ax, format_strings_for_legend(results.conditions(cond_idx)), 'Interpreter', 'none');
title(ax, ttl, 'Interpreter', 'none');
xlabel(ax, 'HRF time bins'); ylabel(ax, 'Response (a.u.)');
hline(0, 'k-');
end

function sig_field = local_signature_field(results, sig)
sig_field = '';
if isfield(results.fits_by_signature, sig)
    sig_field = sig;
    return
end

candidate = matlab.lang.makeValidName(sig);
if isfield(results.fits_by_signature, candidate)
    sig_field = candidate;
    return
end

if isfield(results, 'signature_meta') && isfield(results.signature_meta, 'selected_signatures') && ...
        isfield(results.signature_meta, 'selected_signature_fields')
    names = cellstr(string(results.signature_meta.selected_signatures));
    fields = cellstr(string(results.signature_meta.selected_signature_fields));
    idx = find(strcmp(names, sig), 1);
    if ~isempty(idx) && idx <= numel(fields) && isfield(results.fits_by_signature, fields{idx})
        sig_field = fields{idx};
    end
end
end
