function fits = hrf_fit_all_models(tc, TR, Runc, window_seconds, models, varargin)
%HRF_FIT_ALL_MODELS Fit supported HRF models with a consistent output API.

% Optional name/value:
%   'SuppressWarnings' (default true): temporarily suppress warning spam
%   'DependencyPolicy' (default 'skip'): 'skip' unavailable optional models
%       with a warning, or 'error' if a requested model dependency is missing.
p = inputParser;
p.addParameter('SuppressWarnings', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('DependencyPolicy', 'skip', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
suppress_warnings = logical(p.Results.SuppressWarnings);
dependency_policy = lower(char(p.Results.DependencyPolicy));
if ~ismember(dependency_policy, {'skip', 'error'})
    error('Unknown DependencyPolicy: %s. Use ''skip'' or ''error''.', dependency_policy);
end

models = lower(string(models));
models = local_filter_available_models(models, dependency_policy);
fits = struct();
len = length(tc);

if any(models == "logit")
    [h, fit, e, param] = run_fit(@() Fit_Logit2(tc, TR, Runc, window_seconds, 0), suppress_warnings);
    fits.logit = package_fit(h, fit, e, param, len, 7, tc, TR, Runc);
end
if any(models == "sfir")
    [h, fit, e, param] = run_fit(@() Fit_sFIR(tc, TR, Runc, window_seconds, 1), suppress_warnings);
    fits.sfir = package_fit(h, fit, e, param, len, window_seconds, tc, TR, Runc);
end
if any(models == "canonical")
    [h, fit, e, param, info] = run_fit(@() Fit_Canonical_HRF(tc, TR, Runc, window_seconds, 1), suppress_warnings);
    fits.canonical = package_fit(h, fit, e, param, len, 1, tc, TR, Runc);
    fits.canonical.info = info;
end
if any(models == "spline")
    [h, fit, e, param] = run_fit(@() Fit_Spline(tc, TR, Runc, window_seconds), suppress_warnings);
    fits.spline = package_fit(h, fit, e, param, len, 1, tc, TR, Runc);
end
if any(models == "nlgamma")
    [h, fit, e, param] = run_fit(@() Fit_NLgamma(tc, TR, Runc, window_seconds), suppress_warnings);
    fits.nlgamma = package_fit(h, fit, e, param, len, 1, tc, TR, Runc);
end
end

function models = local_filter_available_models(models, dependency_policy)
valid_models = ["logit", "sfir", "canonical", "spline", "nlgamma"];
models = models(:)';
unknown = setdiff(models, valid_models);
if ~isempty(unknown)
    error('Unknown HRF model(s): %s. Valid models are: %s.', ...
        strjoin(cellstr(unknown), ', '), strjoin(cellstr(valid_models), ', '));
end

keep = true(size(models));
missing_messages = strings(size(models));
for i = 1:numel(models)
    [ok, msg] = local_model_dependency(models(i));
    if ~ok
        keep(i) = false;
        missing_messages(i) = msg;
        switch dependency_policy
            case 'skip'
                local_warn_once(models(i), msg);
            case 'error'
                error('HRF model "%s" is unavailable: %s', models(i), msg);
        end
    end
end

if ~any(keep)
    msg = strjoin(cellstr(unique(missing_messages(missing_messages ~= ""))), ' ');
    error('No requested HRF models are available. %s', msg);
end

models = models(keep);
end

function [ok, msg] = local_model_dependency(model_name)
ok = true;
msg = "";
switch char(model_name)
    case 'canonical'
        if exist('spm_get_bf', 'file') ~= 2
            ok = false;
            msg = "canonical requires SPM on the MATLAB path (missing spm_get_bf).";
        end
    case 'nlgamma'
        if exist('spm_hrf', 'file') ~= 2
            ok = false;
            msg = "nlgamma requires SPM on the MATLAB path (missing spm_hrf).";
        end
    case 'spline'
        if exist('create_bspline_basis', 'file') ~= 2 || exist('eval_basis', 'file') ~= 2
            ok = false;
            msg = "spline requires the FDA package on the MATLAB path (missing create_bspline_basis/eval_basis; see https://github.com/markgewhite/fda).";
        end
end
end

function local_warn_once(model_name, msg)
persistent warned_models
if isempty(warned_models)
    warned_models = strings(0, 1);
end
key = string(model_name);
if any(warned_models == key)
    return
end
warned_models(end + 1, 1) = key;
warning('hrf_fit_all_models:SkippingModel', 'Skipping HRF model "%s": %s', char(key), char(msg));
end

function s = package_fit(h, fit, e, param, len, p, tc, TR, Runc)
s = struct();
s.hrf = h;
s.fit = fit;
s.residual = e;
s.param = param;
s.mse = (1 / (len - 1)) * sum(e .^ 2);
try
    s.mis_modeling_p = ResidScan(e, 4);
catch
    s.mis_modeling_p = NaN;
end
try
    s.power_loss = PowerLoss(e, fit, (len - p), tc, TR, Runc, 0.001);
catch
    s.power_loss = NaN;
end
end


function varargout = run_fit(funh, suppress_warnings)
if suppress_warnings
    ws = warning;
    warning('off', 'all');
    c = onCleanup(@() warning(ws));
    [varargout{1:nargout}] = funh();
else
    [varargout{1:nargout}] = funh();
end
end
