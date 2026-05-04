function fits = hrf_fit_all_models(tc, TR, Runc, window_seconds, models, varargin)
%HRF_FIT_ALL_MODELS Fit supported HRF models with a consistent output API.

% Optional name/value:
%   'SuppressWarnings' (default true): temporarily suppress warning spam
p = inputParser;
p.addParameter('SuppressWarnings', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
suppress_warnings = logical(p.Results.SuppressWarnings);

models = lower(string(models));
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
    [varargout{1:nargout}] = funh(); %#ok<NASGU>
else
    [varargout{1:nargout}] = funh();
end
end
