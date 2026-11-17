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
    fits.logit.uncertainty_source = 'not available for nonlinear logit fit';
end
if any(models == "fir")
    [h, fit, e, param] = run_fit(@() Fit_sFIR(tc, TR, Runc, window_seconds, 0), suppress_warnings);
    fits.fir = package_fit(h, fit, e, param, len, window_seconds, tc, TR, Runc);
    fits.fir = add_linear_uncertainty(fits.fir, 'fir', tc, TR, Runc, window_seconds);
end
if any(models == "sfir")
    [h, fit, e, param] = run_fit(@() Fit_sFIR(tc, TR, Runc, window_seconds, 1), suppress_warnings);
    fits.sfir = package_fit(h, fit, e, param, len, window_seconds, tc, TR, Runc);
    fits.sfir = add_linear_uncertainty(fits.sfir, 'sfir', tc, TR, Runc, window_seconds);
end
if any(models == "canonical")
    [h, fit, e, param, info] = run_fit(@() Fit_Canonical_HRF(tc, TR, Runc, window_seconds, 1), suppress_warnings);
    fits.canonical = package_fit(h, fit, e, param, len, 1, tc, TR, Runc);
    fits.canonical.info = info;
    fits.canonical = add_linear_uncertainty(fits.canonical, 'canonical', tc, TR, Runc, window_seconds);
end
if any(models == "spline")
    [h, fit, e, param] = run_fit(@() Fit_Spline(tc, TR, Runc, window_seconds), suppress_warnings);
    fits.spline = package_fit(h, fit, e, param, len, 1, tc, TR, Runc);
    fits.spline = add_linear_uncertainty(fits.spline, 'spline', tc, TR, Runc, window_seconds);
end
if any(models == "nlgamma")
    [h, fit, e, param] = run_fit(@() Fit_NLgamma(tc, TR, Runc, window_seconds), suppress_warnings);
    fits.nlgamma = package_fit(h, fit, e, param, len, 1, tc, TR, Runc);
    fits.nlgamma.uncertainty_source = 'not available for nonlinear gamma fit';
end
end

function models = local_filter_available_models(models, dependency_policy)
valid_models = ["logit", "fir", "sfir", "canonical", "spline", "nlgamma"];
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
s.time = local_time_vector(size(h, 1), TR);
s.fit = fit;
s.residual = e;
s.param = param;
s.N = len;
s.dfe = [];
s.se = [];
s.t = [];
s.p = [];
s.p_type = '';
s.uncertainty_source = '';
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

function s = add_linear_uncertainty(s, model_name, tc, TR, Runc, window_seconds)
try
    switch model_name
        case 'fir'
            [X, PX, hrf_lift, coef_idx_by_condition] = local_fir_uncertainty_design(Runc, TR, window_seconds);
        case 'sfir'
            [X, PX, hrf_lift, coef_idx_by_condition] = local_sfir_uncertainty_design(Runc, TR, window_seconds);
        case 'canonical'
            [X, PX, hrf_lift, coef_idx_by_condition] = local_canonical_uncertainty_design(Runc, TR, size(s.hrf, 1));
        case 'spline'
            [X, PX, hrf_lift, coef_idx_by_condition] = local_spline_uncertainty_design(Runc, TR, window_seconds, size(s.hrf, 1));
        otherwise
            return
    end

    residual = tc(:) - X * (PX * tc(:));
    dfe = max(size(X, 1) - trace(X * PX), 1);
    mse = sum(residual .^ 2) ./ dfe;
    covb = mse .* (PX * PX');

    n_cond = numel(coef_idx_by_condition);
    n_time = size(s.hrf, 1);
    se = nan(n_time, n_cond);
    for c = 1:n_cond
        idx = coef_idx_by_condition{c};
        covc = covb(idx, idx);
        se(:, c) = sqrt(max(diag(hrf_lift * covc * hrf_lift'), 0));
    end

    tval = s.hrf ./ se;
    pval = 2 * (1 - tcdf(abs(tval), dfe));
    pval(pval == 0) = eps;

    s.se = se;
    s.t = tval;
    s.p = pval;
    s.p_type = sprintf('Two-tailed P-values from %s HRF coefficient SE, dfe = %.3f', model_name, dfe);
    s.dfe = dfe;
    s.N = size(X, 1);
    s.mse = mse;
    s.uncertainty_source = sprintf('%s linear model residual covariance', model_name);
catch err
    s.se = [];
    s.t = [];
    s.p = [];
    s.p_type = '';
    s.dfe = [];
    s.uncertainty_source = sprintf('unavailable: %s', err.message);
end
end

function [X, PX, hrf_lift, coef_idx_by_condition] = local_sfir_uncertainty_design(Runc, TR, window_seconds)
[X, ~, hrf_lift, coef_idx_by_condition] = local_fir_uncertainty_design(Runc, TR, window_seconds);
numstim = numel(Runc);
tlen = size(hrf_lift, 1);

C = (1:tlen)' * ones(1, tlen);
h = sqrt(1 / (7 / TR));
v = 0.1;
sig = 1;
R = v * exp(-h / 2 * (C - C') .^ 2);
RI = R \ eye(size(R));
pen = zeros(numstim * tlen + 1);
for i = 1:numstim
    idx = ((i - 1) * tlen + 1):(i * tlen);
    pen(idx, idx) = sig ^ 2 * RI;
end

PX = (X' * X + pen) \ X';
end

function [X, PX, hrf_lift, coef_idx_by_condition] = local_fir_uncertainty_design(Runc, TR, window_seconds)
numstim = numel(Runc);
t = 1:TR:window_seconds;
tlen = numel(t);
len = numel(Runc{1});
Runs = zeros(len, numstim);
for i = 1:numstim
    Runs(:, i) = Runc{i}(:);
end
X = tor_make_deconv_mtx3(Runs, tlen, 1);
PX = pinv(X);
hrf_lift = eye(tlen);
coef_idx_by_condition = cell(1, numstim);
for i = 1:numstim
    coef_idx_by_condition{i} = ((i - 1) * tlen + 1):(i * tlen);
end
end

function [X, PX, hrf_lift, coef_idx_by_condition] = local_canonical_uncertainty_design(Runc, TR, n_hrf)
numstim = numel(Runc);
len = numel(Runc{1});
h = local_canonical_basis(TR, n_hrf);
Xtask = zeros(len, numstim);
for i = 1:numstim
    v = conv(Runc{i}, h);
    Xtask(:, i) = v(1:len);
end
X = [ones(len, 1) Xtask];
PX = pinv(X);
hrf_lift = h(:);
coef_idx_by_condition = cell(1, numstim);
for i = 1:numstim
    coef_idx_by_condition{i} = i + 1;
end
end

function h = local_canonical_basis(TR, n_hrf)
len = max(round(30 / TR), n_hrf);
xBF.dt = TR;
xBF.length = len;
xBF.name = 'hrf (with time and dispersion derivatives)';
xBF = spm_get_bf(xBF);
h = xBF.bf(1:n_hrf, 1);
h = h ./ max(h);
end

function [X, PX, hrf_lift, coef_idx_by_condition] = local_spline_uncertainty_design(Runc, TR, window_seconds, n_hrf)
numstim = numel(Runc);
len = numel(Runc{1});
t = 1:TR:window_seconds;
tlen = numel(t);
K = 8;
norder = 4;
basis = create_bspline_basis([0, tlen], K + 3, norder);
B = eval_basis((1:tlen), basis);
B = B(:, 3:end-1);
B = B(1:n_hrf, :);

Wi = zeros(len, numstim * K);
for j = 1:numstim
    Wji = tor_make_deconv_mtx3(Runc{j}, tlen, 1);
    Wi(:, (j - 1) * K + 1:j * K) = Wji(:, 1:tlen) * B;
end
X = [ones(len, 1) Wi];
PX = pinv(X);
hrf_lift = B;
coef_idx_by_condition = cell(1, numstim);
for i = 1:numstim
    coef_idx_by_condition{i} = ((i - 1) * K + 2):(i * K + 1);
end
end

function t = local_time_vector(n, TR)
t = 1 + (0:n - 1)' .* TR;
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
