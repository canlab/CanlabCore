function [xhat, info] = hrf_deconv_timeseries(y, hrf, varargin)
%HRF_DECONV_TIMESERIES HRF-informed deconvolution of BOLD timeseries -> neural proxy.
%
% Recovers a neural-activity proxy x(t) from a BOLD timeseries y(t) given a
% region/signature-specific hemodynamic kernel h (e.g. your estimated sFIR
% HRF), by inverting the convolution  y = h * x. This is the step that makes
% downstream Granger/lag analysis interpretable: the single largest confound
% for BOLD effective connectivity is that the HRF differs across regions
% (peak latency varies 1-4 s for purely vascular reasons), which can REVERSE
% the apparent direction of influence. Deconvolving each signal with its OWN
% measured kernel removes that confound before the causal analysis.
%
% Usage
% -----
%   xhat = hrf_deconv_timeseries(y, hrf)
%   xhat = hrf_deconv_timeseries(Y, H, 'RunLengths', [n1 n2 ...])
%   [xhat, info] = hrf_deconv_timeseries(Y, H, 'Method','wiener', 'NSR',0.1)
%
% Inputs
% ------
%   y    - [T x N] timeseries, columns = signals (regions/signatures).
%   hrf  - the hemodynamic kernel(s), sampled at the SAME dt as y (your sFIR
%          lags are already at TR spacing). Either [L x 1] (one kernel applied
%          to every column) or [L x N] (a per-column kernel). NaNs/trailing
%          zeros are trimmed per column.
%
% Optional (name-value)
% ---------------------
%   'Method'     - 'ridge' (default, SVD Tikhonov with GCV-chosen lambda,
%                  parameter-free and robust for short runs) or 'wiener'
%                  (FFT Wiener filter, fast; uses 'NSR').
%   'Lambda'     - ridge: [] (default) => choose by GCV; or a scalar to fix it.
%   'NSR'        - wiener noise-to-signal ratio (default 0.1).
%   'RunLengths' - vector summing to T; deconvolve each run block separately
%                  (the HRF does not span run boundaries). Default [] => one
%                  block.
%   'Normalize'  - peak-normalize each kernel to 1 before solving (default
%                  true). Deconvolution scale is absorbed into x, so this only
%                  sets x's units; GC downstream is scale-invariant.
%   'Detrend'    - linearly detrend each run block of y first (default true).
%
% Outputs
% -------
%   xhat - [T x N] neural proxy, same shape as y.
%   info - struct: .method, .lambda [N x nRuns] (ridge), .nsr, .run_lengths,
%          .burn_in (samples per run dominated by the deconvolution transient;
%          = kernel length, useful to trim before GC).
%
% Notes
% -----
% * 'ridge' builds the T x T lower-triangular Toeplitz convolution matrix and
%   solves via SVD, so lambda is chosen by generalized cross-validation with
%   no free parameter. O(T^3) per (signal,run) -- fine at region/signature
%   level (T a few hundred); use 'wiener' for voxelwise scale.
% * Feed a CLEAN kernel: the group-mean sFIR HRF (hrf_misspec_metrics
%   GroupCurveFirst, or hrf_curve_summaries) is far less noisy than a single
%   subject/run curve.
%
% Example
% -------
%   TR = 0.46; t = (0:TR:30)';
%   h  = spm_hrf(TR);                       % a kernel at TR resolution
%   x  = double(rand(400,1) > 0.9);         % sparse "neural" events
%   y  = conv(x, h); y = y(1:400) + 0.05*randn(400,1);
%   xhat = hrf_deconv_timeseries(y, h, 'TR'); %#ok<NASGU>
%
% See also: hrf_granger_causality, hrf_misspec_metrics, hrf_curve_summaries.

p = inputParser;
p.addRequired('y', @(x) isnumeric(x) && ~isempty(x));
p.addRequired('hrf', @(x) isnumeric(x) && ~isempty(x));
p.addParameter('Method', 'ridge', @(x) ischar(x) || isstring(x));
p.addParameter('Lambda', [], @(x) isempty(x) || isscalar(x));
p.addParameter('NSR', 0.1, @(x) isscalar(x) && x > 0);
p.addParameter('RunLengths', [], @(x) isempty(x) || isvector(x));
p.addParameter('Normalize', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Detrend', true, @(x) islogical(x) || isnumeric(x));
p.parse(y, hrf, varargin{:});
opts = p.Results;

if isvector(y), y = y(:); end
[T, N] = size(y);

% one kernel for all columns, or one per column
if isvector(hrf), hrf = hrf(:); end
if size(hrf, 2) == 1 && N > 1, hrf = repmat(hrf, 1, N); end
if size(hrf, 2) ~= N
    error('hrf_deconv_timeseries:KernelCols', ...
        'hrf must have 1 or N=%d columns; got %d.', N, size(hrf, 2));
end

run_len = opts.RunLengths;
if isempty(run_len), run_len = T; end
run_len = run_len(:)';
if sum(run_len) ~= T
    error('hrf_deconv_timeseries:RunLengths', ...
        'RunLengths sum (%d) must equal T (%d).', sum(run_len), T);
end
run_edges = [0, cumsum(run_len)];
nruns = numel(run_len);

method = lower(char(opts.Method));
xhat = zeros(T, N);
lambda_used = nan(N, nruns);
max_klen = 0;

for n = 1:N
    k = local_clean_kernel(hrf(:, n), opts.Normalize);
    max_klen = max(max_klen, numel(k));
    for r = 1:nruns
        idx = (run_edges(r) + 1):run_edges(r + 1);
        yr = y(idx, n);
        if opts.Detrend, yr = detrend(yr, 1); end
        switch method
            case 'ridge'
                [xr, lam] = local_ridge_deconv(yr, k, opts.Lambda);
            case 'wiener'
                xr = local_wiener_deconv(yr, k, opts.NSR); lam = NaN;
            otherwise
                error('hrf_deconv_timeseries:Method', 'Unknown Method: %s', method);
        end
        xhat(idx, n) = xr;
        lambda_used(n, r) = lam;
    end
end

info = struct('method', method, 'lambda', lambda_used, 'nsr', opts.NSR, ...
    'run_lengths', run_len, 'burn_in', max_klen);
end


% =========================================================================
function k = local_clean_kernel(k, do_norm)
k = k(:);
k(~isfinite(k)) = 0;
% drop trailing zeros (sFIR curves are zero-padded past the window)
last = find(k ~= 0, 1, 'last');
if isempty(last), error('hrf_deconv_timeseries:ZeroKernel', 'A kernel is all zeros.'); end
k = k(1:last);
if do_norm
    pk = max(abs(k));
    if pk > 0, k = k / pk; end
end
end


function [x, lambda] = local_ridge_deconv(y, k, fixed_lambda)
% Tikhonov deconvolution of y = H x via SVD, GCV-chosen lambda.
T = numel(y);
L = min(numel(k), T);
k = k(1:L);
H = toeplitz([k; zeros(T - L, 1)], [k(1); zeros(T - 1, 1)]);  % lower-tri Toeplitz
[U, S, V] = svd(H, 'econ');
s = diag(S);
uty = U' * y;

if ~isempty(fixed_lambda)
    lambda = fixed_lambda;
else
    lambda = local_gcv_lambda(s, uty, T);
end
filt = s ./ (s.^2 + lambda);
x = V * (filt .* uty);
end


function lambda = local_gcv_lambda(s, uty, T)
% Minimize GCV(l) = T*RSS(l) / (T - df(l))^2 over a log grid.
s2 = s.^2;
g = logspace(-6, 3, 60) * (mean(s2) + eps);
best = inf; lambda = g(1);
for i = 1:numel(g)
    l = g(i);
    df = sum(s2 ./ (s2 + l));
    % residual energy: components attenuated by l/(s^2+l)
    rss = sum(((l ./ (s2 + l)) .* uty).^2);
    denom = (T - df)^2;
    if denom <= 0, continue; end
    gcv = T * rss / denom;
    if gcv < best, best = gcv; lambda = l; end
end
end


function x = local_wiener_deconv(y, k, nsr)
% FFT Wiener filter: X = conj(K).*Y ./ (|K|^2 + nsr).
T = numel(y);
K = fft(k, T);
Y = fft(y, T);
X = (conj(K) .* Y) ./ (abs(K).^2 + nsr);
x = real(ifft(X));
end
