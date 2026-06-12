function [Z, stats] = normalize_gm_shift_scale(Y, gm_mask, wm_mask, csf_mask, varargin)
% Normalize gray-matter voxel intensities across subjects by:
% (1) Removing a subject-specific additive shift estimated from CSF/WM medians
% (2) Correcting a subject-specific multiplicative scale estimated from MADs
%
% This function implements a compartment-based shift and scale normalization
% using robust statistics (medians and MAD) in CSF, WM, and GM.
%
% USAGE
%   [Z, stats] = normalize_gm_shift_scale(Y, gm_mask, wm_mask, csf_mask, ...
%                                         'do_scale', true, ...
%                                         'log_scale', false, 'trim_pct', 5);
%
% INPUTS
%   Y        : V x S numeric matrix of voxel intensities
%              - V = number of voxels
%              - S = number of subjects
%
%   gm_mask  : V x 1 logical vector; true for GM voxels
%   wm_mask  : V x 1 logical vector; true for WM voxels
%   csf_mask : V x 1 logical vector; true for CSF voxels
%
% OPTIONAL INPUTS (name/value pairs)
%   'do_scale'  : logical (default = true)
%                 - If true, estimate and apply multiplicative scale
%                   normalization.
%                 - If false, skip scale estimation and apply only the
%                   additive shift normalization.
%
%   'log_scale' : logical (default = false)
%                 - If true, fit scale model on log(MAD) values:
%                   log(r_GM) ~ log(r_CSF) + log(r_WM)
%                 - If false, fit linear model:
%                   r_GM ~ r_CSF + r_WM
%
%   'trim_pct'  : scalar (default = 5)
%                 - Percentage of values to trim at both lower and upper tails
%                   within each tissue for robust median/MAD estimation.
%                 - For example, trim_pct = 5 removes bottom and top 5%.
%
%   'use_gpu'   : logical (default = false)
%                 - If true, attempts to run the numeric kernel on a GPU
%                   using gpuArray. Results are gathered before return.
%                   If no GPU is available, or a GPU operation is not
%                   supported, the function warns once and falls back to CPU.
%
% OUTPUTS
%   Z     : V x S numeric matrix, shift- and scale-normalized image data
%           - Only GM voxels are modified. Non-GM voxels are copied from Y.
%
%   stats : structure with summary statistics and regression parameters
%       .m_gm, .m_wm, .m_csf   : S x 1 medians per subject/tissue
%       .r_gm, .r_wm, .r_csf   : S x 1 MAD-based scales per subject/tissue
%       .beta_mean             : 3 x 1 coefficients for:
%                                m_GM ~ [1, m_CSF, m_WM]
%       .mu_nuis               : S x 1 estimated nuisance shift for each subject
%       .do_scale              : logical flag indicating whether scale
%                                normalization was enabled
%       .beta_scale or .lambda : 3 x 1 coefficients for scale model
%                                (present only when do_scale is true)
%       .sigma_hat             : S x 1 fitted GM scale per subject
%       .sigma_ref             : scalar reference GM scale
%       .scale_factor          : S x 1 scaling factor (sigma_ref ./ sigma_hat),
%                                or ones(S,1) when do_scale is false
%       .log_scale             : logical flag indicating model type
%       .trim_pct              : trimming percentage used
%
% EXAMPLE
%   % Y: V x S matrix of beta images; gm_mask/wm_mask/csf_mask: V x 1 logical
%   [Z, stats] = normalize_gm_shift_scale(Y, gm_mask, wm_mask, csf_mask, ...
%                                         'log_scale', true, 'trim_pct', 5);
%
%   % Use Z(gm_mask, :) in downstream analyses
%
% NOTES
%   - This function assumes Y contains images for a single contrast/condition,
%     already coregistered across subjects.
%   - GM is assumed to contain signal of interest; CSF and WM are assumed to
%     contain no true task/contrast signal.
%   - Shift and scale parameters are estimated across subjects using tissue
%     summaries and then applied voxelwise to whole image
%   - values within GM are adjusted for individual differences in mean and scale of WM+CSF.
%
%   Author:  Tor Wager + ChatGPT5.2
%   Date:    2025-12-09
%

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

p = inputParser;
p.addRequired('Y', @(x) isnumeric(x) && ismatrix(x));
p.addRequired('gm_mask', @(x) islogical(x) && isvector(x));
p.addRequired('wm_mask', @(x) islogical(x) && isvector(x));
p.addRequired('csf_mask', @(x) islogical(x) && isvector(x));

[use_gpu_default, varargin] = parse_use_gpu_flag(varargin);

p.addParameter('do_scale', true, @(x) islogical(x) && isscalar(x));
p.addParameter('log_scale', false, @(x) islogical(x) && isscalar(x));
p.addParameter('trim_pct', 5, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 50);
p.addParameter('use_gpu', use_gpu_default, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));

p.parse(Y, gm_mask, wm_mask, csf_mask, varargin{:});

do_scale  = p.Results.do_scale;
log_scale = p.Results.log_scale;
trim_pct  = p.Results.trim_pct;
use_gpu   = logical(p.Results.use_gpu);

[Y, gm_mask, wm_mask, csf_mask] = check_and_resolve_dims(Y, gm_mask, wm_mask, csf_mask);

if use_gpu
    [gpu_available, gpu_message] = local_gpu_available();

    if gpu_available
        try
            [Z, stats] = normalize_gm_shift_scale_core(gpuArray(Y), gpuArray(gm_mask), ...
                gpuArray(wm_mask), gpuArray(csf_mask), do_scale, log_scale, trim_pct);
            Z = gather(Z);
            stats = gather_stats_struct(stats);
            return
        catch gpu_err
            warning('normalize_gm_shift_scale:GPUFallback', ...
                ['GPU path failed (%s). Falling back to CPU for ', ...
                 'normalize_gm_shift_scale.'], gpu_err.message);
        end
    else
        warning('normalize_gm_shift_scale:GPUUnavailable', ...
            'GPU requested, but no usable GPU is available (%s). Falling back to CPU.', gpu_message);
    end
end

[Z, stats] = normalize_gm_shift_scale_core(Y, gm_mask, wm_mask, csf_mask, ...
    do_scale, log_scale, trim_pct);

end % function normalize_gm_shift_scale


function [Z, stats] = normalize_gm_shift_scale_core(Y, gm_mask, wm_mask, csf_mask, do_scale, log_scale, trim_pct)

[~, S] = size(Y);

% -------------------------------------------------------------------------
% Preallocate summary statistics
% -------------------------------------------------------------------------

m_gm  = local_nan(S, 1, Y); % Median values per subject
m_wm  = local_nan(S, 1, Y);
m_csf = local_nan(S, 1, Y);

r_gm  = local_nan(S, 1, Y); % Spatial MAD per subject
r_wm  = local_nan(S, 1, Y);
r_csf = local_nan(S, 1, Y);

mad_scale_const = 1.4826;

% -------------------------------------------------------------------------
% Compute trimmed medians and MAD-based scales per tissue and subject
% -------------------------------------------------------------------------

for s = 1:S

    ys = Y(:, s);

    % GM
    [m_gm(s), r_gm(s)] = tissue_robust_summary(ys, gm_mask, trim_pct, mad_scale_const);

    % WM
    [m_wm(s), r_wm(s)] = tissue_robust_summary(ys, wm_mask, trim_pct, mad_scale_const);

    % CSF
    [m_csf(s), r_csf(s)] = tissue_robust_summary(ys, csf_mask, trim_pct, mad_scale_const);

end

% -------------------------------------------------------------------------
% Shift model: m_GM ~ 1 + m_CSF + m_WM
% -------------------------------------------------------------------------

X_mean = [local_ones(S, 1, Y), m_csf, m_wm];  % S x 3
y_mean = m_gm;                       % S x 1

beta_mean = X_mean \ y_mean;         % 3 x 1 OLS

% Nuisance shift component: gamma_1 * m_CSF + gamma_2 * m_WM
mu_nuis = beta_mean(2) * m_csf + beta_mean(3) * m_wm;   % S x 1, fitted mean gray matter per subject

% -------------------------------------------------------------------------
% Scale model: r_GM ~ r_CSF + r_WM  or  log(r_GM) ~ log(r_CSF) + log(r_WM)
% -------------------------------------------------------------------------

if do_scale
    eps_val = 1e-12;  % For log scale, avoid log(0)

    if ~log_scale

        % Linear scale model
        X_scale = [local_ones(S, 1, Y), r_csf, r_wm];
        y_scale = r_gm;

        beta_scale = X_scale \ y_scale;          % 3 x 1
        sigma_hat  = X_scale * beta_scale;       % fitted r_GM, S x 1

        scale_param = beta_scale;                % store under stats

    else

        % Log-scale model
        rc = max(r_csf, eps_val);
        rw = max(r_wm,  eps_val);
        rg = max(r_gm,  eps_val);

        X_log = [local_ones(S, 1, Y), log(rc), log(rw)];
        y_log = log(rg);

        lambda = X_log \ y_log;                  % 3 x 1
        log_sigma_hat = X_log * lambda;
        sigma_hat     = exp(log_sigma_hat);      % S x 1, positive by construction

        scale_param = lambda;                    % store under stats

    end

    % Reference scale: median fitted GM scale
    % Group median will be normalized to this value
    sigma_ref = median(sigma_hat(~isnan(sigma_hat) & isfinite(sigma_hat)));

    scale_factor = sigma_ref ./ sigma_hat;       % S x 1
else
    sigma_hat    = local_nan(S, 1, Y);
    sigma_ref    = NaN;
    scale_factor = local_ones(S, 1, Y);
end

% -------------------------------------------------------------------------
% Apply shift and scale normalization to GM voxels
% -------------------------------------------------------------------------

Z = Y;  % initialize output

% gm_idx = gm_mask(:);

for s = 1:S

    if local_scalar_true(isnan(scale_factor(s)) | ~isfinite(scale_factor(s)))
        % If something went wrong for this subject, leave as-is
        warning('normalize_gm_shift_scale: Scale factor is infinite for image %3.0f!', s)
        continue;
    end

    y_s = Y(:, s);

    % Shift-correct to adjust GM values (whole image)
    y_s = y_s - mu_nuis(s);

    % Scale-correct within GM
    y_s = scale_factor(s) * y_s;

    Z(:, s) = y_s;

end

% -------------------------------------------------------------------------
% Collect stats
% -------------------------------------------------------------------------

stats = struct();
stats.m_gm    = m_gm;
stats.m_wm    = m_wm;
stats.m_csf   = m_csf;
stats.r_gm    = r_gm;
stats.r_wm    = r_wm;
stats.r_csf   = r_csf;

stats.beta_mean   = beta_mean;
stats.mu_nuis     = mu_nuis;
stats.do_scale    = do_scale;

if do_scale && ~log_scale
    stats.beta_scale = scale_param;
elseif do_scale
    stats.lambda = scale_param;
end

stats.sigma_hat    = sigma_hat;
stats.sigma_ref    = sigma_ref;
stats.scale_factor = scale_factor;

stats.log_scale    = log_scale;
stats.trim_pct     = trim_pct;
stats.mad_const    = mad_scale_const;

end % function normalize_gm_shift_scale_core

% =========================================================================
% Helper functions
% =========================================================================

function [use_gpu_default, args] = parse_use_gpu_flag(args)
% Support a bare 'use_gpu' flag in addition to the name/value form.

use_gpu_default = false;

for i = 1:numel(args)
    if (ischar(args{i}) || isstring(args{i})) && strcmpi(char(args{i}), 'use_gpu')
        if i == numel(args) || ischar(args{i + 1}) || isstring(args{i + 1})
            use_gpu_default = true;
            args(i) = [];
            return
        end
    end
end

end


function [tf, message] = local_gpu_available()
% Check for a usable MATLAB GPU without making GPU support a hard dependency.

tf = false;
message = 'gpuDeviceCount/gpuDevice not available';

try
    if exist('gpuDeviceCount', 'file') ~= 2 && exist('gpuDeviceCount', 'builtin') ~= 5
        return
    end

    if gpuDeviceCount < 1
        message = 'gpuDeviceCount returned 0';
        return
    end

    gpuDevice;
    tf = true;
    message = '';

catch err
    message = err.message;
end

end


function stats = gather_stats_struct(stats)
% Gather every numeric/logical gpuArray field in the stats struct.

fn = fieldnames(stats);

for i = 1:numel(fn)
    if isa(stats.(fn{i}), 'gpuArray')
        stats.(fn{i}) = gather(stats.(fn{i}));
    end
end

end


function x = local_nan(m, n, exemplar)

if isa(exemplar, 'gpuArray')
    x = gpuArray(nan(m, n));
else
    x = nan(m, n);
end

end


function x = local_ones(m, n, exemplar)

if isa(exemplar, 'gpuArray')
    x = gpuArray(ones(m, n));
else
    x = ones(m, n);
end

end


function tf = local_scalar_true(x)

if isa(x, 'gpuArray')
    x = gather(x);
end

tf = logical(x);

end

function [m_t, r_t] = tissue_robust_summary(y, mask, trim_pct, mad_scale_const)
% Compute trimmed median and MAD-based robust scale for one tissue

vals = y(mask);

if isempty(vals)
    m_t = NaN;
    r_t = NaN;
    return;
end

% Trim extremes if requested
if trim_pct > 0
    lo = prctile(vals, trim_pct);
    hi = prctile(vals, 100 - trim_pct);
    vals = vals(vals >= lo & vals <= hi);
end

if isempty(vals)
    m_t = NaN;
    r_t = NaN;
    return;
end

m_t = median(vals);
mad_raw = median(abs(vals - m_t));
r_t = mad_scale_const * mad_raw;

end


function [Y_out, gm_mask_out, wm_mask_out, csf_mask_out] = check_and_resolve_dims(Y, gm_mask, wm_mask, csf_mask)
% Basic dimension checks and reshaping if needed

[V, S] = size(Y); %#ok<ASGLU>

gm_mask = gm_mask(:);
wm_mask = wm_mask(:);
csf_mask = csf_mask(:);

if length(gm_mask) ~= V || length(wm_mask) ~= V || length(csf_mask) ~= V
    error('Masks must have length equal to size(Y,1).');
end

% Logical checks
if ~islogical(gm_mask),  gm_mask  = gm_mask  ~= 0; end
if ~islogical(wm_mask),  wm_mask  = wm_mask  ~= 0; end
if ~islogical(csf_mask), csf_mask = csf_mask ~= 0; end

Y_out        = Y;
gm_mask_out  = gm_mask;
wm_mask_out  = wm_mask;
csf_mask_out = csf_mask;

end
