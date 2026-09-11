function out = canlab_power_allocation_curves(varargin)
% Normative power curves for allocating scanner hours between subjects and scan time, as a function of within- and between-subject variance
%
% This is a "what if" calculator that needs no data. For a grid of
% within-subject noise levels (per image volume) and between-subject
% standard deviations, and a fixed budget of scanner hours, it computes
% the expected group t-statistic and power for every candidate number of
% subjects N, and reports the N (and minutes of scanning per subject) that
% maximizes power. The same time-allocation model as canlab_effect_size_map
% is used (see canlab_scan_time_allocation), so results are comparable to
% the data-driven analysis; use the within- and between-subject standard
% deviations printed by canlab_effect_size_map to locate your study on the
% grid.
%
% Model: for a contrast of magnitude effect, N subjects, n images per
% subject, within-subject std per image sigma_w, and between-subject std
% sigma_b, the expected one-sample t-statistic is
%
%     t_expected = effect * sqrt(N) / sqrt(sigma_b^2 + sigma_w^2 / n)
%
% and power at a two-tailed alpha is 1 - nctcdf(tinv(1 - alpha/2, N-1), N-1, t_expected).
%
% :Usage:
% ::
%
%     out = canlab_power_allocation_curves([optional inputs])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2010, 2026 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Optional Inputs:
%
%   **'within_std', [vector]:**
%        Within-subject standard deviation per image volume (residual std
%        times the square root of the per-scan design variance; the
%        "std_w" reported by canlab_effect_size_map). Default = 12:72.
%
%   **'between_std', [vector]:**
%        Between-subject standard deviation of the contrast ("std_b").
%        Default = 1:0.1:10.
%
%   **'effect', [scalar]:**
%        Contrast magnitude (mean effect across subjects). Default = 1.
%        Only the ratios effect/std matter, so you can set effect = 1 and
%        express the stds in units of the effect.
%
%   **'alpha', [scalar]:**
%        Two-tailed uncorrected alpha level. Default = 0.001.
%
%   **'hours', 'hours_per_session', 'TR', 'N_range', 'setup_min_first', 'setup_min_repeat', 'min_images', 'df_shrink_factor':**
%        Passed to canlab_scan_time_allocation; see its help. Defaults:
%        60 hours, 1.5 hours/session, TR = 2, N_range = 3:100, 30 and 15
%        minutes setup, min_images = 30, df_shrink_factor = 0.84.
%
%   **'power_method', ['noncentral' | 'shift']:**
%        'noncentral' (default) uses the noncentral t distribution;
%        'shift' uses the faster approximation 1 - tcdf(u - t_expected, df).
%
%   **'doplot', [logical flag]:**
%        Line plots of power vs. N and contour plots of the optimal N and
%        minutes per subject. Default = true. 'noplot' to turn off.
%
%   **'verbose', [logical flag]:**
%        Print a short summary. Default = true. 'noverbose' to turn off.
%
% :Outputs:
%
%   **out:**
%        Structure with fields:
%
%        - .alloc                   scan-time allocation per candidate N
%        - .within_std, .between_std, .effect, .alpha
%        - .expected_t              [within x N x between] expected t
%        - .power                   [within x N x between] power
%        - .max_power               [within x between] best power over N
%        - .optimal_N               [within x between] N with the best power
%        - .optimal_min_per_subject [within x between] functional minutes
%                                   per subject at the optimal N
%
% :Examples:
% ::
%
%    % Default normative grid, 60 scanner hours
%    out = canlab_power_allocation_curves;
%
%    % A quick look at a few noise levels with 100 hours and a 1-s TR
%    out = canlab_power_allocation_curves('within_std', [20 40 60], ...
%        'between_std', 1:5, 'hours', 100, 'TR', 1);
%    disp(out.optimal_N)
%
%    % Find the optimal N for your own study's variance components
%    % (std_w and std_b as printed by canlab_effect_size_map)
%    out = canlab_power_allocation_curves('within_std', 35, 'between_std', 2.5, ...
%        'effect', 1.2, 'doplot', false);
%    fprintf('Optimal N = %d, %d min of scanning per subject\n', out.optimal_N, out.optimal_min_per_subject);
%
% :References:
%   Mumford, J. A. & Nichols, T. E. (2008). Power calculation for group
%   fMRI studies accounting for arbitrary design and temporal
%   autocorrelation. NeuroImage, 39, 261-268.
%
% :See also:
%   - canlab_effect_size_map, canlab_scan_time_allocation, power_from_variance
%

% ..
%    Programmers' notes:
%    2010: Tor Wager. Normative section of the legacy effect_size_map.m
%          script (used for figures in teaching materials).
%    2026-09: Tor Wager. Converted to a function with inputParser; computes
%          the whole within x between grid per N in one vectorized step;
%          replaced deprecated legend(..., 7) syntax; noncentral t power.
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

% Parse special command keywords and remove them before inputParser

doplot = true;
plot_idx = strcmpi(varargin, 'noplot');
if any(plot_idx)
    doplot = false;
    varargin(plot_idx) = [];   % remove so inputParser doesn't see it
end
plot_idx = strcmpi(varargin, 'plot');
if any(plot_idx)               % Override: omit 'doplot' key/value pair
    doplot = true;
    varargin(plot_idx) = [];
end

verbose = true;
verbose_idx = strcmpi(varargin, 'noverbose');
if any(verbose_idx)
    verbose = false;
    varargin(verbose_idx) = [];   % remove so inputParser doesn't see it
end

% Use inputParser to parse key/value pairs

valfcn_posvector = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'vector', 'positive'});
valfcn_posscalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', 'positive'});
valfcn_logical = @(x) islogical(x) || isnumeric(x);

p = inputParser;
p.addParameter('within_std', 12:72, valfcn_posvector);
p.addParameter('between_std', 1:.1:10, valfcn_posvector);
p.addParameter('effect', 1, valfcn_posscalar);
p.addParameter('alpha', 0.001, @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0, '<', 1}));

p.addParameter('hours', 60, valfcn_posscalar);
p.addParameter('hours_per_session', 1.5, valfcn_posscalar);
p.addParameter('TR', 2, valfcn_posscalar);
p.addParameter('N_range', 3:100, @(x) validateattributes(x, {'numeric'}, {'nonempty', 'vector', 'positive', 'integer'}));
p.addParameter('setup_min_first', 30, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p.addParameter('setup_min_repeat', 15, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p.addParameter('min_images', 30, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
p.addParameter('df_shrink_factor', 0.84, @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0, '<=', 1}));
p.addParameter('power_method', 'noncentral', @(x) any(strcmpi(x, {'noncentral', 'shift'})));

% Special key/value pairs that we have potentially set with optional keywords
p.addParameter('doplot', doplot, valfcn_logical);
p.addParameter('verbose', verbose, valfcn_logical);

% process inputs
p.parse(varargin{:});
ARGS = p.Results;

doplot = logical(ARGS.doplot);
verbose = logical(ARGS.verbose);
power_method = lower(char(ARGS.power_method));

within_std = double(ARGS.within_std(:));      % column: rows of the grid
between_std = double(ARGS.between_std(:)');   % row: columns of the grid
effect = ARGS.effect;
alpha = ARGS.alpha;

% -------------------------------------------------------------------------
% Scan-time allocation for each candidate N
% -------------------------------------------------------------------------

alloc = canlab_scan_time_allocation(ARGS.hours, 'N_range', ARGS.N_range, ...
    'hours_per_session', ARGS.hours_per_session, 'TR', ARGS.TR, ...
    'setup_min_first', ARGS.setup_min_first, 'setup_min_repeat', ARGS.setup_min_repeat, ...
    'df_shrink_factor', ARGS.df_shrink_factor, 'min_images', ARGS.min_images);

nW = numel(within_std);
nB = numel(between_std);
nN = numel(alloc.N);

% -------------------------------------------------------------------------
% Expected t and power on the within x N x between grid
% -------------------------------------------------------------------------

[power, expected_t] = deal(zeros(nW, nN, nB));

for i = 1:nN
    N = alloc.N(i);
    df = N - 1;

    % Within-subject contribution to the group-level variance for this
    % many images per subject: [nW x 1]; between: [1 x nB]
    sig2wi = within_std .^ 2 ./ alloc.images_per_subject(i);
    sig2b = between_std .^ 2;

    ncp = effect .* sqrt(N) ./ sqrt(sig2b + sig2wi);       % [nW x nB], implicit expansion
    u = tinv(1 - alpha / 2, df);                            % two-tailed critical t

    switch power_method
        case 'noncentral'
            pw = 1 - nctcdf(u, df, ncp);
        case 'shift'
            pw = tcdf(u - ncp, df, 'upper');
    end

    expected_t(:, i, :) = reshape(ncp, nW, 1, nB);
    power(:, i, :) = reshape(pw, nW, 1, nB);
end

% Best N for each (within, between) combination
[max_power, wh_best] = max(power, [], 2);
max_power = reshape(max_power, nW, nB);
wh_best = reshape(wh_best, nW, nB);

optimal_N = reshape(alloc.N(wh_best), nW, nB);
optimal_min = reshape(alloc.functional_min_per_subject(wh_best), nW, nB);

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------

out = struct();
out.alloc = alloc;
out.within_std = within_std;
out.between_std = between_std;
out.effect = effect;
out.alpha = alpha;
out.power_method = power_method;
out.expected_t = expected_t;
out.power = power;
out.max_power = max_power;
out.optimal_N = optimal_N;
out.optimal_min_per_subject = optimal_min;

if verbose
    fprintf('Normative power curves: %3.0f scanner hours, effect = %3.2f, alpha = %g (two-tailed)\n', ARGS.hours, effect, alpha);
    fprintf('Within-subject std per image: %d levels from %3.1f to %3.1f\n', nW, min(within_std), max(within_std));
    fprintf('Between-subject std:          %d levels from %3.1f to %3.1f\n', nB, min(between_std), max(between_std));
    fprintf('Candidate N: %d to %d (%d values)\n', min(alloc.N), max(alloc.N), nN);
    fprintf('Optimal N ranges from %d to %d across the grid; best power from %3.2f to %3.2f\n', ...
        min(optimal_N(:)), max(optimal_N(:)), min(max_power(:)), max(max_power(:)));
end

% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------

if doplot
    plot_power_lines(out);
    if nW > 1 && nB > 1
        plot_optimal_contours(out);
    end
end

end % main function



% =========================================================================
% Subfunctions
% =========================================================================

function plot_power_lines(out)
% Power vs. N at low / middle / high within-subject noise, one line per
% between-subject std level (up to 10 levels, light gray = high sigma_b)

alloc = out.alloc;
nW = numel(out.within_std);
nB = numel(out.between_std);

wi_index = unique([1 round(nW / 2) nW]);
wh_b = unique(round(linspace(1, nB, min(nB, 10))));
grays = linspace(0, .7, numel(wh_b));

create_figure('Normative power curves', 1, numel(wi_index));

for wi = 1:numel(wi_index)
    subplot(1, numel(wi_index), wi);

    legstr = cell(1, numel(wh_b));
    for b = 1:numel(wh_b)
        plot(alloc.N, squeeze(out.power(wi_index(wi), :, wh_b(b))), 'Color', grays(b) * [1 1 1], 'LineWidth', 2);
        legstr{b} = sprintf('\\sigma_B = %3.1f', out.between_std(wh_b(b)));
    end

    axis tight
    set(gca, 'YLim', [0 1]);
    plot_vertical_line(alloc.session_change_N);

    if wi == numel(wi_index)
        legend(legstr, 'Location', 'best');
    end

    if wi == round(numel(wi_index) / 2)
        title(sprintf('Allocation of %3.0f scan hours\n\\sigma_w = %3.0f', alloc.hours, out.within_std(wi_index(wi))));
        xlabel('Number of subjects');
    else
        title(sprintf('\\sigma_w = %3.0f', out.within_std(wi_index(wi))));
    end
    if wi == 1, ylabel('Power'); end
end
end


function plot_optimal_contours(out)
% Contour maps of the optimal N and minutes per subject over the
% (between std, within std) grid

[X, Y] = meshgrid(out.between_std, out.within_std);

create_figure('Optimal allocation contours', 1, 2);

subplot(1, 2, 1);
[C, H] = contourf(X, Y, out.optimal_N);
axis tight
hh = clabel(C, H);
set(hh, 'FontWeight', 'bold');
cm = colormap(gray);
colormap(cm(15:end, :));
xlabel('Between-subjects error: \sigma_B');
ylabel('Within-subjects error: \sigma_w');
title('Optimal sample size (N)');

subplot(1, 2, 2);
[C, H] = contourf(X, Y, out.optimal_min_per_subject);
axis tight
hh = clabel(C, H);
set(hh, 'FontWeight', 'bold');
xlabel('Between-subjects error: \sigma_B');
title('Minutes of scanning per subject');
end
