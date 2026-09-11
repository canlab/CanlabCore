function alloc = canlab_scan_time_allocation(hours, varargin)
% Allocate a fixed budget of scanner hours across candidate sample sizes (more subjects vs. more scan time per subject)
%
% Given a total number of scanner hours, this function works out, for each
% candidate number of subjects N, how many scanning sessions each subject
% would need, how much of each subject's time is usable functional
% scanning once per-session setup overhead is subtracted, and how many
% functional image volumes (TRs) that yields per subject. It is the
% "design tradeoff" engine used by canlab_effect_size_map and
% canlab_power_allocation_curves, but can also be used on its own to plan
% a study.
%
% The logic is:
%
%   hours_per_subject     = hours / N
%   sessions_per_subject  = ceil(hours_per_subject / hours_per_session)
%   functional_hours      = hours_per_subject - setup_min_first/60
%                           - setup_min_repeat/60 * (sessions_per_subject - 1)
%   images_per_subject    = functional_hours * 3600 / TR
%   effective_images      = images_per_subject * df_shrink_factor
%
% Candidate N for which effective_images < min_images (not enough data to
% estimate the first-level model) are dropped from the output.
%
% :Usage:
% ::
%
%     alloc = canlab_scan_time_allocation(hours, [optional inputs])
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2026 Tor Wager
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
% :Inputs:
%
%   **hours:**
%        Total scanner hours available for the whole study (scalar).
%
% :Optional Inputs:
%
%   **'N_range', [vector of integers]:**
%        Candidate numbers of subjects to evaluate. Default = 3:100.
%
%   **'hours_per_session', [scalar]:**
%        Maximum length of a single scanning session, in hours.
%        Default = 1.5.
%
%   **'TR', [scalar]:**
%        Repetition time in seconds (or an "effective TR" if you want to
%        count acquisitions differently). Default = 2.
%
%   **'setup_min_first', [scalar]:**
%        Minutes lost to setup, localizers, and structural scans in a
%        subject's first session. Default = 30.
%
%   **'setup_min_repeat', [scalar]:**
%        Minutes lost to setup in each additional session. Default = 15.
%
%   **'df_shrink_factor', [scalar in (0, 1]]:**
%        Ratio of effective (autocorrelation-corrected) degrees of freedom
%        to the number of image volumes. Used only for the minimum-data
%        feasibility check. The default of 0.84 comes from a TR = 2 s
%        dataset with an AR(2) noise model, where 215 images translated to
%        about 180 effective df (Wager et al., 2009). If you have subjects'
%        SPM.mat files, SPM.xX.erdf / sum(SPM.nscan) gives your own value.
%
%   **'min_images', [scalar]:**
%        Minimum number of effective images per subject for a first-level
%        model to be estimable. Candidate N that fall below this are
%        removed. A sensible value is the number of columns in the
%        first-level design matrix. Default = 30.
%
%   **'doplot', [logical flag]:**
%        Plot sessions and functional hours as a function of N.
%        Default = false. The keyword 'plot' also turns plotting on.
%
%   **'verbose', [logical flag]:**
%        Print a short table of the allocation. Default = false.
%        'noverbose' also turns this off.
%
% :Outputs:
%
%   **alloc:**
%        Structure with one entry per feasible candidate N:
%
%        - .hours                          total scanner hours (input)
%        - .N                              candidate numbers of subjects
%        - .hours_per_subject              scanner hours per subject
%        - .sessions_per_subject           number of sessions per subject
%        - .functional_hours_per_subject   usable functional scan time (hours)
%        - .functional_min_per_subject     same, in minutes
%        - .images_per_subject             functional image volumes per subject
%        - .effective_images_per_subject   images * df_shrink_factor
%        - .session_change_N               N values at which the number of
%                                          sessions per subject changes
%                                          (useful for marking plots)
%        - .infeasible_N                   candidate N that were dropped
%        - .parameters                     the option values used
%
% :Examples:
% ::
%
%    % 60 scanner hours, 1.5-hour sessions, TR = 2 s
%    alloc = canlab_scan_time_allocation(60, 'doplot', true, 'verbose', true);
%
%    % A design with 40 columns needs at least ~40 effective images/subject
%    alloc = canlab_scan_time_allocation(100, 'TR', 1, 'min_images', 40);
%    plot(alloc.N, alloc.images_per_subject); xlabel('N'); ylabel('Images per subject');
%
% :References:
%   Wager, T. D., et al. (2009). Brain mediators of cardiovascular
%   responses to social threat, Part I. NeuroImage, 47, 821-835.
%
% :See also:
%   - canlab_effect_size_map, canlab_power_allocation_curves,
%     power_from_variance
%

% ..
%    Programmers' notes:
%    2026-09: Tor Wager. Extracted from the legacy effect_size_map.m script
%    and converted to a documented function with inputParser.
% ..

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

% Parse special command keywords and remove them before inputParser

doplot = false;
plot_idx = strcmpi(varargin, 'plot');
if any(plot_idx)               % Override: omit 'doplot' key/value pair
    doplot = true;
    varargin(plot_idx) = [];   % remove so inputParser doesn't see it
end

verbose = false;
verbose_idx = strcmpi(varargin, 'noverbose');
if any(verbose_idx)
    verbose = false;
    varargin(verbose_idx) = [];   % remove so inputParser doesn't see it
end

% Use inputParser to parse key/value pairs

valfcn_posscalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', 'positive'});
valfcn_nonnegscalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', 'nonnegative'});
valfcn_logical = @(x) islogical(x) || isnumeric(x);

p = inputParser;
p.addRequired('hours', valfcn_posscalar);
p.addParameter('N_range', 3:100, @(x) validateattributes(x, {'numeric'}, {'nonempty', 'vector', 'positive', 'integer'}));
p.addParameter('hours_per_session', 1.5, valfcn_posscalar);
p.addParameter('TR', 2, valfcn_posscalar);
p.addParameter('setup_min_first', 30, valfcn_nonnegscalar);
p.addParameter('setup_min_repeat', 15, valfcn_nonnegscalar);
p.addParameter('df_shrink_factor', 0.84, @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', '>', 0, '<=', 1}));
p.addParameter('min_images', 30, valfcn_nonnegscalar);

% Special key/value pairs that we have potentially set with optional keywords
p.addParameter('doplot', doplot, valfcn_logical);
p.addParameter('verbose', verbose, valfcn_logical);

% process inputs
p.parse(hours, varargin{:});
ARGS = p.Results;

doplot = logical(ARGS.doplot);
verbose = logical(ARGS.verbose);

% -------------------------------------------------------------------------
% Allocation
% -------------------------------------------------------------------------

N = double(ARGS.N_range(:)');                       % row vector of candidate N
hours_per_subject = hours ./ N;

% How many sessions does each subject need to use their share of hours?
sessions_per_subject = ceil(hours_per_subject ./ ARGS.hours_per_session);

% Subtract setup overhead: a larger chunk in the first session
% (structural scans, localizers, instructions), a smaller chunk for
% each repeat session.
functional_hours = hours_per_subject ...
    - ARGS.setup_min_first / 60 ...
    - (ARGS.setup_min_repeat / 60) .* (sessions_per_subject - 1);

% Convert usable functional time to image volumes
images_per_subject = functional_hours .* 3600 ./ ARGS.TR;

% Effective (autocorrelation-corrected) images, for the feasibility check
effective_images = images_per_subject .* ARGS.df_shrink_factor;

% Remove candidate N where there is not enough time per subject to
% estimate a first-level model at all.
wh_infeasible = functional_hours <= 0 | effective_images < ARGS.min_images;

if all(wh_infeasible)
    error('canlab_scan_time_allocation:NoFeasibleN', ...
        'No candidate N is feasible: %3.1f hours is not enough to collect %3.0f effective images per subject for any N in N_range.', ...
        hours, ARGS.min_images);
end

infeasible_N = N(wh_infeasible);

N(wh_infeasible) = [];
hours_per_subject(wh_infeasible) = [];
sessions_per_subject(wh_infeasible) = [];
functional_hours(wh_infeasible) = [];
images_per_subject(wh_infeasible) = [];
effective_images(wh_infeasible) = [];

% N at which the number of sessions per subject changes (first N with the
% new session count). Handy for marking plots.
wh_change = find(diff(sessions_per_subject) ~= 0) + 1;
session_change_N = N(wh_change);

% -------------------------------------------------------------------------
% Output structure
% -------------------------------------------------------------------------

alloc = struct();
alloc.hours = hours;
alloc.N = N;
alloc.hours_per_subject = hours_per_subject;
alloc.sessions_per_subject = sessions_per_subject;
alloc.functional_hours_per_subject = functional_hours;
alloc.functional_min_per_subject = round(functional_hours .* 60);
alloc.images_per_subject = images_per_subject;
alloc.effective_images_per_subject = effective_images;
alloc.session_change_N = session_change_N;
alloc.infeasible_N = infeasible_N;

alloc.parameters = rmfield(ARGS, {'doplot', 'verbose'});

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------

if verbose
    fprintf('Allocation of %3.1f scanner hours (%3.1f-hour sessions, TR = %3.2f s)\n', ...
        hours, ARGS.hours_per_session, ARGS.TR);
    fprintf('%6s %10s %12s %14s %10s\n', 'N', 'Hrs/subj', 'Sessions', 'Func. min', 'Images');

    wh_show = unique([1 round(linspace(1, numel(N), min(numel(N), 10)))]);
    for i = wh_show
        fprintf('%6d %10.2f %12d %14.0f %10.0f\n', N(i), hours_per_subject(i), ...
            sessions_per_subject(i), alloc.functional_min_per_subject(i), images_per_subject(i));
    end

    if ~isempty(infeasible_N)
        fprintf('Dropped %d candidate N (too little functional time per subject): N >= %d\n', ...
            numel(infeasible_N), min(infeasible_N));
    end
end

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------

if doplot
    create_figure('Scan time allocation');

    plot(N, sessions_per_subject, 'k', 'LineWidth', 3);
    plot(N, functional_hours, 'k:', 'LineWidth', 3);
    legend({'Sessions per subject' 'Functional scan hours per subject'});
    axis tight
    plot_vertical_line(session_change_N);
    xlabel('Number of subjects');
    title(sprintf('Allocation of %3.0f scanner hours', hours));
end

end % main function
