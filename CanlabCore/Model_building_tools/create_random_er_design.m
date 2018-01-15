function [X, e, ons] = create_random_er_design(TR, ISI, eventduration, freqConditions, HPlength, dononlin)
% [X, e, ons] = create_random_er_design(TR, ISI, eventduration, freqConditions, HPlength, dononlin)
%
% Create and plot design for single event
% Enter all values in sec:
% TR = time repetition for scans
% ISI = inter-stimulus interval
% eventduration = duration in sec of events
% freqConditions = vector of frequencies of each event type, e.g. [.2 .2 .2] for 3 events at 20% each (remainder is rest)
%   - do not have to sum to one
%
% HPlength is high-pass filter length in sec, or Inf or [] for no filter.
%
% Examples:
% create_figure;
% [X, e] = create_random_er_design(1, 1.3, 1, [.2 .2], 180, 0);
% axis tight


% ----------------------------------------------------------------
% * Fixed parameters
% ----------------------------------------------------------------
scanLength = 200; % in sec

conditions = 1:length(freqConditions);

if ~isempty(HPlength) && ~isinf(HPlength), dohpfilt = 1; else dohpfilt = 0; end

if dononlin
    nonlinstr = 'nonlinsaturation';
else
    nonlinstr = 'nononlin';
end

nconditions = length(conditions);
len = ceil(scanLength .* TR);
ons = cell(1, nconditions);


% Create onsets
% ----------------------------------------------------------------

num_stim_slots = scanLength ./ ISI;
possible_onset_times = [0:ISI:(num_stim_slots - 1)]';
total_n_events = length(possible_onset_times);

for i = 1:nconditions
    
    n = length(possible_onset_times);
    
    nevents = floor(freqConditions(i) .* total_n_events);
    
    rnd = randperm(n);
    
    ons{i} = sort(possible_onset_times(rnd(1:nevents)));
    
    ons{i}(:, 2) = eventduration;
    
    % remove these onset times from list of available onset times
    possible_onset_times(rnd(1:nevents)) = [];
    
end

% Create contrasts
% ----------------------------------------------------------------

% there should be one column per condition in your design, *including* the
% intercept. One row per contrast.

contrasts = create_orthogonal_contrast_set(nconditions);
contrasts(:, end+1) = 0; % for intercept

% there should be one column per condition in your design, *including* the
% intercept. One row per contrast.



% Build Design Matrix
% ----------------------------------------------------------------


X = onsets2fmridesign(ons, TR, len, 'hrf', nonlinstr);

% high-pass filtering
if dohpfilt
    
    X(:, 1:end-1) = hpfilter(X(:, 1:end-1), TR, HPlength, scanLength);
    
end

plotDesign(ons, [], TR, 'samefig', 'durs', eventduration, nonlinstr);

set(gca, 'XLim', [0 scanLength], 'XTick', round(linspace(0, scanLength, 10)));

% Test efficiency
% ----------------------------------------------------------------

e = calcEfficiency(ones(1, size(contrasts, 1)), contrasts, pinv(X), []);

% related to (same as without contrasts, autocorr, filtering):
% 1 ./ diag(inv(X' * X))


end % function




