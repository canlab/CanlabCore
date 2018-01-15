function [X, e] = create_design_single_event(TR, ISI, eventduration, HPlength, dononlin)
% [X, e] = create_design_single_event(TR, ISI, eventduration, HPlength, dononlin)
%
% Create and plot design for single event
% Enter all values in sec:
% TR = time repetition for scans
% ISI = inter-stimulus interval
% eventduration = duration in sec of events
%
% HPlength is high-pass filter length in sec, or Inf or [] for no filter.

% ----------------------------------------------------------------
% * Fixed parameters
% ----------------------------------------------------------------
scanLength = 200; % in TRs
%dononlin = 0;
conditions = 1;
contrasts = 1;


if ~isempty(HPlength) && ~isinf(HPlength), dohpfilt = 1; else dohpfilt = 0; end

if dononlin
    nonlinstr = 'nonlinsaturation';
else
    nonlinstr = 'nononlin';
end

nconditions = length(conditions);
len = ceil(scanLength .* TR);
ons = cell(1, nconditions);

% there should be one column per condition in your design, *including* the
% intercept. One row per contrast.

% onsets and lengths in seconds
    ons{1} = [0:ISI:len]';
    ons{1}(:, 2) = eventduration;
    
    X = onsets2fmridesign(ons, TR, len, 'hrf', nonlinstr);

    % high-pass filtering
    if dohpfilt
        
        X(:, 1:end-1) = hpfilter(X(:, 1:end-1), TR, HPlength, scanLength);
        
    end
    
    plotDesign(ons, [], TR, 'samefig', 'durs', eventduration, nonlinstr, 'overlapping');
    
    set(gca, 'XLim', [0 scanLength * TR], 'XTick', round(linspace(0, scanLength * TR, 10)));
    
    % test efficiency
    
    e = calcEfficiency(ones(1, size(contrasts, 1)), [1 0], pinv(X), []);
    
    % related to (same as without contrasts, autocorr, filtering):
    % 1 ./ diag(inv(X' * X))

    
end % function




