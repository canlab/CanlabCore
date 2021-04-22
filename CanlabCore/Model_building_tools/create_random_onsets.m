function ons = create_random_onsets(scanLength, ISI, freqConditions, eventduration)
% ons = create_random_onsets(scanLength, ISI, freqConditions, eventduration)
%
% scanLength: length in seconds
% ISI: min inter-stimulus interval in seconds
% freqConditions: vector of probability values for each condition, 
%    e.g.,[.4 .4]. Should sum to <= 1.  If <1, rest intervals will be
%    included.
% eventduration: duration in seconds
%
% ons: Cell array of onset times for each condition in seconds
%
% tor wager

% Create onsets
% ----------------------------------------------------------------
conditions = 1:length(freqConditions);

nconditions = length(conditions);
ons = cell(1, nconditions);

% num_stim_slots = scanLength ./ ISI;
possible_onset_times = [0:ISI:scanLength-1]';

num_stim_slots = length(possible_onset_times);
    
for i = 1:nconditions

    % this is updated as we cycle through conditions...
    n = length(possible_onset_times);
    
    % Num of events for this condition
    nevents = round(freqConditions(i) .* num_stim_slots);
    
    rnd = randperm(n);
    
    ons{i} = sort(possible_onset_times(rnd(1:nevents)));
    
    ons{i}(:, 2) = eventduration;
    
    % remove these onset times from list of available onset times
    possible_onset_times(rnd(1:nevents)) = [];
    
end

end