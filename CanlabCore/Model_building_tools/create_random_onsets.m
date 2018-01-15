function ons = create_random_onsets(scanLength, ISI, freqConditions, eventduration)
% ons = create_random_onsets(scanLength, ISI, freqConditions, eventduration)
% tor wager

% Create onsets
% ----------------------------------------------------------------
conditions = 1:length(freqConditions);

nconditions = length(conditions);
ons = cell(1, nconditions);


num_stim_slots = scanLength ./ ISI;
possible_onset_times = [0:ISI:(num_stim_slots - 1)]';


for i = 1:nconditions
    
    n = length(possible_onset_times);
    
    nevents = round(freqConditions(i) .* num_stim_slots);
    
    rnd = randperm(n);
    
    ons{i} = sort(possible_onset_times(rnd(1:nevents)));
    
    ons{i}(:, 2) = eventduration;
    
    % remove these onset times from list of available onset times
    possible_onset_times(rnd(1:nevents)) = [];
    
end

end