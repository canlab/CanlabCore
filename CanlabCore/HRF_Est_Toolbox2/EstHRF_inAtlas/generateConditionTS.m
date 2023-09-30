function Condition=generateConditionTS(fmri_d, conditions, onsets, durations)
    % Helper script to create condition vectors for hrf_fit_one_voxel()
    % Michael Sun, Ph.D.
    % - Takes fmri_data() object or the number of TRs in a 4D object.
    % - conditions: cell-vector of cellstrings for each condition e.g., {'hot', 'warm', 'imagine'}
    % - onsets: SPM-style onsets cell array in seconds from first-level model, e.g., {{1,2,3}, {4, 5, 6}, {7, 8, 9}}
    % - duration: SPM-style duration cell array in seconds from first-level model, e.g., {{12, 12 ,12}, {12, 12, 12}, {12, 12, 12}} 
    %
    % *Usage:
    % ::
    %    Condition = generateConditionTS(image_obj, {'hot','warm','imagine'}, onsets, durations})
    %
    % Note 1: Preset a SPIKES or SPIKETRAINS variable in order to toggle the
    % generation of Spikes (single 1s) or Spiketrains (a train of 1s) to
    % represent each event.
    %
    % Note 2: If there was slice-timing correction performed, then you will
    % want to correct for it by adding 0.5TRs to all of your times.

    % Default to modeling single spike instead of duration of events.
    if ~exist('SPIKES', 'var') &&  ~exist('SPIKETRAINS', 'var')
        SPIKES=1;
        SPIKETRAINS=0;
    end
    
    if strcmp(class(fmri_d), 'fmri_data')
        n=size(fmri_d.dat,2);
    elseif isnumeric(fmri_d)
        n=fmri_d;
    end


    % Number of conditions
    nconds = length(conditions);
    
    % Initialize the Condition cell array
    for c = 1:nconds
        Condition{c} = zeros(n,1);
    end
    
    % Loop through each condition to process their times and update Condition
    for c = 1:nconds
        % Extract the times for the current condition
        % current_times = eval([conditions{c} '_times{sub}{j}']);
        current_times=onsets{c};
        current_dur=durations{c};
        
        % Remove the onset times that didn't get recorded. 
        % Note: Add 0.5TRs to correct for slice-timing correction if needed
        % current_times = ceil(current_times(current_times < n) + 0.5);
        current_times = ceil(current_times(current_times < n));
        
        % Model events as a single impulse
        if SPIKES == 1
            Condition{c}(current_times) = 1;
        end
        
        % Model event epochs as trains of events
        if SPIKETRAINS == 1
            for i = 1:numel(current_times)    
                Condition{c}(current_times(i):current_times(i)+(current_dur(i))) = 1;
            end
        end
    end
end