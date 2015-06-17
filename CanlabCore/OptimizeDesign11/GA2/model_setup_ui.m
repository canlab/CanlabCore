% model_setup_ui
% Run this to specify a design for optimization with GA2
%
% The script creates two struct variables, mspec and conditions, 
% which are then used as inputs to ga2.m
%
% mspec and conditions contain fields, accessed with mspec.fieldname,
% that contain information about the design.  
% Numbers in these variables describe timing info in TRs! (except the runtime and TR fields, which are in s)
%
% This function is good for:
%   1 - creating a design matrix with combinations of convolved, unconvolved, and shape-free elements
%   2 - creating designs AT THE RESOLUTION OF THE TR
%   3 - creating designs for optimization of temporal jitter and estimation of parts of a single trial
%   4 - creating designs with completely independent conditions (can overlap, etc.)
%
% This function is NOT as good for:
%   1 - optimizing the ordering of stimulus presentation for multiple trial types
%   2 - creating designs that have 1 or 2 events in some conditions.  confuses the program.
%   3 - creating designs in which trial types are mutually exclusive (present A or B at time t, but not both)

% scanlength, tr

% event, epoch, block
% convolve
% regular / random

% each condition is a trial type
% each trial type can be
%   event
%   epoch
%   multi-part trial

% each condition has these fields:
%   onsets      vector, for fixed onsets, times in s
%               [onset range; trial length range] for random onsets, in s
%               0 is first second of run
%
%   parts       epochs within a trial, for multipart trials
%               for each part, specify onsets (vector or range)
%               onsets can be an integer, for fixed onsets rel to trial or prev part.

% for multi-part
%   no of parts
%   

clear mspec, clear conditions

conditions(1).onsets = 1;
i = 1;

disp(['Model setup user-interface: for specifying a design matrix.'])
disp(['by Tor Wager, version a1'])
disp(['___________________________________________________________'])

mspec.TR = input('Enter TR in s ');
mspec.runtime = input('Enter run length in s ');
mspec.hrf = spm_hrf(mspec.TR); mspec.hrf = mspec.hrf ./ max(mspec.hrf);
mspec.numframes = mspec.runtime ./ mspec.TR;

tmp = input('All trials at least N TRs apart? (Enter N or return to skip.) ');
if ~isempty(tmp), mspec.toverlap = [1 tmp];,end

disp(['Condition setup:'])
disp(['Press return when asked for onsets for a condition to exit.'])
disp(['___________________________________________________________'])
  


while ~isempty(conditions(i).onsets)
    
    conditions(i).onsets = input(['Onsets for condition ' num2str(i) ' in s ']);
    conditions(i).onsets = conditions(i).onsets ./ mspec.TR;
    
    if isempty(conditions(i).onsets), break, end
    
    conditions(i).parts = input(['Parts for condition ' num2str(i) ': enter 1 or num parts ']);
    % if this is a multi-part trial, parts > 1
        
        for j = 1:conditions(i).parts
            
            rel_to = input(['Onsets for condition ' num2str(i) ' part ' num2str(j) ' relative to: 1 for start of trial, 2 for prev. part ']);
            if rel_to == 1, str = ' start of trial ';, elseif rel_to == 2, str = ' previous part ';, else, error('Choose 1 or 2'), end
            
            if rel_to == 2 & j == 1, error('Cannot specify onset of 1st trial part relative to previous part.'),end
            conditions(i).subcond(j).rel_to = rel_to;
            conditions(i).subcond(j).onsets = input(['Onsets in s for condition ' num2str(i) ' part ' num2str(j) ' from' str]);
            if length(conditions(i).subcond(j).onsets) == 1, 
                conditions(i).subcond(j).onsets = repmat(conditions(i).subcond(j).onsets,1,length(conditions(i).onsets));
            end
            conditions(i).subcond(j).onsets = conditions(i).subcond(j).onsets ./ mspec.TR;
            
            conditions(i).subcond(j).stimlength = input(['Stimulation length in TRs for condition ' num2str(i) ' part ' num2str(j) ' ']);
       
            conditions(i).subcond(j).convolve = input(['Convolve condition ' num2str(i) ' part ' num2str(j) ' with canonical HRF? 1/0 ']);
            if ~conditions(i).subcond(j).convolve
                conditions(i).subcond(j).hrfest = input(['Estimate shape-free HRF for condition ' num2str(i) ' part ' num2str(j) ' ? 1 or TP to estimate ']);
            else
                conditions(i).subcond(j).hrfest = [];
            end
            
        end
            
    if conditions(i).parts == 1, conditions(i).stimlength = conditions(i).subcond(1).stimlength;, end
        
    i = i + 1;
    conditions(i).onsets = 1;
end

conditions = conditions(1:end-1);