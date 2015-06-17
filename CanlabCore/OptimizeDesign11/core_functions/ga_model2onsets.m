function onsets = ga_model2onsets(M)
% This function takes a GA output M structure and returns onsets for each condition/event type in seconds
%
% onsets = ga_model2onsets(M)
%
% Zero is the start of the first event (i.e., the start of the run, after
% discarded acquistions)
%
% Onsets are stored in a cell array, one cell per condition, in order of
% the integer condition numbers stored in M.stimlist
%
% The onsets cell array is compatible with the format SPM5 (and 8?) uses
% for its "multiple conditions" input in the first-level experimental design
% specification
%
% You will have to add "names" and "durations" variables, but then you
% should be set to import the results of the GA into SPM
%
% Example:
% load model1_8-11-2009.mat
% onsets = ga_model2onsets(M)


    onsets = {};
    conditions = unique(M.stimlist(M.stimlist > 0));

    n = length(conditions);
    
    for i = 1:n
        
        ons = find(M.stimlist == conditions(i)) - 1; % onsets in TRs, zero is start of first onset
        
        ons = ons .* M.ga.ISI; % onsets in seconds
        
        onsets{i} = ons;
        
    end
    
    
        
end

