% D = canlab_dataset(varargin);
%
% Constructs a new, empty instance of a canlab_dataset object.
% The fields of this object can subsequently be assigned data.
%
% The dataset (D) The dataset (D) is a way of collecting behavioral data 
% and/or meta-data about fMRI (and other) studies in a standardized format.
% It has entries for experiment description, Subject-level
% variables, Event-level variables (e.g., for fMRI events or trials in a
% behavioral experiment) and other fields for sub-event level  data (e.g., 
% physio or continuous ratings within trials).
% 
% See D.Description for some basic info.
%
%
% varargin options: 'fmri'  -  Adds standard fmri variable names to the
%                               dataset under Event_Level
%
% Dataset methods include:
%
% write_text    -> writes a flat text file (data across all subjects)
% concatenate   -> returns flat matrix of data across all subjects to the workspace
% get_var       -> get event-level data from one variable and return in matrix and cell formats
% scatterplot   -> standard or multi-line scatterplot of relationship between two variables
% scattermatrix -> Scatterplot matrix of pairwise event-level variables
%
% - type methods(canlab_dataset) for a list of all methods
%
% Copyright Tor Wager, 2013
%
% Conventions for creating objects:
% --------------------------------------------------------------------------------------
% Follow these so your dataset object will work with functions!!
% If you have n subjects and k variables, 
% 1. obj.Subj_Level.id is an n-length cell vector of strings with subject IDs
% 2. obj.Subj_Level.type is a k-length cell vector with 'text' or 'numeric' for each variable
% 3. obj.Subj_Level.data is n x k numeric data, with columns of NaNs if the
%     variable is a text variable
% 4. obj.Subj_Level.textdata is an n x k cell array of text data, with columns of empty cells if the
%     variable is a numeric variable
% 5. If you have fMRI onsets and durations, etc., use the 'fmri' option to
%    create standardized variable names.

classdef canlab_dataset
    
    properties
        
        % Build a dataset containing the key variables
        
        
        Description = struct('Experiment_Name', [], 'Missing_Values', [], 'Subj_Level', [], 'Event_Level', []);
        
        Subj_Level = struct('id', cell(1), 'names', cell(1), 'type', cell(1), 'units', cell(1), 'descrip', cell(1), 'data', cell(1), 'textdata', cell(1));
        
        Event_Level = struct('names', cell(1), 'type', cell(1), 'units', cell(1), 'descrip', cell(1), 'data', cell(1), 'textdata', cell(1));
        
        Sub_Event_Level = struct('names', cell(1), 'type', cell(1), 'units', cell(1), 'descrip', cell(1), 'data', cell(1), 'textdata', cell(1));
        
        Continuous = struct('names', cell(1), 'type', cell(1), 'units', cell(1), 'descrip', cell(1), 'data', cell(1), 'textdata', cell(1));
        
        wh_keep = struct();
        
    end % properties
    
    %% Subject-level data
    
    methods
        
        % Class constructor
        function obj = canlab_dataset(varargin)
            
            % Subject Level Descriptions
            obj.Description.Subj_Level{1, 1} = 'id: Unique text identifier for each subject (e.g., subject code)';
            
            obj.Description.Subj_Level{2, 1} = 'names: cell array containing names for each variable in .data and .textdata.';
            
            obj.Description.Subj_Level{3, 1} = 'type: cell array identifying variable type; used to index .data or .textdata.';
            obj.Description.Subj_Level{4, 1} = '      Acceptable values are either ''text'' or ''numeric''';
            
            obj.Description.Subj_Level{5, 1} = 'units: cell array identifying units used for each variable (e.g. seconds, count, text)';
            
            obj.Description.Subj_Level{6, 1} = 'descrip: cell array providing detailed description of each variable in .names';
            
            obj.Description.Subj_Level{7, 1} = 'data: cell array, one cell per subject, with rectangular data matrix.';
            obj.Description.Subj_Level{8, 1} = '      Each row is one subject, each column a variable, with names in .names';
            obj.Description.Subj_Level{9, 1} = '      Numeric data only, text data is assigned as NaN';
            
            obj.Description.Subj_Level{10, 1} = 'textdata: cell array, one cell per subject, with cell matrix.';
            obj.Description.Subj_Level{11, 1} = '      Each row is one subject, each column a variable, with names in .names';
            obj.Description.Subj_Level{12, 1} = '      Text data only, numeric data is assigned as ''''';
            
            %Event Level Descriptions
            obj.Description.Event_Level{1, 1} = 'names: cell array containing names for each variable in .data and .textdata.';
            
            obj.Description.Event_Level{2, 1} = 'type: cell array identifying variable type; used to index .data or .textdata.';
            obj.Description.Event_Level{3, 1} = '      Acceptable values are either ''text'' or ''numeric''';
            
            obj.Description.Event_Level{4, 1} = 'units: cell array identifying units used for each variable (e.g. seconds, count, text)';
            
            obj.Description.Event_Level{5, 1} = 'descrip: cell array providing detailed description of each variable in .names';
            
            obj.Description.Event_Level{6, 1} = 'data: cell array, one cell per subject, with rectangular data matrix.';
            obj.Description.Event_Level{7, 1} = '      Each row is one event, each column a variable, with names in .names';
            obj.Description.Event_Level{8, 1} = '      Numeric .type data only, text .type data is assigned as NaN';
            
            obj.Description.Event_Level{9, 1} = 'textdata: cell array, one cell per subject, with cell matrix.';
            obj.Description.Event_Level{10, 1} = '      Each row is one event, each column a variable, with names in .names';
            obj.Description.Event_Level{11, 1} = '      Text .type data only, numeric .type data is assigned as ''''';
            
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch(varargin{i})
                        case 'fmri'
                            obj.Event_Level.names = {'SessionNumber' 'RunName' 'RunNumber' 'TaskName' 'TrialNumber' 'EventName' 'EventOnsetTime' 'EventDuration'};
                            obj.Event_Level.type = {'numeric' 'text' 'numeric' 'text' 'numeric' 'text' 'numeric' 'numeric'};
                            obj.Event_Level.units = {'count' 'text' 'count' 'text' 'count' 'text' 'seconds' 'seconds'};
                            obj.Event_Level.descrip = {'Within-study session number (e.g. 1,2,3)'
                                                       'A descriptive name of the fMRI run (e.g. PreRevealPlacebo)'
                                                       'Within-session run number (e.g. 1,2,3)'
                                                       'The name of the task being performed during this event (e.g. Compassion)'
                                                       'The trial number within the task for this event (e.g. 1,2,3)'
                                                       'The name of this event (e.g. HighPain, Cue, Fixation)'
                                                       'The run time when this event starts'
                                                       'The length of time this event lasts'}';
                            varargin{i} = [];
                        otherwise, warning(['Unknown input string operation:' varargin{i}])
                    end
                end
            end
            
            
        end % constructor function
        
    end % methods
    
    
end % classdef



