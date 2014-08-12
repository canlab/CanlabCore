% D = canlab_dataset();
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

classdef canlab_dataset
    
    properties
        
        % Build a dataset containing the key variables
        
        
        Description = struct('Experiment_Name', [], 'Missing_Values', [], 'Subj_Level', [], 'Event_Level', []);
        
        Subj_Level = struct('id', cell(1), 'names', cell(1), 'type', cell(1), 'units', cell(1), 'descrip', cell(1), 'data', []);
        
        Event_Level = struct('names', cell(1), 'type', cell(1), 'units', cell(1), 'descrip', cell(1), 'data', cell(1));
        
        Sub_Event_Level = struct('names', cell(1), 'type', cell(1), 'units', cell(1), 'descrip', cell(1), 'data', cell(1));
        
        Continuous = struct('names', cell(1), 'type', cell(1), 'units', cell(1), 'descrip', cell(1), 'data', cell(1));
        
        wh_keep = struct();
        
    end % properties
    
    %% Subject-level data
    
    methods
        
        % Class constructor
        function obj = canlab_dataset(varargin)
            
            obj.Description.Subj_Level{1, 1} = 'id: Unique text identifier for each subject (e.g., subject code)';
            
            obj.Description.Subj_Level{2, 1} = 'data: cell array, one cell per subject, with rectangular data matrix.';
            obj.Description.Subj_Level{3, 1} = '      Each row is one event, each column a variable, with names in .names';
            
            obj.Description.Event_Level{1, 1} = 'data: cell array, one cell per subject, with rectangular data matrix.';
            obj.Description.Event_Level{2, 1} = '      Each row is one event, each column a variable, with names in .names';
            
            
        end % constructor function
        
    end % methods
    
    
end % classdef



