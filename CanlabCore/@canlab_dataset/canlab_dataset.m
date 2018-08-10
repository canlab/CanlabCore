% D = canlab_dataset(varargin);
%
% Constructs a new, empty instance of a canlab_dataset object.
% The fields of this object can subsequently be assigned data.
%
% The dataset (D) is a way of collecting behavioral data 
% and/or meta-data about fMRI (and other) studies in a standardized format.
% It has entries for experiment description, Subject-level
% variables, Event-level variables (e.g., for fMRI events or trials in a
% behavioral experiment) and other fields for sub-event level  data (e.g., 
% physio or continuous ratings within trials).
% 
%
% varargin options: 'fmri'  -  Adds standard fmri variable names to the
%                               dataset under Event_Level
%
% The canlab_dataset object contains data at five potential levels of analysis:
% 1. Experiment-level
% 2. Subject_level
% 3. Event_Level
% 4. Sub_Event_Level
% 5. Continuous
%
% 1. Experiment-level describes entire dataset 
%
% 2. Subject-level data contains one value per person per variable.
%    It is entered in DAT.Subj_Level.data and DAT.Subj_Level.textdata
%    In a matrix of N subjects x V variables
%
% 3. Event_Level data contains one value per event/trial per subject per variable.
%    It is entered in DAT.Event_Level.data and DAT.Event_Level.textdata
%    In a cell array with one cell per subject.
%    Each cell contains a matrix of T trials/events x V variables
%
% 4. Sub_Event_Level data contains one or more values per event/trial per subject per variable.
%    It is entered in DAT.Sub_Event_Level.data and DAT.Sub_Event_Level.textdata
%    In a cell array with one cell per subject.
%    Each cell contains a variable-size matrix of K observations x V variables
%
% 5. Continuous data is intended to contain time-series of K samples per subject per variable.
%    It is entered in DAT.Continuous.data and DAT.Continuous.textdata
%    In a cell array with one cell per subject.
%    Each cell contains a variable-size matrix of K samples x V variables
%
% The most important levels for many analyses are the Subject_level and
% Event_level, which provide sufficient data to run multi-level mixed
% effects models.
%
% Documenting meta-data:
% -----------------------------------------------------------------------
% All data levels have fields for variable names, stored in cell arrays:
% DAT.Subj_Level.names
% DAT.Event_Level.names
% etc.
% All data levels also have fields for units of analysis and data type.
% Data type is either 'numeric' or 'text'
% All data levels have fields for long-form descriptions of variables,
% stored in .descrip, in cell arrays with one cell per variable.
% A complete dataset will make use of these important fields.
%
% For more information, see:
%   help(canlab_dataset)
%   disp(DAT.Description.Subj_Level)
%   disp(DAT.Description.Event_Level)
%
% Dataset methods include:
%
% write_text    -> writes a flat text file (data across all subjects)
% concatenate   -> returns flat matrix of data across all subjects to the workspace
% get_var       -> get event-level data from one variable and return in matrix and cell formats
% scatterplot   -> standard or multi-line scatterplot of relationship between two variables
% scattermatrix -> Scatterplot matrix of pairwise event-level variables
%
% More methods:
% add_var             get_var             mediation           scattermatrix       write_text          
% bars                glm                 plot_var            scatterplot         
% canlab_dataset      glm_multilevel      print_summary       spm2canlab_dataset  
% concatenate         histogram           read_from_excel     ttest2  
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
%
% fmri-type object
% --------------------------------------------------------------------------------------
% This has standard variable names, which enable the creation of first-level fMRI designs
% from canlab_dataset objects.
% These are the variables you need at minimum for an fmri-type object:
%
%         Var_name                                      Description                               
%     ________________    ________________________________________________________________________
% 
%     'SessionNumber'     'Within-study session number (e.g. 1,2,3)'                              
%     'RunName'           'A descriptive name of the fMRI run (e.g. PreRevealPlacebo)'            
%     'RunNumber'         'Within-session run number (e.g. 1,2,3)'                                
%     'TaskName'          'The name of the task being performed during this event (e.g. Compassi?'
%     'TrialNumber'       'The trial number within the task for this event (e.g. 1,2,3)'          
%     'EventName'         'The name of this event (e.g. HighPain, Cue, Fixation)'                 
%     'EventOnsetTime'    'The start time for each event in sec from start of the run'            
%     'EventDuration'     'The length of time this event lasts in sec'  
%
% Examples:
% --------------------------------------------------------------------------------------
% 
% See "canlab_dataset_basic_usage.html" in the Canlab_help_examples
% repository.  Additional code snippets are below.
%
% % Load sample data file (this is in Sample_datasets in Core - make sure it's on your path)
% load('Jepma_IE2_single_trial_canlab_dataset.mat'); % loads object called DAT 
%
% % Print a summary:      
% print_summary(DAT);     % see also list_variables(DAT)
%
% % Extract expectation and pain ratings, two continuous event-level variables from the object
% [~, expect] = get_var(DAT, 'expected pain');
% [~, pain] = get_var(DAT, 'reported pain');
%
% % Make a multi-level scatterplot 
% create_figure('lines_expect_pain', 1, 2);
% [han, Xbin, Ybin] = line_plot_multisubject(expect, pain, 'n_bins', 5, 'group_avg_ref_line', 'MarkerTypes', 'o', 'colors', custom_colors([1 .7 .4], [1 .7 .4], 100), 'gcolors', {[.4 .2 0] [.8 .4 0]});
% axis tight
% xlabel('Expectations'); ylabel('Pain'); set(gca, 'FontSize', 24);
%
% % Now extract expectation and pain variables conditional on cue value:
% [~, expect_cshigh] = get_var(DAT, 'expected pain', 'conditional', {'cue valence' 1});
% [~, expect_csmed] = get_var(DAT, 'expected pain', 'conditional', {'cue valence' 0});
% [~, expect_cslow] = get_var(DAT, 'expected pain', 'conditional', {'cue valence' -1});
% [~, pain_cshigh] = get_var(DAT, 'reported pain', 'conditional', {'cue valence' 1});
% [~, pain_csmed] = get_var(DAT, 'reported pain', 'conditional', {'cue valence' 0});
% [~, pain_cslow] = get_var(DAT, 'reported pain', 'conditional', {'cue valence' -1});
% 
% % Now make the plot, without individual subject lines
% % First set colors and condition names
% color1 = {[.9 .4 .2] [.6 .3 0]};
% color2 = {[.5 .5 .5] [.2 .2 .2]};
% color3 = {[.4 .8 .4] [.2 .7 .2]};
% condition_names = {'High Cue' 'Medium Cue' 'Low Cue'};
%
% % Now plot:
% subplot(1, 2, 2)
% [han1, Xbin, Ybin] = line_plot_multisubject(expect_cshigh, pain_cshigh, 'n_bins', 5, 'group_avg_ref_line', 'MarkerTypes', 'o', 'colors', color1, 'center', 'noind', 'nolines');
% [han2, Xbin, Ybin] = line_plot_multisubject(expect_csmed, pain_csmed, 'n_bins', 5, 'group_avg_ref_line', 'MarkerTypes', 'o', 'colors',color2, 'center', 'noind', 'nolines');
% [han3, Xbin, Ybin] = line_plot_multisubject(expect_cslow, pain_cslow, 'n_bins', 5, 'group_avg_ref_line', 'MarkerTypes', 'o', 'colors',color3, 'center', 'noind', 'nolines');
% xlabel('Expectations'); ylabel('Pain'); set(gca, 'FontSize', 24);
% legend([han1.grpline_handle(1) han2.grpline_handle(1) han3.grpline_handle(1)], condition_names);
% drawnow, snapnow

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
            
            % Enforce cells, which we could not assign before
            [obj.Subj_Level.id, obj.Subj_Level.names, obj.Subj_Level.type, obj.Subj_Level.units, obj.Subj_Level.descrip, obj.Subj_Level.textdata] = deal({});

            [obj.Event_Level.names, obj.Event_Level.type, obj.Event_Level.units, obj.Event_Level.descrip, obj.Event_Level.data, obj.Event_Level.textdata] = deal({});
            
            [obj.Sub_Event_Level.names, obj.Sub_Event_Level.type, obj.Sub_Event_Level.units, obj.Sub_Event_Level.descrip, obj.Sub_Event_Level.data, obj.Sub_Event_Level.textdata] = deal({});

            [obj.Continuous.names, obj.Continuous.type, obj.Continuous.units, obj.Continuous.descrip, obj.Continuous.data, obj.Continuous.textdata] = deal({});

            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch(varargin{i})
                        case 'fmri'
                            
                            obj.Subj_Level.names = {'Sex' 'Age' 'Race' 'Hispanic' 'Handedness'};
                            obj.Subj_Level.type = {'text' 'numeric' 'text' 'numeric' 'text' 'text'};
                            obj.Subj_Level.descrip = {'M or F, male/female'
                                                       'Participant age'
                                                       'Af_Am/White/Native_Am/Asian/PacificIslander/Other'
                                                       'Hispanic: 1 = yes or 0 = no'
                                                       'L or R, left or right-handed'}';
                            obj.Subj_Level.units = {'text' 'years' 'text' 'binary' 'text'};
                                                   
                                                   
                            obj.Event_Level.names = {'SessionNumber' 'RunName' 'RunNumber' 'TaskName' 'TrialNumber' 'EventName' 'EventOnsetTime' 'EventDuration'};
                            obj.Event_Level.type = {'numeric' 'text' 'numeric' 'text' 'numeric' 'text' 'numeric' 'numeric'};
                            obj.Event_Level.units = {'count' 'text' 'count' 'text' 'count' 'text' 'seconds' 'seconds'};
                            obj.Event_Level.descrip = {'Within-study session number (e.g. 1,2,3)'
                                                       'A descriptive name of the fMRI run (e.g. PreRevealPlacebo)'
                                                       'Within-session run number (e.g. 1,2,3)'
                                                       'The name of the task being performed during this event (e.g. Compassion)'
                                                       'The trial number within the task for this event (e.g. 1,2,3)'
                                                       'The name of this event (e.g. HighPain, Cue, Fixation)'
                                                       'The start time for each event in sec from start of the run'
                                                       'The length of time this event lasts in sec'}';
                            varargin{i} = [];
                        otherwise, warning(['Unknown input string operation:' varargin{i}])
                    end
                end
            end
            
        end % constructor function
        
    end % methods
    
    
end % classdef



