function print_summary(D, varargin)
% Prints summaries for every variable, or specified variables
%
% :Usage:
% ::
%
%    print_summary(D, [optional inputs]) 
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2013 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
%
% :Inputs:
%
%   **D:**
%        a canlab_dataset object
%
% :Optional Inputs:
%
%   **subj:**
%        followed by a cell array of subject level var names, to only see those vars
%
%   **event:**
%        followed by a cell array of event level var names, to only see those vars
%
%   if either varargin is unspecified, all variables will be printed
%

    fprintf('\n\n --------- DATASET VARS -------- \n\n');
  
    fprintf('Experiment name:\n%s\n\n', D.Description.Experiment_Name);
    
    fprintf('%d subjects, %d subject-level vars, %d event-level vars\n', ...
        length(D.Subj_Level.id), length(D.Subj_Level.names), length(D.Event_Level.names));
    
    % Names and descriptions
    % --------------------------------------------------------------------
    [subj_descriptions, event_descriptions] = list_variables(D, varargin{:});
    
    % Descriptive stats for each variable
    % --------------------------------------------------------------------  
    [subj_descriptives, event_descriptives] = get_descriptives(D);
    
    
    fprintf('\n\n --------- SUBJECT LEVEL VARS: data_obj.Subj_Level -------- \n\n');
    
    disp(subj_descriptives)

    
    fprintf('\n\n --------- EVENT LEVEL VARS: data_obj.Event_Level -------- \n\n');
   
    disp(event_descriptives)
    
end % main function
