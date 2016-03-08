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
  
    fprintf('%d subjects, %d subject-level vars, %d event-level vars\n', ...
        length(D.Subj_Level.id), length(D.Subj_Level.names), length(D.Event_Level.names));
    
    fprintf('\n\n --------- SUBJECT LEVEL VARS -------- \n\n');
    
    subj_varnames = D.Subj_Level.names;
    svars = find(strcmp('subj', varargin));
    if ~isempty(svars), subj_varnames = varargin{svars+1}; end

    
    for i=1:length(subj_varnames)
        
        vname = subj_varnames{i};
        [var,~,~,descrip] = get_var(D, vname);
        
        if iscell(var)
            % istext
            fprintf('%s (%s): Text. Unique values: %d\t\n', ...
            vname, descrip, length(unique(var)));
        
        else
            % isnumeric
        fprintf('%s (%s): min:%3.2f\t max:%3.2f\t mean:%3.2f\t sd:%3.2f NaNs:%d\n', ...
            vname, descrip, min(var), max(var), nanmean(var), nanstd(var), sum(isnan(var)));
        end
        
    end
    
    
    fprintf('\n\n --------- EVENT LEVEL VARS -------- \n\n');
   
    event_varnames = D.Event_Level.names;
    evars = find(strcmp('event', varargin));
    if ~isempty(evars), event_varnames = varargin{evars+1}; end
    
    for i=1:length(event_varnames)
        vname = event_varnames{i};
        [var,~,~,descrip] = get_var(D, vname);
        fprintf('%s (%s): min:%3.2f\t max:%3.2f\t mean:%3.2f NaNs:%d\n', ...
            vname, descrip, min(min(var)), max(max(var)), nanmean(nanmean(var)), sum(isnan(isnan(var))));
    end
end
