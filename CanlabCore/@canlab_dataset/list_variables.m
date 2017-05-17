function [subj_descriptions, event_descriptions] = list_variables(D, varargin)
% Prints a list of variable names and descriptions for every variable, or specified variables
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

subj_descriptions = table;
event_descriptions = table;

subj_varnames = D.Subj_Level.names;
if iscolumn(subj_varnames), subj_varnames = subj_varnames'; end % need row

svars = find(strcmp('subj', varargin));
if ~isempty(svars), subj_varnames = varargin{svars+1}; end

% Names and descriptions - Subject level
% --------------------------------------------------------------------
Var_name = subj_varnames';

clear Description

for i=1:length(subj_varnames)
    
    vname = subj_varnames{i};
    [var,~,~,Description{i, 1}] = get_var(D, vname);
    
end

if length(subj_varnames) > 0
    subj_descriptions = table(Var_name, Description);
end

fprintf('\n\n --------- SUBJECT LEVEL VARS: data_obj.Subj_Level -------- \n\n');
disp(subj_descriptions);


% Names and descriptions - Event level
% --------------------------------------------------------------------
event_varnames = D.Event_Level.names;
if iscolumn(event_varnames), event_varnames = event_varnames'; end % need row

evars = find(strcmp('event', varargin));
if ~isempty(evars), event_varnames = varargin{evars+1}; end

Var_name = event_varnames';

clear Description

for i=1:length(event_varnames)
    
    vname = event_varnames{i};
    [var,~,~,Description{i, 1}] = get_var(D, vname);
    
end

if length(event_varnames) > 0
    event_descriptions = table(Var_name, Description);
end

fprintf('\n\n --------- EVENT LEVEL VARS: data_obj.Event_Level -------- \n\n');
disp(event_descriptions);
