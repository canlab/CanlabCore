function [subj_descriptives, event_descriptives] = get_descriptives(D, varargin)
% Returns descriptive stats for every variable, or specified variables
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
% Examples:
%
% % Get descriptive stats in tables all variables in canalab_dataset object D
% [subj_descriptives, event_descriptives] = get_descriptives(D, 'event', 'temp')
%
% Get descriptive stats in tables for the Event_Level variable 'temp'
% [subj_descriptives, event_descriptives] = get_descriptives(D, 'event', 'temp')
%
% Get descriptive stats in tables for the Event_Level variables 'temp' and 'ratings'
%[subj_descriptives, event_descriptives] = get_descriptives(D, 'event', {'temp' 'ratings'})

subj_descriptives = table;
event_descriptives = table;

% Default: All variables
subj_varnames = D.Subj_Level.names;
if iscolumn(subj_varnames), subj_varnames = subj_varnames'; end % need row

event_varnames = D.Event_Level.names;
if iscolumn(event_varnames), event_varnames = event_varnames'; end % need row

svars = find(strcmp('subj', varargin));
evars = find(strcmp('event', varargin));

if ~isempty(svars) || ~isempty(evars)
    % Replace default lists with input variables only
    [subj_varnames, event_varnames] = deal({}); 
end

    
% Subject level
% -------------------------------------------------------------------------

if ~isempty(svars), subj_varnames = varargin{svars+1}; end
if ~isempty(svars) && ~iscell(subj_varnames), subj_varnames = {subj_varnames}; end

Var_name = subj_varnames';
clear Min_value Max_value Mean_value St_Dev IQR NaN_count

for i=1:length(subj_varnames)
    
    vname = subj_varnames{i};
    
    [Min_value(i, 1) Max_value(i, 1) Mean_value(i, 1) St_Dev(i, 1) IQR(i, 1) NaN_count(i, 1), text_vals] = get_descriptives_var(D, vname);
    
end

if length(subj_varnames) > 0
    subj_descriptives = table(Var_name, Min_value, Max_value, Mean_value, St_Dev, IQR, NaN_count);
end

% Event level
% -------------------------------------------------------------------------
if ~isempty(evars), event_varnames = varargin{evars+1}; end
if ~isempty(evars) && ~iscell(event_varnames), event_varnames = {event_varnames}; end

Var_name = event_varnames';
clear Min_value Max_value Mean_value St_Dev IQR NaN_count

for i=1:length(event_varnames)
    
    vname = event_varnames{i};
    
    [Min_value(i, 1) Max_value(i, 1) Mean_value(i, 1) St_Dev(i, 1) IQR(i, 1) NaN_count(i, 1), text_vals] = get_descriptives_var(D, vname);
    
end

if length(event_varnames) > 0
    event_descriptives = table(Var_name, Min_value, Max_value, Mean_value, St_Dev, IQR, NaN_count);
end

end % main function



function [Min_value, Max_value, Mean_value, St_Dev, myIQR, NaN_count, text_vals] = get_descriptives_var(D, vname)

% intialize
[Min_value, Max_value, Mean_value, St_Dev, myIQR, NaN_count] = deal(NaN);
text_vals = [];

[var, varcell, ~, descrip] = get_var(D, vname);

if isempty(var)
    return
end

if ismatrix(var)
    % Happens when this is Event_Level var with same number of events for
    % all subjects; or padded with NaNs.
    var = var(:);
end

if ischar(var)
    % This happens for event-level variables with different numbers of
    % observations for different cells. This issue may be removed in the
    % future, however.  
    
    if iscell(varcell{1})
        fprintf('%s (%s): Cell array. Length of first cell: %d\t\n', ...
        vname, descrip, length(varcell{1}));
    
    else
    
        var = double(cat(1, varcell{:}));
        
    end
end

if iscell(var)  % SHOULD BE DEPRECATED - need to handle text differently
    % istext
    fprintf('%s (%s): Text. Unique values: %d\t\n', ...
        vname, descrip, length(unique(var)));
     
else 
    % Standard Subj-level and Event-level vars with same no. of obs.
    % isnumeric
    Min_value = min(var);
    Max_value = max(var);
    Mean_value = nanmean(var);
    St_Dev = nanstd(var);
    myIQR = iqr(var);
    NaN_count = sum(isnan(var));
end

end


