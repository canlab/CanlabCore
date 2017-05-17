function [names, ids, dat, descrips] = concatenate(D, varargin)
% Concatenates Subject-level and Event-level data across all subjects
%
% :Usage:
% ::
%
%    [names ids dat] = concatenate(D, [optional inputs])
%
% :Inputs:
%
%   **D:**
%        canlab_dataset object
%
% :Optional Inputs:
%   **a logical array**
%       a vector of 1/0 values to use as wh_keep
%
%
% :Outputs:
%
%   **names:**
%        cell array of variable names
%
%   **ids:**
%        subject IDs matching data rows in dat
%
%   **dat:**
%        subjects*events x variables matrix
%           - subject number, event number are included
%           - all subject-level and event-level data are included
%           - this format appropriate for, e.g., SAS/HLM
%
%   **descrips:**
%        cell array of variable descriptions
%
% :Examples:
% ::
%
%    [names, ids, flatdat] = concatenate(D);
%    id_numbers = flatdat(:, 1);
%
%    wh_subjs = true(size(D.Subj_Level.id));
%    wh_subjs([13 18 19]) = false;
%    [names, ids, dat] = concatenate(D, wh_subjs);
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

wh_ids = true(length(D.Subj_Level.id), 1); % select subjects
if length(varargin) > 0, wh_ids = logical(varargin{1}); end
whfind = find(wh_ids);

names = {'subj_num' D.Subj_Level.names{:} 'Event_number' D.Event_Level.names{:}};
descrips = {'subj_num' D.Subj_Level.descrip{:} 'Event_number' D.Event_Level.descrip{:}};


n = length(D.Subj_Level.id(wh_ids));

x = D.Event_Level.data(wh_ids);


for i = 1:n

    e = size(D.Event_Level.data{whfind(i)}, 1);
    
    ids{i} = repmat(D.Subj_Level.id{whfind(i)}, e, 1);
    
    %           subj ID                     subj level data              trialID  event data
    x{i} = [repmat(i, e, 1) repmat(D.Subj_Level.data(whfind(i), :), e, 1) (1:e)' x{i}];
    
end

ids = strvcat(ids{:});

dat = cat(1, x{:});

end % function
