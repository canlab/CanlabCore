function [names, ids, dat, descrips] = concatenate(D, varargin)
%
% Concatenates Subject-level and Event-level data across all subjects
%
% [names ids dat] = concatenate(D)
% 
%  INPUT:
%     D:  canlab_dataset
%     varargin:  currently accepts a wh_keep (logical array)
%
%  OUTPUT:
% Names: cell array of variable names
%
% Descrip: cell array of variable descriptions
%
% ids:   subject IDs matching data rows in dat
%
% dat:   subjects*events x variables matrix
%
%           - subject number, event number are included
%           - all subject-level and event-level data are included
%           - this format appropriate for, e.g., SAS/HLM
%
% Examples:
% [names, ids, flatdat] = concatenate(D);
% id_numbers = flatdat(:, 1);
%
% wh_subjs = true(size(D.Subj_Level.id));
% wh_subjs([13 18 19]) = false;
% [names, ids, dat] = concatenate(D, wh_subjs);
%
% % Copyright Tor Wager, 2013

% select subjects
wh_ids = true(length(D.Subj_Level.id), 1);
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