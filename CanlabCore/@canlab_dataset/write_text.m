function [headername, dataname, fid] = write_text(D, varargin)
% "Flatten" dataset and write text files with header and data
% For all Event-level and Subject-level data.  Files are created in the
% current working directory.
%
% :Usage:
% ::
%
%    function [headername, dataname, fid] = write_text(D, [optional inputs])
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
%    the first varargin parameter is the delimiter. Comma-delimited by default.
%    the second varargin parameter is a logical vector of subjects to use. All subjects by default.
%
%
% :Outputs:
%
%   **headername:**
%        filename of header output file
%
%   **dataname:**
%        filename of data output file
%
%   **fid:**
%        file ID (currently does not look to be used)
%


delim = ','; % flesh this out later as needed
wh_keep = true(size(D.Subj_Level.id));

if nargin > 1, delim = varargin{1}; end
if nargin > 2, wh_keep = varargin{2}; end



% Checks
% ----------------------------------------------------------------------
for i = 1:length(D.Subj_Level.names)
    
    if length(D.Subj_Level.descrip) < i
        D.Subj_Level.descrip{i} = 'No description provided';
    end
    
end

for i = 1:length(D.Event_Level.names)
    
    if length(D.Event_Level.descrip) < i
        D.Event_Level.descrip{i} = 'No description provided';
    end
    
end

if isempty(D.Event_Level.data)
    D.Event_Level.data = cell(1,length(D.Subj_Level.id));
end

if isempty(D.Event_Level.names)
    D.Event_Level.names = {};
end



% Open Files
% ----------------------------------------------------------------------
headername = [D.Description.Experiment_Name '_info_' scn_get_datetime '.txt'];
dataname = [D.Description.Experiment_Name '_data_' scn_get_datetime '.csv'];
fid = fopen(headername, 'w');


% Write Header
% ----------------------------------------------------------------------

u = '______________________________________________________________';

fprintf(fid, 'Experiment: %s\n', D.Description.Experiment_Name);

fprintf(fid, '\n%d subjects\n', length(D.Subj_Level.id));
fprintf(fid, 'Missing values coded with: %f\n', D.Description.Missing_Values);
fprintf(fid, '%s\n', u);

fprintf(fid, 'Subject Level\n\n');
fprintf(fid, 'Description:\n');

for i = 1:length(D.Description.Subj_Level)
    fprintf(fid, '\t%s\n', D.Description.Subj_Level{i});
end

fprintf(fid, 'Names:\n');
for i = 1:length(D.Subj_Level.names)
    fprintf(fid, '\t%s\t%s\n', D.Subj_Level.names{i}, D.Subj_Level.descrip{i});
end

fprintf(fid, '%s\n', u);

fprintf(fid, 'Event Level\n\n');
fprintf(fid, 'Description:\n');

for i = 1:length(D.Description.Event_Level)
    fprintf(fid, '\t%s\n', D.Description.Event_Level{i});
end

fprintf(fid, 'Names:\n');
for i = 1:length(D.Event_Level.names)
    fprintf(fid, '\t%s\t%s\n', D.Event_Level.names{i}, D.Event_Level.descrip{i});
end

fprintf(fid, '%s\n', u);

fclose(fid);

%% Write subj-level & event-level data
% -----------------------------------------------------------------------

fid = fopen(dataname, 'w');

if isempty(D.Event_Level.names)
    names = {'id' D.Subj_Level.names{:}};
else
    names = {'id' D.Subj_Level.names{:} 'Event_number' D.Event_Level.names{:}};
end

n = length(D.Subj_Level.id);

%slevels = length(D.Subj_Level.names);

printcell = @(x) fprintf(fid, ['%s' delim], x);
cellfun(printcell, names);
fprintf(fid, '\n');

for i = 1:n  % for each subject
    
    if ~wh_keep(i), continue; end % skip this subj

    e = size(D.Event_Level.data{i}, 1);
    
    if e==0 % no events
        fprintf(fid, ['%s' delim], D.Subj_Level.id{i});  % ID, can be char
        datarow = [D.Subj_Level.data(i, :)];
        printrow(fid, datarow, delim);        
    else     
        for j = 1:e  % for events within subject
            fprintf(fid, ['%s' delim], D.Subj_Level.id{i});  % ID, can be char
            datarow = [D.Subj_Level.data(i, :) j D.Event_Level.data{i}(j, :)];
            printrow(fid, datarow, delim);

        end  % events
    end
end  % subjects

fclose(fid);


end % function

function printrow(fid, datarow, delim)
    % Could switch by datatype in a better way here!!  need datatype codes in dataset.
    for k = 1:length(datarow)
        if datarow(k) == round(datarow(k))
            fprintf(fid, ['%d' delim], datarow(k));
        else
            fprintf(fid, ['%3.3f' delim], datarow(k));  
        end
    end % row

    fprintf(fid, '\n');
end
