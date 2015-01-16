function [headername, dataname, fid] = write_text(D, varargin)
% function [headername, dataname, fid] = write_text(D)
%
% "Flatten" dataset and write text files with header and data
% For all Event-level and Subject-level data.  Files are created in the
% current working directory.
%
% first varargin parameter is the delimiter.  Comma-delimited by default
%
% % Copyright Tor Wager, 2013


% flesh this out later as needed
delim = ',';
for i=1:length(varargin)
    delim = varargin{1};
end


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