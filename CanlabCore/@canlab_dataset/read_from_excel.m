function dat = read_from_excel(dat, ExperimentFileName, SubjectFileList, varargin)
% Read from datafile into canlab_dataset format - currently requires file
% extensions .xls or .xlsx, but in the future will use importdata to take
% .csv or .txt extensions as well.
%
% Datafiles require column headers
%     - Design file requires: id, names, units, descrip
%               - other columns can be added
%               - ONLY between subject columns identified in the 'names'
%                   column are added.
%               - See Sample_canlab_dataset_experiment_level.xlsx for an
%                   example design file
%     - Subject files require no specific column headers, but all column
%               headers must be identical across all subjects.
%               - Enter NaN for data field in file if no value for that column within a specific Event
%               - ALL columns of subject files are written to canlab_dataset
% 
% :Usage:
% ::
%
%    dat = read_from_excel(dat, ExperimentFileName, SubjectFileList, [optional inputs])
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
%   **dat:**
%        a canlab_dataset object
%
%   **ExperimentFileName:**
%        the absolute path of the experiment data file
%
%   **SubjectFileList:**
%        list of absolute paths for individual subject files
%             - plays well with filenames()
%
% :Optional Inputs:
%
%   **fmri:**
%        Indicates construction of canlab_dataset object using 'fmri'
%        code. Suppresses overwrite warnings specific to 'fmri' inputs.
%
%
% :Outputs:
%
%   **dat:**
%        canlab_dataset object with uploaded values
%
%
% :Examples:
% ::
%
%    % To output a file into a raw fmri dataset
%
%    DesignFile = fullfile(pwd,'Sample_canlab_dataset_experiment_level.xlsx');
%    SubjectFiles = filenames(fullfile(pwd,'Sample_canlab_dataset_subject*.xlsx'));
%    dat = canlab_dataset('fmri');
%    dat = read_from_excel(dat,DesignFile,SubjectFiles,'fmri');
%




fmri = 0;
warnflag = 1;
newvars = {};

for i = 1:numel(varargin)
    switch varargin{i}
        case 'fmri'
            fmri = 1;
        otherwise
            warning(['Unknown argin: ',varargin{i}]);
    end
end

[~,~,Xdata] = xlsread(ExperimentFileName); %raw eXperimentdata

colheaders = Xdata(1,:);
data = Xdata(2:end,:);

dat.Subj_Level.id = data(:,strcmp(colheaders,'id')); 
nsubj = numel(dat.Subj_Level.id);

% choose rows that are not NaNs and the column identified by 'names'.
% Similar behavior for other vars.
wh_col = strcmp(colheaders,'names');
mydat = data(~cellfun(@nanremove,data(:, wh_col)), wh_col);
if iscolumn(mydat), mydat = mydat'; end % enforce row
dat.Subj_Level.names = mydat;

mydat = data(~cellfun(@nanremove,data(:,strcmp(colheaders,'units'))),strcmp(colheaders,'units'));
if iscolumn(mydat), mydat = mydat'; end % enforce row
dat.Subj_Level.units = mydat;

dat.Subj_Level.descrip = data(~cellfun(@nanremove,data(:,strcmp(colheaders,'descrip'))),strcmp(colheaders,'descrip'));


%loop over names and import data into canlab_dataset object
blankdata = nan(nsubj,1);
blanktext = cell(nsubj,1);
for i = 1:nsubj
    blanktext{i,1} = '';
end


for i = 1:length(dat.Subj_Level.names)
    datacol = data(:,strcmp(colheaders,dat.Subj_Level.names{i}));
    if sum(cellfun(@ischar, datacol)) == 0 %numeric data
        dat.Subj_Level.data(:,i) = cellfun(@(x) x, datacol);
        dat.Subj_Level.textdata(:,i) = blanktext;
        dat.Subj_Level.type{i,1} = 'numeric';
    else %text data
        %turn all data to text to avoid potential loss of data
        notchar = find(~cellfun(@ischar,datacol));
        for j = 1:numel(notchar)
            if isnan(datacol{notchar(j)})
                datacol{notchar(j)} = '';
            else
                datacol{notchar(j)} = num2str(datacol{notchar(j)});
            end
        end
        dat.Subj_Level.data(:,i) = blankdata;
        dat.Subj_Level.textdata(:,i) = datacol;
        dat.Subj_Level.type{i,1} = 'text';
    end
end

% Loop over subject files and import Event Level data
for subjidx = 1:length(SubjectFileList)
    
    [~,~,subdata] = xlsread(SubjectFileList{subjidx});
    subcolheaders = subdata(1,:);
    subdata = subdata(2:end,:);
    
    if fmri
        eventidx = numel(dat.Event_Level.names)+1; %used for non-standard entries into dataset
        for i = 1:numel(subcolheaders)
            switch subcolheaders{i}
                case {'SessionNumber' 'RunName' 'RunNumber' 'TaskName' 'TrialNumber' 'EventName' 'EventOnsetTime' 'EventDuration'}
                    dat = update_field(dat,subdata,i,subjidx, subcolheaders);
                    
                otherwise
                    [dat,warnflag,eventidx, newvars] = write_new_field(dat,subdata,warnflag,i,eventidx,subjidx, newvars, subcolheaders);
            
            end %switch
        end
        
    elseif subjidx == 1 %not a standard fmri canlab_dataset object & first subject run
        
        eventidx = numel(dat.Event_Level.names)+1; 
        for i = 1:numel(subcolheaders)
            [dat,warnflag,eventidx, newvars] = write_new_field(dat,subdata,warnflag,i,eventidx,subjidx, newvars, subcolheaders);
        end 
        
    else %a different subject - column names have already been initialized. 
        eventidx = numel(dat.Event_Level.names)+1; %used for non-standard entries into dataset - these should not exist at this point and function will correctly error if used.
        warnflag = 0; %If it were a problem it would already trigger at subjidx == 1. This is set in case building a new canlab_dataset object, at which point warnflag would be silly.
        for i = 1:numel(subcolheaders)
            [dat,warnflag] = write_new_field(dat,subdata,warnflag,i,eventidx,subjidx, newvars, subcolheaders); %eventidx and newvars don't need to be updated now - error is thrown if attempted.
        end        

    end %fmri if
end %for loop over subjects

end % main function

function [dat,warnflag,eventidx,newvars] = write_new_field(dat,subdata,warnflag,idx,eventidx,subjidx,newvars,subcolheaders)
%this function 

if any(strcmp(subcolheaders{idx},dat.Event_Level.names))
    if warnflag && ~any(strcmp(subcolheaders{idx},newvars)) %Haven't been warned, and overwritten var hadn't just been created
        reply = input(['Warning, variable name ' subcolheaders{idx} ' exists. Overwrite (Y/N)? Selecting Y will propagate choice for *all* subsequent queries (SO BE SURE): '], 's');
        if strcmp(reply,'Y')
            fprintf('Continuing.\n')
            %variable name exists, put it in the correct spot
            dat = update_field(dat,subdata,idx,subjidx, subcolheaders);
            warnflag = 0;
        else
            error('You have chosen not to overwrite the canlab_dataset object. No changes made.') 
        end
    else
        dat = update_field(dat,subdata,idx,subjidx, subcolheaders);
    end
else
    %variable is new, output error if this is not the first subject
    if subjidx ~= 1
        error(['Tried to illegally generate Event_Level name: ' subcolheaders{idx} '; This did not exist in previous data objects.' ...
            ' Best practice is to have the same variables in all datafiles. Fill it with NaNs if that variable is not used for a given Event.']);
    end
    dat.Event_Level.names{eventidx} = subcolheaders{idx};
    newvars{end+1} = subcolheaders{idx};
    eventidx = eventidx + 1;
    dat = update_field(dat,subdata,idx,subjidx,subcolheaders);
end
end %write_new_field function

function dat = update_field(dat,subdata,idx,subjidx,subcolheaders)
%adds data to field in canlab_dataset object
loc = find(strcmp(dat.Event_Level.names,subcolheaders{idx}));

ntrials = size(subdata,1);
blankdata = nan(ntrials,1);
blanktext = cell(ntrials,1);
for i = 1:ntrials
    blanktext{i,1} = '';
end

datacol = subdata(:,idx);
if sum(cellfun(@ischar, datacol)) == 0 %numeric data
    dat.Event_Level.data{subjidx}(:,loc) = cellfun(@(x) x, datacol);
    dat.Event_Level.textdata{subjidx}(:,loc) = blanktext;
    dat.Event_Level.type{loc} = 'numeric';
else %text data
    %turn all data to text to avoid potential loss of data
    notchar = find(~cellfun(@ischar,datacol));
    for j = 1:numel(notchar)
        if isnan(datacol{notchar(j)})
            datacol{notchar(j)} = '';
        else
            datacol{notchar(j)} = num2str(datacol{notchar(j)});
        end
    end
    dat.Event_Level.data{subjidx}(:,loc) = blankdata;
    dat.Event_Level.textdata{subjidx}(:,loc) = datacol;
    dat.Event_Level.type{loc} = 'text';
    dat.Event_Level.units{loc} = 'text';
end

end %update_field function



