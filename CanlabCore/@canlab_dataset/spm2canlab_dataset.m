function obj = spm2canlab_dataset(obj, subject, spm)
% Extract Event_Level data from subjects' SPM.mat files to add data to
% canlab_dataset object. 
%
% :Usage:
% ::
%
%    obj = spm2canlab_dataset(obj, subject, spm)
%
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2015  Wani Woo
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
%   **obj:**
%        Canlab_dataset object (see canlab_dataset)
%
%   **subject:**
%        Subject list (it could be one subject [in a string 
%        format], or it could be multiple subjects in cell array)
%
%   **spm:**
%        This could be loaded SPM (struct), or one path for one 
%        subject's SPM.mat file (string), or multiple loaded SPM or
%        paths in cell array
%
% :Outputs:
%
%   **obj:**
%        Canlab_dataset object with new data
%
%
% :Examples: 
% ::
%
%    subj = {'dpsp002','dpsp003'};
%    spm = {'dpsp002_SPM.mat', 'dpsp003_SPM.mat'};
%
%    D = canlab_dataset; % if D doesn't exist yet
%    D = spm2canlab_dataset(D, subj, spm);
%
%
% :See also:
%   canlab_dataset
%   spm_mat2batchinput
%

% ..
%    Programmers' notes:
%       After running this, please check D.Subj_Level.data and
%       D.Event_level.data to see if there are NaNs. NaNs could be there when 
%       information cannot be extracted from the given SPM.mat files. 
% ..


subj = parseSUBJ(subject); % Check 1: parsing subjects

% Check 2: parsing SPM input
SPM = parseSPM(spm);

% Adding data to the canlab_dataset obj
for subj_i = 1:numel(subj)
    
    % check if the subject exists in obj
    [isexist_sub, wh_subj] = exist_data(obj.Subj_Level.id, subj{subj_i});
    
    % if the subject is new, add the subject to obj.Subj_Level
    if ~isexist_sub
        
        obj.Subj_Level.id{wh_subj} = subj{subj_i};
        
        if ~isempty(obj.Subj_Level.names)
            
            str = ['D.Subj_Level.data for subject ' subj{subj_i} ' are now NaNs. Please enter the data.'];
            warning(str);
            
            obj.Subj_Level.data(wh_subj,:) = NaN(1,size(obj.Subj_Level.data,2));
            
        end
    end
    
    % get information from SPM
    out = spm2canlab_dataset_format(SPM{subj_i});
    names = out.names;
    
    for fields_i = 1:numel(names)
        
        % check the same information exists already
        [isexist_data, wh_col] = exist_data(obj.Event_Level.names, names{fields_i});
        
        % if the info is new, the existing subjects data for that column should be NaNs.
        if ~isexist_data && numel(obj.Subj_Level.id) ~= 1
            for i = 1:numel(obj.Event_Level.data)
                obj.Event_Level.data{i}(:,wh_col) = NaN(size(obj.Event_Level.data{i},1),1);
            end
            warning('For existing subjects, please check NaNs in D.Event_Level.data.');
        end
        
        % add the data
        obj.Event_Level.names{wh_col} = names{fields_i};
        obj.Event_Level.type{wh_col} = out.type{fields_i};
        obj.Event_Level.descrip{wh_col} = out.descrip{fields_i};
        obj.Event_Level.data{wh_subj}(:,wh_col) = out.data(:,fields_i);
        
        wh_cols(fields_i) = wh_col;
            
    end
    
    % if the subject is new and there were existing data, the data that the
    % new subject doesn't have should be NaNs
    if ~isexist_sub
        wh_old = ~ismember(1:size(obj.Event_Level.data{wh_subj},2),wh_cols);
        if any(wh_old)
            obj.Event_Level.data{wh_subj}(:,wh_old) = ...
                NaN(size(obj.Event_Level.data{wh_subj}(:,wh_old)));
            warning('For new subjects, please check NaNs in D.Event_Level.data.');
        end
    end
    
end

end

% ============ SUB-FUNCTIONS ============ 

function subj = parseSUBJ(subject)

if ischar(subject)
    subj{1} = subject;
else
    subj = subject;
end

end

function SPM = parseSPM(spm)

if ischar(spm)
    temp = load(spm);
    SPM{1} = temp.SPM; 
else
    if iscell(spm)
        for i = 1:numel(spm)
            if ischar(spm{i})
                temp = load(spm{i});
                SPM{i} = temp.SPM;
            elseif isstruct(spm{i})
                SPM{i} = spm{i};
            else
                error('SPM is unknown format.');
            end
        end
    elseif isstruct(spm)
        SPM{1} = spm;
    else
        error('SPM is unknown format.');
    end
end

end

function [isexist, wh] = exist_data(obj_field, field)

isexist = any(strcmp(obj_field, field));

if isexist
    wh = find(strcmp(obj_field, field));
else
    wh = numel(obj_field) + 1;
end

end

% ---- MAJOR SUBFUNCTION ---- 

function out = spm2canlab_dataset_format(SPM)

% Extract Event_Level data from one subject's SPM.mat file to enter data to
% canlab_dataset object. This function is called by "spm2canlab_dataset".
%
% Usage:
% -------------------------------------------------------------------------
% out = spm2canlab_dataset_format(SPM)
%
% Inputs:
% -------------------------------------------------------------------------
% SPM            SPM saved in SPM.mat of the model directory. 
%
% Outputs:
% -------------------------------------------------------------------------
% out.names      Event_Level names in cell array 
%                (RunNumber, EventName, Onsets, Durations)
% out.descrip    Event_Level descript in cell array
% out.type       Event_Level variable types of data
% out.data       Event_Level data
%
% See also spm2canlab_dataset, canlab_dataset, spm_mat2batchinput
%
% Examples: 
% -------------------------------------------------------------------------
% load SPM.mat;
% out = spm2dataset(SPM);
%
% % output:
% %   names: {'RunNumber'  'EventNames'  'Onsets'  'Durations'}
% % descrip: {'Run 1-8'  {12x1 cell}  'Onsets'  'Duration'}
% %    type: {'Categorical'  'Categorical'  'Continuous'  'Continuous'}
% %    data: [192x4 double]
%
% -------------------------------------------------------------------------
% Copyright (C) 2015  Wani Woo

% Programmers' notes:


% Parsing SPM to know the total number of Event types

[~, ~, ~, names] = spm_mat2batchinput(SPM);

% 1. Event_Level variable names
out.names = {'RunNumber','EventNames','Onsets','Durations'};

% 2. Event_Level descript
out.descrip{1} = ['Run 1-' num2str(numel(SPM.Sess))];
for i = 1:numel(unique(names))
    out.descrip{2}{i,1} = [num2str(i) ': ' deblank(names{i})];
end
out.descrip{3} = 'Onsets';
out.descrip{4} = 'Duration';

% 3. Event_Level variable type
out.type{1} = 'Categorical';
out.type{2} = 'Categorical';
out.type{3} = 'Continuous';
out.type{4} = 'Continuous';

% 4. Event_Level variable type
out.data = [];
for run_i = 1:numel(SPM.Sess)
    
    row_i = 0;
    for eventtype_i = 1:numel(SPM.Sess(run_i).U)
        
        EventName = find(strcmp(unique(names), SPM.Sess(run_i).U(eventtype_i).name));
        
        for event_i = 1:numel(SPM.Sess(run_i).U(eventtype_i).ons)
            
            row_i = row_i+1;
            Durations = SPM.Sess(run_i).U(eventtype_i).dur(event_i);
            Onsets = SPM.Sess(run_i).U(eventtype_i).ons(event_i);
            Ddata_run(row_i,:) = [run_i EventName Onsets Durations];
            
        end
    end
    
    % sorting events using onset information
    [~, idx] = sort(Ddata_run(:,3));
    out.data = [out.data; Ddata_run(idx,:)];
    
end

end
