function CLDAT = add_vars(CLDAT, dat_to_add, wh_level, varargin)
%
% Add a matrix, table (table object), or Excel file of variables to a canlab_dataset object
%
% - Enter matrix data, data in table object format, or Excel file name in dat_to_add
% - Replace existing variables with same names, or add new variables if names are new
% - text data not supported yet.
%
% :Usage:
% ::
%
% CLDAT = add_vars(CLDAT, dat_to_add, wh_level, ['names', names_to_add], ['descrip', descrip_to_add])
%
% :Inputs:
%
%   **wh_level:**
%        Must be 'Subj_Level' or 'Event_Level'
%
%   **dat_to_add:**
%        Matrix or table
%        Num observations must match existing number in CLDAT, unless CLDAT is empty.
%        For tables, variable names are used as names to add unless these
%        are overwritten by a manual entry of names.
%
% :Optional Inputs:
%
%   **'names', names_to_add:**
%        'names' followed by cell array of names to add
%        number of elements of names_to_add must match number of variables
%        added
%        Names is obligatory, and must either be entered as an input
%        argument or be obtained automatically from table object names
%
%   **'descrip', descrip_to_add:**
%        'descrip' followed by cell array of names to add
%        number of elements of names_to_add must match number of variables
%        added
%
%   **'subjectid', subject_id:**
%        'subjectid' followed by string with subject id to add
%        subject_id must match an existing entry in CLDAT.Subj_Level.id
%        This is used only for Event_Level data, otherwise ignored
%        For Event_Level data it is obligatory
%
% :Outputs:
%
%   **CLDAT:**
%        Output canlab_dataset object with variables added
%
% Notes:
% You can add variables manually
%
% Example:
% DesignFile = which('Sample_subject_level_data_to_add.xlsx');
% dat = canlab_dataset();
% dat = add_vars(dat, DesignFile, 'Subj_Level');
% print_summary(dat)
%

%     Author and copyright information:
%
%     Copyright (C) 2017 Tor Wager
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

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

[names_to_add, descrip_to_add]  = deal({});
subject_id = '';

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'names', names_to_add = varargin{i+1}; varargin{i+1} = [];
            case 'descrip', descrip_to_add = varargin{i+1}; varargin{i+1} = [];
                
            case 'subjectid', subject_id = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if exist(dat_to_add, 'file')
    
    myfile = dat_to_add;
    dat_to_add = [];
    
    fprintf('Importing from file: %s\n', myfile);
    
    [~,~,Xdata] = xlsread(myfile); %raw eXperimentdata

    names_to_add = Xdata(1,:);

    filedata = Xdata(2:end,:);
    dat_to_add = NaN .* ones(size(filedata));
    
    for i = 1:length(names_to_add)
        
        mydat = filedata(:, i);
        
        if all(cellfun(@isnumeric, mydat))
            
            mydat = cat(1, mydat{:});
            
            dat_to_add(:, i) = mydat;
            
        else
            % text - add text data import later...
            dat_to_add(:, i) = NaN;
        end
        
        descrip_to_add{i} = ['Imported from ' myfile];
        
    end % names

    
end % file

    
if ~isempty(names_to_add)
    names_to_add = enforce_var_names(names_to_add);
end


% -------------------------------------------------------------------------
% Convert from table if needed
% -------------------------------------------------------------------------

if istable(dat_to_add)
    
    names_to_add = dat_to_add.Properties.VariableNames;
    descrip_to_add = dat_to_add.Properties.VariableDescriptions;
    dat_to_add = table2array(dat_to_add);
    
    %units_to_add = dat_to_add.Properties.VariableUnits;
    
end

n = length(names_to_add);

% -------------------------------------------------------------------------
% Error checking
% -------------------------------------------------------------------------

if length(names_to_add) ~= size(dat_to_add, 2)
    error('Names missing/wrong length? Enter ''names'' followed by cell array of names for each column to add');
end

% Check for ID - and create if needed
if isempty(CLDAT.Subj_Level.id)
    
    CLDAT.Subj_Level.id = {};
    
    nsubj = size(dat_to_add, 1);
    
    for i = 1:nsubj
        CLDAT.Subj_Level.id{i, 1} = sprintf('TEMPORARY Subj %d - ENTER dat.Subj_Level.id in cells!', i);
    end
    
end

switch wh_level
    case 'Subj_Level'
        
        if ~isempty(CLDAT.(wh_level).data) && size(dat_to_add, 1) ~= size(CLDAT.(wh_level).data, 1)
            error(sprintf('Num observations in %s.data field does not match input dat_to_add.'), wh_level)
        end
        
    case 'Event_Level'
        
        if ischar(subject_id) % this should really be char according to format, but accept numbers here too for convenience
            wh_subj = strcmp(subject_id, CLDAT.Subj_Level.id);
        else
            wh_subj = subject_id == cat(1, CLDAT.Subj_Level.id{:});
        end
        
        if ~any(wh_subj)
            error('For Event_Level data, enter ''subjectid'' followed by string of a subject id in CLDAT.Subj_Level.id');
        end
        
        if ~isempty(CLDAT.(wh_level).data{wh_subj}) && size(dat_to_add, 1) ~= size(CLDAT.(wh_level).data{wh_subj}, 1)
            error(sprintf('Num observations in %s.data field does not match input dat_to_add.'), wh_level)
        end
        
        
    otherwise
        error('wh_level input not valid. Must be ''Subj_Level'' or ''Event_Level''');
end

% check for existing names

if strcmp(wh_level, 'Subj_Level')
    % Don't print for event_level as we'll often replace
    for i = 1:n
        
        wh = strcmp(names_to_add{i}, CLDAT.(wh_level).names);
        
        if any(wh)
            % Don't print for event_level as we'll often replace
            fprintf('Var %s exists. Replacing.');
        end
        
    end
    
end


% -------------------------------------------------------------------------
% Add variables
% -------------------------------------------------------------------------

for i = 1:n
    
    wh = strcmp(names_to_add{i}, CLDAT.(wh_level).names);
    
    if any(wh)
        wh_var = find(wh);
    else
        wh_var = length(CLDAT.(wh_level).names) + 1; % Empty var names may be entered. Keep track this way.
    end
    
    switch wh_level
        case 'Subj_Level'
            
            if any(wh)
                CLDAT.(wh_level).data(:, wh_var) = dat_to_add(:, i);
                
            else
                wh_var = length(CLDAT.(wh_level).names) + 1; % Empty var names may be entered. Keep track this way.
                
                CLDAT.(wh_level).data(:, wh_var) = dat_to_add(:, i);
                CLDAT.(wh_level).names(wh_var) = names_to_add(i);
            end
            
        case 'Event_Level'
            
            if any(wh)
                CLDAT.(wh_level).data{wh_subj}(:, wh_var) = dat_to_add(:, i);
                
            else
                
                CLDAT.(wh_level).data{wh_subj}(:, wh_var) = dat_to_add(:, i);
                CLDAT.(wh_level).names(wh_var) = names_to_add(i);
            end
            
            
        otherwise
            error('wh_level input not valid. Must be ''Subj_Level'' or ''Event_Level''');
    end
    
    % Description, if entered
    
    if ~isempty(descrip_to_add) && ~isempty(descrip_to_add{i})
        
        CLDAT.(wh_level).descrip{wh_var, 1} = descrip_to_add{i};
        
    end
    
    % Var type
    CLDAT.(wh_level).type{wh_var} = 'numeric';
    
end % variable to add



end % function



function names_to_add = enforce_var_names(names_to_add)
% Enforce valid names - these will become variable names

wh_nan = cellfun(@(x) any(isnan(x)), names_to_add);
names_to_add(wh_nan) = {'Noname'};

names_to_add = cellfun(@(x) strrep(x, ' ', '_'), names_to_add, 'UniformOutput', false);
names_to_add = cellfun(@(x) strrep(x, ',', '_'), names_to_add, 'UniformOutput', false);
names_to_add = cellfun(@(x) strrep(x, ':', '_'), names_to_add, 'UniformOutput', false);
names_to_add = cellfun(@(x) strrep(x, '!', '_'), names_to_add, 'UniformOutput', false);
names_to_add = cellfun(@(x) strrep(x, '?', '_'), names_to_add, 'UniformOutput', false);
names_to_add = cellfun(@(x) strrep(x, ';', '_'), names_to_add, 'UniformOutput', false);
names_to_add = cellfun(@(x) strrep(x, '=', '_'), names_to_add, 'UniformOutput', false);

end
