function [DAT, StudyTable] = struct2canlab_dataset(data_struct)
% Take an ad hoc data structure data_struct and return a Table format
% object (StudyTable) and a canlab_dataset object DAT.
%
% - field names become variable names
% - Column vectors of numbers or text becomes Subj_Level data
% - Cell vectors containing column vectors of trial-level numeric data become Event_level data
% - Number of elements in vector stored in each field should match num. observations (subjects) for Subj_level array data
% - Trial-level data should contain same-length vectors for all variables for a given subject
%
% See also struct2table and table2canlab_dataset

N = fieldnames(data_struct);

% Empty structures to sort into
CELLS = struct();
COLUMNS = struct();

for i = 1:length(N)
    
    myvar = data_struct.(N{i});
    
    % Try to enforce columnar orientation
    if ~iscolumn(myvar), myvar = myvar'; end
    
    % If still not a column, skip this (matrix?)
    if ~iscolumn(myvar), continue; end
    
    if iscell(myvar)
        
        if all(cellfun(@length, myvar) == 1) && ~any(cellfun(@ischar, myvar))
            % all single values, really numbers
            myvar = cat(1, myvar{:});
            COLUMNS.(N{i}) = myvar;
        else
            CELLS.(N{i}) = myvar;
        end
        
    else
        % Text or event-level data
        COLUMNS.(N{i}) = myvar;
        
    end
    
end

%% Get sizes and remove odd-sized

CELLS = remove_odd_sized(CELLS);

COLUMNS = remove_odd_sized(COLUMNS);

% ------------------------------------------------------------------------


%% Make table

StudyTable = struct2table(COLUMNS);

%% Initialize canlab_dataset
% Deal with subject-level data in matrix first

DAT = table2canlab_dataset(StudyTable);

%% Add text data

DAT = add_text_data_from_CELLS(DAT, CELLS);

%% Add event-level data from cells

CELLS = remove_text_data_from_CELLS_struct(DAT, CELLS);

N = fieldnames(CELLS);
nobs = length(CELLS.(N{1}));

% Build array of event-level data
% Deal data from each of nobs into a matrix of vars for each obs
XX = cell(1, nobs);

for i = 1:length(N)
    
    myvar = CELLS.(N{i});
    XXnames(1, i) = N(i);
    
    for j = 1:nobs
        
        datatoadd =  myvar{j};
        
        if ~iscolumn(datatoadd), datatoadd = datatoadd'; end
        if ~iscolumn(datatoadd), continue; end
        
        if ~isempty(XX{j}) && size(XX{j}, 1) ~= size(datatoadd, 1)
            continue
        end
        
        XX{j}(:, i) = datatoadd;
        
    end
    
end

DAT.Event_Level.names = XXnames;
DAT.Event_Level.data = XX;
DAT.Event_Level.type = repmat({'numeric'}, 1, length(XXnames));

end % function


function CELLS = remove_odd_sized(CELLS)

N = fieldnames(CELLS);

for i = 1:length(N)
    
    cellsize(i) = size(CELLS.(N{i}), 1);
    
end

cellmode = mode(cellsize);
cellmode = cellmode(1);

wh_omit = cellsize ~= cellmode;

for i = 1:length(N)
    
    if wh_omit(i)
        
        CELLS = rmfield(CELLS, N{i});
        
    end
    
end

end




function DAT = add_text_data_from_CELLS(DAT, CELLS)

nobs = 0;
if ~isempty(DAT.Subj_Level.data), nobs = size(DAT.Subj_Level.data, 1); end

N = fieldnames(CELLS);

for i = 1:length(N)
    
    cellsize = size(CELLS.(N{i}), 1);
    
    % Skip if Subj_level nobs is defined and this does not match
    if nobs && cellsize ~= nobs
        continue
    end
    
    istext = all(cellfun(@ischar, CELLS.(N{i})));
    
    if istext
        
        nextvar = max(size(DAT.Subj_Level.data, 2), size(DAT.Subj_Level.textdata, 2));
        
        DAT.Subj_Level.textdata(:, nextvar) = CELLS.(N{i});
        DAT.Subj_Level.names(nextvar) = N(i);
        DAT.Subj_Level.type(nextvar) = {'text'};
        
    end
    
    
end

end



function CELLS = remove_text_data_from_CELLS_struct(DAT, CELLS)

nobs = 0;
if ~isempty(DAT.Subj_Level.data), nobs = size(DAT.Subj_Level.data, 1); end

N = fieldnames(CELLS);

for i = 1:length(N)
    
    cellsize = size(CELLS.(N{i}), 1);
    
    % Remove if Subj_level nobs is defined and this does not match
    if nobs && cellsize ~= nobs
        
        CELLS = rmfield(CELLS, N{i});
        
    end
    
    istext = all(cellfun(@ischar, CELLS.(N{i})));
    
    if istext
        
        CELLS = rmfield(CELLS, N{i});
        
    end
    
    
end

end


