% Create a canlab_dataset object from a Matlab Table object
%
% The canlab_dataset object provides a standard data format, and allows the use of more methods,
% such as writing the data object to text for backup/export to other
% software, plotting, analysis tools, and more.
%
% Input ::
%
%   StudyTable          - n x d table object of text and numeric data to 
%                           store in dataset object.
%
%   'labels'            - (optional arg) followed by a n x 1 vector of 
%                           subject labels. If provided DAT.Event_Level 
%                           will be populated with data from StudyTable as 
%                           well, with one cell per subject.
%
% Output ::
% 
%   DAT                 - canlab_dataset object
%
function DAT = table2canlab_dataset(StudyTable, varargin)

labels = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'labels'
                if length(varargin{i+1}) ~= size(StudyTable,1);
                    error('''labels'' length does not match StudyTable rows.');
                else
                    labels = varargin{i+1};
                end
            otherwise
                warning('Did not understand argument %s. Ignoring.',varargin{i});
        end
    end
end

% Create canlab_dataset object
% -------------------------------------------------------------------------

DAT = canlab_dataset;

% Fill in Experiment-level information
% ------------------------------------------------------------------------

DAT.Description.Experiment_Name = [];
DAT.Description.Missing_Values = NaN;

DAT.Description.Subj_Level = {};

DAT.Description.Subj_Level{1} = 'Imported from structure' ;

% Names of subject-level data we want to store
% -------------------------------------------------------------------------

DAT.Subj_Level.names = StudyTable.Properties.VariableNames;

% Data from table
% -------------------------------------------------------------------------

studyCell = table2cell(StudyTable);
studyRow = studyCell(1,:);
numeric = ones(size(studyRow));
character = ones(size(studyRow));
for i = 1:length(studyRow)
    if ~isnumeric(studyRow{i});
        numeric(i) = 0;
    end
    if ~ischar(studyRow{i})
        character(i) = 0;
    end
end

if any(numeric==0) % upload numeric data to rows where available, nans elsewhere, upload text to textdata.
    warning(sprintf('Table is a cell array. Normally entire table gets assigned to Subj_Level.data. Will additionally assign text data to Subj_Level.textdata. Unexpected things may happen downstream.'));
    if sum(numeric) + sum(character) ~= size(studyRow,2)
        error('Non numeric and non character data is no supported by table2canlab_dataset. Please remove nontext and nonnumeric classes from table and try again.');
    end
    DAT.Subj_Level.data = zeros(size(studyCell));
    DAT.Subj_Level.data(:) = nan;
    DAT.Subj_Level.data(:,logical(numeric)) =  cell2mat(studyCell(:,logical(numeric)));
    
    DAT.Subj_Level.textdata = cell(size(studyCell));
    DAT.Subj_Level.textdata(:,logical(character)) = studyCell(:,logical(character));
else
    DAT.Subj_Level.data = table2array(StudyTable);
end
DAT.Subj_Level.descrip{1} = 'Imported from Table object. Add descriptions.';
DAT.Subj_Level.type = repmat({'numeric'}, 1, length(DAT.Subj_Level.names));

% Subject ID - put in cell array of strings

DAT.Subj_Level.id = {};
nobs = size(DAT.Subj_Level.data, 1);

for i = 1:nobs
    DAT.Subj_Level.id{i} = sprintf('Obs%3.0f', i);
end

%% populate Event_Level property
if ~isempty(labels)
    DAT.Event_Level.names = DAT.Subj_Level.names;
    DAT.Event_Level.descrip = 'Imported from table.';
    DAT.Event_Level.units = repmat({'Imported from table.'},1,length(DAT.Event_Level.names));
    DAT.Event_Level.type = cell(1,length(numeric));
    DAT.Event_Level.type(logical(numeric)) = repmat({class(DAT.Subj_Level.data)},1,sum(numeric));
    DAT.Event_Level.type(logical(character)) = repmat({'char'},1,sum(character));
    uniq_lbls = unique(labels);
    n = length(uniq_lbls);
    for i = 1:n
        this_lbl = uniq_lbls(i);
        lbl_idx = find(labels == this_lbl); 
        DAT.Event_Level.data{end+1} = DAT.Subj_Level.data(lbl_idx,:);
        if ~isempty(DAT.Subj_Level.textdata)
            DAT.Event_Level.textdata{end+1} = DAT.Subj_Level.textdata(lbl_idx,:);
        end
    end
end
    
end % function