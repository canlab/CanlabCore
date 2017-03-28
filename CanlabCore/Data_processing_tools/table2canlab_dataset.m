function DAT = table2canlab_dataset(StudyTable)
% Create a canlab_dataset object from a Matlab Table object
%
% The canlab_dataset object provides a standard data format, and allows the use of more methods,
% such as writing the data object to text for backup/export to other
% software, plotting, analysis tools, and more.

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

DAT.Subj_Level.data = table2array(StudyTable);
DAT.Subj_Level.descrip{1} = 'Imported from Table object. Add descriptions.';
DAT.Subj_Level.type = repmat({'numeric'}, 1, length(DAT.Subj_Level.names));

% Subject ID - put in cell array of strings

DAT.Subj_Level.id = {};
nobs = size(DAT.Subj_Level.data, 1);

for i = 1:nobs
    DAT.Subj_Level.id{i} = sprintf('Obs%3.0f', i);
end
    
end % function