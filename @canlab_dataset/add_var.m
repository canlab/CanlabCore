% Not complete yet. Please edit me
% 
% Function for adding a variable to a dataset in a systematic way.
% - Checks IDs of subjects to make sure data is added in the correct order.
% - Values for missing data are coded as missing values, as specified in dat.Description.Missing_Values
% - Handles Subject_Level or Event_Level data

varname = 'ValenceType';
var_to_add = SETUP.data.X;

% get indices
[subjn, ia, ib] = intersect(dat.Subj_Level.id, subjname2, 'stable');

wh_subjs = false(length(dat.Subj_Level.id));
wh_subjs(ia) = true;

n = length(dat.Subj_Level.id);
    

% set missing value
missingval = dat.Description.Missing_Values;
if ischar(missingval), missingval = str2num(missingval);
    

nmissing = length(dat.Subj_Level.id) - length(ia);
fprintf('Adding variable %s:\n\tSubjects with missing data: %d\n', varname, nmissing);

nextra = length(subjname2) - length(ia);
fprintf('\tSubjects not in dataset (will not be added): %d\n', nextra);

is_subj_level = iscell(var_to_add);

% Subject level
% simplest case
% --------------------------------------------------------------------
if is_subj_level
    
    newvar = repmat(missingval, n, 1);
    
    newvar(ia) = var_to_add(ib);
    
    return
    
end

% Event level
% --------------------------------------------------------------------

for i = 1:n
    
    
    
end


str2num('NaN')

