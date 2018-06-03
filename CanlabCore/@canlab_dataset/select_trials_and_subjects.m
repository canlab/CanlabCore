function [trial_indx, subject_indx, datavar_matrix, datavar_cell] = select_trials_and_subjects(DAT, conditioning_text_var_name, text_pattern, varargin)
% Filter an event-level variable:
% - Return a logical vector of which events and subjects have values on a
% string variable (conditioning_text_var_name) that match a string pattern
% (text_pattern)
%
% - If the name of a second variable is entered as an optional argument,
% extract data from that variable, selecting only events that match the
% filter. [var name to extract data] in the command string below is the
% optional argument.
%
% conditioning_text_var_name is the name of an Event_Level variable in your dataset
% text_pattern is a string pattern to match with regexp
%
% [trial_indx, subject_indx, datavar_matrix, datavar_cell] = select_trials_and_subjects(DAT, conditioning_text_var, text_pattern, [var name to extract data])
%
% datavar_matrix is matrix of subjects x trials with extracted data
% datavar_cell is cell array of length subjects containing column vectors with extracted data

[datavar_matrix, datavar_cell] = deal([]);

[~, block] = get_var(DAT, conditioning_text_var_name);


% This function turns a string into a vector of which trials match that string.
subj_indx = @(subjcell) ~cellfun(@isempty, regexp(subjcell, text_pattern, 'match'));

% This function tests for all empty or nans in a cell.
nanfun = @(block) all(cellfun(@(x) all(isempty(x) | isnan(x)), block{1}));
%nanfun(block(1))

% This function tests whether any trials match the string
anyfun = @(x) any(subj_indx(x{1}));

% This function operates on one cell, for use with cellfun wrapper
onesubjfun = @(x) subj_indx(x{1});

% subj_indx(block{5})
% onesubjfun(block(5))

% doesn't work for unknown reasons - actually, can't handle nans.
% indic = cellfun(onesubjfun, block, 'UniformOutput', false)
% indic = cellfun(@(x) any(x), block)

clear trial_indx subject_indx
for i = 1:length(block)
    
    if nanfun(block(i))
        trial_indx{i} = false(size(block{i}));
        subject_indx(i, 1) = false;
    else
        trial_indx{i} = subj_indx(block{i});
        subject_indx(i, 1) = anyfun(block(i));
        
    end
    
end

% retrieve target data variable, removing NaNs and concatenating
if length(varargin) > 0
    
    var_name_to_extract = varargin{1};
    
    [~, datavar_cell, ~, ~, ~, textflag] = get_var(DAT, var_name_to_extract);

    
    
    for i = 1:length(trial_indx)
       
        datavar_cell{i} = datavar_cell{i}(trial_indx{i});
        
    end
    
    % pad to longest and concatenate into matrix
    
    %[~, textflag] = get_varlevel(DAT, var_name_to_extract);
    
    d = pad_cells_with_nan(datavar_cell, textflag);
    
    datavar_matrix = cat(2, d{:})';
    
end


end % function




function d = pad_cells_with_nan(d, textflag)

slen = max(cellfun(@length, d)); % max length for any subject

refvec = ones(slen, 1);

for i = 1:length(d)
    
    if textflag % text var
        reftext = cell(1, abs(length(d{i}) - length(refvec)))';
        [reftext{:}] = deal('');
        
        d{i} = [d{i}; reftext];
        
    else % numeric
        d{i} = padwithnan(d{i}, refvec, 1);
    end
    
end

end


