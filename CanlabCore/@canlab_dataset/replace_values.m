function D = replace_values(D, varname, function_handle, replace_with)
% Replace values identified with your custom input function with NaN
%
% D = replace_values(D, function_handle to select values to replace, replace_with, varname)
%
% function_handle is an anonymous function that selects values to replace.
% It should return a logical vector of true/false values. e.g.,  @(x) x <= 0
%
% e.g.,
% Take the variable 'ratings' in canlab_dataset object D
% Replace all ratings <= 0 with NaN
%
% D = replace_values(D, 'ratings', @(x) x <= 0, NaN)

[dat, datcell, wh_level, descrip, wh_indx] = get_var(D, varname);  % to list variables

switch wh_level
    case 1  % Subject
        
        wh = function_handle(dat); % select observations to replace, logical
        
        dat(wh) = replace_with;
        
        D.Subj_Level.data(:, wh_indx) = dat;
        
    case 2  % Event
        
        n = length(datcell);  % for each subject
        
        for i = 1:n
            
            % Do not use datcell because NaNs are removed there
            
            mydat = D.Event_Level.data{i}(:, wh_indx);
            
            wh = function_handle(mydat); % select observations to replace, logical
            
            mydat(wh) = replace_with;
            
            D.Event_Level.data{i}(:, wh_indx) = mydat;
            
        end % subject
        
end

end % function