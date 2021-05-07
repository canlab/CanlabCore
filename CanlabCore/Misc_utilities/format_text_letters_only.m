function [str_array, warnings] = format_text_letters_only(str_array, varargin)
% Make sure string or cell vector of strings is valid as a variable name or table column
% Replaces all non-letter characters in a string or cell array of strings
% with '_'.  Good for saving variables with str_array, and other things.
%
% str_array = format_text_letters_only(str_array)
%
% Optional arguments:
% 'numbers'
% 'squeeze'  replace with no character instead of _
% 'cleanup'  Replace leading and trailing _'s with no char, and replace double __ with no char
%            Good for formatting for use as variable names!
%
% str_array = format_text_letters_only(str_array, 'numbers')  % numbers ok too, do not replace
%
% Best for variable names and table names:
% str_array = format_text_letters_only(str_array, 'numbers', 'cleanup', 'squeeze', 'underscore_ok');



search_criteria = 'alpha';
replace_with = '_';
docleanup = false;
warnings = {};
underscore_ok = false;

if any(strcmp(varargin, 'numbersok')) || any(strcmp(varargin, 'numbers'))
    % numbers are ok
    search_criteria = 'alphanum';
end

if any(strcmp(varargin, 'squeeze'))
    % replace with no character
    replace_with = '';
end

if any(strcmp(varargin, 'cleanup'))
    % replace with no character
    docleanup = true;
end

if any(strcmp(varargin, 'underscore_ok'))
    % replace with no character
    underscore_ok = true;
end


% -----------------------------------------------------------------

% Handles cell and non-cell input
is_letter = isstrprop(str_array, search_criteria);

if iscell(is_letter)
    
    for i = 1:length(str_array)
        
        wh = ~is_letter{i};
        
        if underscore_ok
            wh(str_array{i} == '_') = false;
        end
        
        if any(wh)
            warnings{end + 1} = 'Warning: String has special characters that will be dropped.';
        end
        
        if ~isempty(replace_with)
            str_array{i}(wh) = replace_with;
        else
            str_array{i}(wh) = '';
        end
        
        if docleanup
            
            str_array{i} = clean_up(str_array{i});
            
        end
        
        if isempty(str_array{i})
            warnings{end + 1} = 'Warning: After dropping special characters, string is empty. Invalid string.';
        end
        
    end
    
else % not a cell
    
    wh = ~is_letter;
    
    if underscore_ok
        wh(str_array{i} == '_') = false;
    end
        
    if any(wh)
        warnings{end + 1} = 'Warning: String has special characters that will be dropped.';
    end
        
    if ~isempty(replace_with)
        str_array(wh) = replace_with;
    else
        str_array(wh) = '';
    end
    
    if docleanup
        
        str_array = clean_up(str_array);
        
    end
    
    if isempty(str_array)
        warnings{end + 1} = 'Warning: After dropping special characters, string is empty. Invalid string.';
    end
    
end

end % main function



function str_array = clean_up(str_array)

% Recursively eliminate double __

do_double = true;
while do_double
    
    str_array = strrep(str_array, '__', '_');
    
    if ~any(strfind(str_array, '__')), do_double = false;
        
    end
    
end
    
% Get rid of leading and trailing non-chars, including _
%
not_letter = ~isstrprop(str_array, 'alpha');  % get rid of leading/trailing numbers always.

% remove trailing
try not_letter(2:end-1) = false; catch, end  % try because may be too short

% remove leading by appending
if not_letter(1)
    str_array = ['X_' str_array];
end

% replace leading char with letter if not letter


end % function
