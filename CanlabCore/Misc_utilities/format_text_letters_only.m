function str_array = format_text_letters_only(str_array, varargin)
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

search_criteria = 'alpha';
replace_with = '_';
docleanup = false;

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

% -----------------------------------------------------------------

% Handles cell and non-cell input
is_letter = isstrprop(str_array, search_criteria);

if iscell(is_letter)
    
    for i = 1:length(str_array)
        
        wh = ~is_letter{i};
        
        str_array{i}(wh) = replace_with;
        
        if docleanup
            
            str_array{i} = clean_up(str_array{i});
            
        end
    end
    
else % not a cell
    
    wh = ~is_letter;
    
    str_array(wh) = replace_with;
    
    if docleanup
        
        str_array = clean_up(str_array);
        
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
