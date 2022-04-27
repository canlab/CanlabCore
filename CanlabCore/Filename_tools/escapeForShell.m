function escapedString = escapeForShell(string)
% Returns converted string for use with the current OS's shell. Hence,
% it's operation will be different, depending on where it's executed.
%
% :Usage:
% ::
%
%     [escapedString] = escapeForShell(string)
%
    % Michael: This function wouldn't accept an actual string so here I
    % converted the string to a character array.
    if(isstring(string))
        string = convertStringsToChars(string);
    end
    
    if(ischar(string))
        for i=1:size(string,1)
            escapedString{i} = efsReplace(string(i,:));
        end
        escapedString = char(escapedString);
    elseif(iscellstr(string))
        for i=1:length(string)
            escapedString{i} = efsReplace(string{i});
        end
    else
        error('Input is not a string or cellstr.');
    end
end

function e = efsReplace(s)
    e = s;
    if(~ispc())
        e = strrep(e, '\', '\\');
    end
    e = strrep(e, '"', '\"');
    e = strrep(e, '''', '\''');
    e = strrep(e, ' ', '\ ');
    e = strrep(e, '(', '\(');
    e = strrep(e, ')', '\)');
end
