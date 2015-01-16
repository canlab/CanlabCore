% function erase_string(str)
%   Backspaces the number of chars of the length of str

function erase_string(str)
    fprintf(1, repmat('\b', 1, length(str)));
end