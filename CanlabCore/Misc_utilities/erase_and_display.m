function erase_and_display(old_str, new_str)
% Backspaces the number of chars of the length of old_str, then prints new_str

    fprintf([repmat('\b', 1, length(old_str)) new_str]);
end
