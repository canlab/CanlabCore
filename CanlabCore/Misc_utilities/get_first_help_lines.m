function helptext = get_first_help_lines(functionname, maxlines)
% Returns the first n lines of the help for a function in a cell array
%
% :Usage:
% ::
%
%     helptext = get_first_help_lines(functionname, [maxlines])
%
%     helptext = get_first_help_lines(functionname, maxlines)
%
% :Examples:
% ::
%
%    helptext = get_first_help_lines('fmri_data.apply_mask', 5);
%    char(helptext{:})
%
% ..
%    Tor Wager, Oct 2012
% ..

if nargin < 2
    maxlines = 5;
end


h = help(functionname);

isnewline = find(h == sprintf('\n'));

nlines = min(maxlines, length(isnewline));

h = h(1:isnewline(nlines));

isnewline = [1 isnewline];

for i = 1:nlines

    helptext{i} = h(isnewline(i):isnewline(i+1)-1);
   
    helptext{i}(helptext{i} == sprintf('\n')) = [];
    
end

end % function
