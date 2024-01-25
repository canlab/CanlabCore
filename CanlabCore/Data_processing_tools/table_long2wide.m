function [out_table, out_matrix, colnames] = table_long2wide(dat, outcomevarname, varname1, varname2, sidname)
% Converts long-format table to wide-format, averaging values within-participant for each combo of conditions in varname1 and varname2  
%
% [out_table, out_matrix, colnames] = table_long2wide(dat, outcomevarname, varname1, varname2, sidname)
%
% dat               A table object
% outcomevarname:   Name of table variable (column) that defines the outcome to extract; numeric
% varname1:         Name of table variable (column) that defines a repeated measure (condition) for output columns; must be string
% varname2:         Name of table variable (column) that defines a repeated measure (condition) for output columns; must be string
% sidname:          Name of table variable (column) that defines subject ID; must be numeric
% 
% This is a stub without documentation -- developers please flesh this out
% using documentation_template.m for formatting
% And: this could be extended to:
% 1 - handle variables with different formats
% 2 - make it more flexible. e.g., handle 1, 2, or more condition variables
% 3 - Could also have other functions besides mean as an input
% 4 - Could allow entry of names of levels as varargin to order them
%
% Tor Wager, Jan 2024

y = dat.(outcomevarname);
v1 = dat.(varname1);
v2 = dat.(varname2);
svar = dat.(sidname);

funhan = @(x) nanmean(x);

u1 = unique(v1);
u2 = unique(v2);

sid = unique(svar, 'stable');

for i1 = 1:length(u1)

    for i2 = 1:length(u2)

        whcol = length(u2) * (i1 - 1) + i2;

        colnames{whcol} = [u1{i1} '_' u2{i2}];

        for s = 1:length(sid)

            wh = strcmp(v1, u1(i1)) &  strcmp(v2, u2(i2)) & svar == sid(s);

            out_matrix(s, whcol) = funhan(y(wh));

            out_ycell{whcol}(s, 1) =  out_matrix(s, whcol); % for table later
        end

    end

end

out_table = table(out_ycell{:}, 'VariableNames', colnames);

end % function

% 
% 
% tmp = unstack(dat, ratingname, {cuename_string stimnname_string})
% tmp2 = unstack(tmp, ratingname, stimnname_string)
% tmp2 = rowfun(@(x) nanmean(x), tmp, 'InputVariables', {'high_cue'}, 'GroupingVariables', {'src_subject_id'})
