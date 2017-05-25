function idx = cellregexp(cell_in,regexpstring)


    idx = logical(1 - cellfun(@isempty, regexp(cell_in, regexpstring)));
end