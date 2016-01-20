function cellarray = remove_str_from_cell(cellarray,str)

    wh = cellfun(@isempty, strfind(cellarray, str));
    cellarray = cellarray(wh);

end
