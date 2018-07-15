function is_empty = isempty(obj)


is_empty = length(obj) == 0 || isempty(obj(1).XYZ) || length(obj(1).XYZ) == 0;


end

