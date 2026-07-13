function s = labeltable2struct(T)
% labeltable2struct Convert a label_table (MATLAB table) to a struct array.
%
% Inverse of struct2labeltable: turns the fmri_surface_data.label_table MATLAB
% table (variables key, name, rgba) back into the struct array (.key, .name,
% .rgba) that the native CIFTI/GIFTI writers expect. Passes an existing struct
% array through unchanged, so callers can accept either format.
%
% :See also: struct2labeltable, write, apply_parcellation

if isstruct(T)
    s = T;
    return
end

if isempty(T) || height(T) == 0
    s = struct('key', {}, 'name', {}, 'rgba', {});
    return
end

n = height(T);
s = struct('key', cell(1, n), 'name', cell(1, n), 'rgba', cell(1, n));
for i = 1:n
    s(i).key = T.key(i);
    s(i).name = char(T.name(i));
    s(i).rgba = T.rgba(i, :);
end
end
