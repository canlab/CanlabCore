function T = struct2labeltable(s)
% struct2labeltable Convert a CIFTI/GIFTI label struct array to a MATLAB table.
%
% The native readers return a label table as a struct array with fields
% .key (integer), .name (char), .rgba (1x4). The fmri_surface_data.label_table
% property stores this as a MATLAB table with variables key (double), name
% (string), rgba (Nx4 double). This helper converts struct -> table (and passes
% an existing table through unchanged).
%
% :See also: labeltable2struct, fmri_surface_data, canlab_read_cifti

if istable(s)
    T = s;
    return
end

if isempty(s)
    T = table('Size', [0 3], 'VariableTypes', {'double', 'string', 'double'}, ...
        'VariableNames', {'key', 'name', 'rgba'});
    return
end

n = numel(s);
key = reshape([s.key], [], 1);
name = strings(n, 1);
rgba = nan(n, 4);
for i = 1:n
    name(i) = string(s(i).name);
    c = s(i).rgba(:)';
    rgba(i, 1:numel(c)) = c;
end
T = table(key, name, rgba, 'VariableNames', {'key', 'name', 'rgba'});
end
