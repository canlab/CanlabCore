function isdiff = compare_space(obj, obj2)
% compare_space Compare the grayordinate spaces of two fmri_surface_data objects.
%
% :Usage:
% ::
%     isdiff = compare_space(obj, obj2)
%
% Surface analogue of image_vector.compare_space, preserving the same integer
% return contract so callers (cat, apply_mask) can branch on it. Compares the
% surface_space tag and the brain_model layout (per-model structure, count, and
% surface vertex count), then the in-data selection (vertlist / voxlist).
%
% :Outputs:
%   **isdiff:** integer code:
%       - 0  same space and same in-data grayordinates
%       - 1  different spaces (tag or model layout differ)
%       - 2  no brain_model for one or more objects
%       - 3  same space, but different in-data grayordinates (vertlist/voxlist
%            or row count differ)
%
% :See also: image_vector.compare_space, resample_space, cat, apply_mask

if isempty(obj.brain_model) || isempty(obj2.brain_model)
    isdiff = 2;
    return
end

% Space tag
if ~strcmp(char(obj.surface_space), char(obj2.surface_space))
    isdiff = 1;
    return
end

m1 = obj.brain_model.models;
m2 = obj2.brain_model.models;
if numel(m1) ~= numel(m2)
    isdiff = 1;
    return
end

% Model layout: structure name, in-data count, and surface vertex count
for i = 1:numel(m1)
    same = strcmp(m1{i}.struct, m2{i}.struct) ...
        && isequaln(m1{i}.numvert, m2{i}.numvert) ...
        && strcmp(m1{i}.type, m2{i}.type);
    if ~same
        isdiff = 1;
        return
    end
end

% Same layout: now check the exact in-data selection
isdiff = 0;
if size(obj.dat, 1) ~= size(obj2.dat, 1)
    isdiff = 3;
    return
end
for i = 1:numel(m1)
    if m1{i}.count ~= m2{i}.count ...
            || ~isequal(m1{i}.vertlist, m2{i}.vertlist) ...
            || ~isequal(m1{i}.voxlist, m2{i}.voxlist)
        isdiff = 3;
        return
    end
end
end
