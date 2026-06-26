function obj = cat(obj, varargin)
% cat Concatenate fmri_surface_data objects along the image (map) dimension.
%
% :Usage:
% ::
%     combined = cat(obj1, obj2, obj3, ...)
%
% Horizontally concatenates the [grayordinates x maps] data of two or more
% fmri_surface_data objects (e.g. combining subjects/contrasts into one dataset),
% along with their per-map annotations (X, Y, covariates, image_names,
% metadata_table). All objects must share the same grayordinate space
% (compare_space == 0); otherwise an error is raised. Surface analogue of
% fmri_data.cat. No replace_empty step is needed (.dat is always full, D5b).
%
% :Inputs:
%   **obj, varargin:** two or more fmri_surface_data objects on the same space.
%
% :Outputs:
%   **obj:** a single fmri_surface_data with all maps concatenated.
%
% :Examples:
% ::
%     all_subs = cat(sub1, sub2, sub3);
%     t = ttest(all_subs);
%
% :See also: horzcat, compare_space, fmri_surface_data

for i = 1:numel(varargin)
    o2 = varargin{i};
    if ~isa(o2, 'fmri_surface_data')
        error('fmri_surface_data:cat:type', 'All inputs must be fmri_surface_data objects.');
    end
    if compare_space(obj, o2) ~= 0
        error('fmri_surface_data:cat:space', ...
            ['Objects are not on the same grayordinate space (compare_space ~= 0). ' ...
             'Resample to a common space before concatenating.']);
    end

    obj.dat = [obj.dat, o2.dat];

    obj.X                = local_vcat(obj.X, o2.X);
    obj.Y                = local_vcat(obj.Y, o2.Y);
    obj.covariates       = local_vcat(obj.covariates, o2.covariates);
    obj.image_names      = local_catcell(obj.image_names, o2.image_names);
    obj.metadata_table   = local_cattable(obj.metadata_table, o2.metadata_table);
    if ~isempty(obj.images_per_session) || ~isempty(o2.images_per_session)
        obj.images_per_session = [obj.images_per_session(:); o2.images_per_session(:)]';
    end
end

obj.removed_images = false(size(obj.dat, 2), 1);
obj.removed_voxels = false(size(obj.dat, 1), 1);
obj.history{end+1} = sprintf('cat: concatenated to %d maps', size(obj.dat, 2));
end


% -------------------------------------------------------------------------
function v = local_vcat(a, b)
if isempty(a) && isempty(b), v = []; return; end
v = [a; b];          % per-map (rows = maps) fields
end

function c = local_catcell(a, b)
if isempty(a) && isempty(b), c = {}; return; end
a = reshape(a, [], 1); b = reshape(b, [], 1);
c = [a; b];
end

function t = local_cattable(a, b)
if isempty(a) && isempty(b), t = []; return; end
if isempty(a), t = b; return; end
if isempty(b), t = a; return; end
t = [a; b];
end
