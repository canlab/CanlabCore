function [parcel_means, parcel_labels, parcel_table] = apply_parcellation(obj, parcels, varargin)
% apply_parcellation Average grayordinate data within parcels of a surface atlas.
%
% :Usage:
% ::
%     parcel_means = apply_parcellation(obj, parcels)
%     [parcel_means, labels, tbl] = apply_parcellation(obj, parcels, 'area')
%
% Computes the mean of each map over the grayordinates in each parcel of a
% surface/grayordinate parcellation, mirroring image_vector.apply_parcellation.
% Parcels with key 0 (and NaN) are treated as background / medial wall and are
% excluded. Operates directly on .dat (no resampling needed when on the same
% space).
%
% :Inputs:
%   **obj:**     an fmri_surface_data object ([nGrayordinates x nMaps]).
%   **parcels:** the parcellation, one of
%       - an fmri_surface_data with integer label keys (e.g. a `.dlabel`),
%         on the SAME grayordinate space (compare_space == 0); its label_table
%         provides parcel names.
%       - an integer vector [nGrayordinates x 1] of label keys.
%
% :Optional Inputs:
%   **'area':** area-weight the average using per-vertex surface area (cortex)
%               instead of an unweighted mean. (Subcortical voxels use unit weight.)
%
% :Outputs:
%   **parcel_means:** [nMaps x nParcels] mean value per parcel.
%   **parcel_labels:** 1 x nParcels cell of parcel names (or 'parcel_<key>').
%   **parcel_table:**  table with columns key, label, n_grayordinates, (area).
%
% :Examples:
% ::
%     atl = fmri_surface_data(which('Gordon333.32k_fs_LR_Tian_Subcortex_S2.dlabel.nii'));
%     s   = fmri_surface_data(which('transcriptomic_gradients.dscalar.nii'));
%     pm  = apply_parcellation(s, atl);          % [nMaps x nParcels]
%
% :See also: fmri_surface_data, reparse_contiguous, surface_region, condf2indic

use_area = any(strcmpi(varargin, 'area'));

% ---- Resolve parcel keys + names ----
names_by_key = containers.Map('KeyType', 'double', 'ValueType', 'char');
if isa(parcels, 'fmri_surface_data')
    if compare_space(obj, parcels) ~= 0
        error('fmri_surface_data:apply_parcellation:space', ...
            'Parcellation is not on the same grayordinate space (compare_space ~= 0).');
    end
    keys = round(double(parcels.dat(:, 1)));
    if ~isempty(parcels.label_table)
        for i = 1:numel(parcels.label_table)
            names_by_key(parcels.label_table(i).key) = parcels.label_table(i).name;
        end
    end
elseif isnumeric(parcels)
    keys = round(double(parcels(:)));
    if numel(keys) ~= size(obj.dat, 1)
        error('fmri_surface_data:apply_parcellation:length', ...
            'Parcel key vector length (%d) must equal the number of grayordinates (%d).', ...
            numel(keys), size(obj.dat,1));
    end
else
    error('fmri_surface_data:apply_parcellation:type', ...
        'parcels must be an fmri_surface_data or an integer vector.');
end

% ---- Indicator matrix over positive keys (exclude 0 / NaN) ----
ukeys = unique(keys(keys > 0 & ~isnan(keys)));
nP = numel(ukeys);
if nP == 0
    error('fmri_surface_data:apply_parcellation:noparcels', 'No positive parcel keys found.');
end
indic = double(keys == ukeys');                 % [nGray x nP] 0/1

% ---- Per-grayordinate weights ----
if use_area
    w = local_vertex_areas(obj);                % [nGray x 1]
else
    w = ones(size(obj.dat, 1), 1);
end

D = double(obj.dat);                            % [nGray x nMaps]
wsum = (w' * indic);                            % [1 x nP] total weight per parcel
sums = (D .* w)' * indic;                       % [nMaps x nP]
parcel_means = sums ./ wsum;                    % [nMaps x nP]

% ---- Labels + table ----
parcel_labels = cell(1, nP);
counts = sum(indic, 1);                         % grayordinates per parcel
for i = 1:nP
    if names_by_key.isKey(ukeys(i)) && ~isempty(names_by_key(ukeys(i)))
        parcel_labels{i} = names_by_key(ukeys(i));
    else
        parcel_labels{i} = sprintf('parcel_%d', ukeys(i));
    end
end

if nargout >= 3
    parcel_table = table(ukeys(:), parcel_labels(:), counts(:), wsum(:), ...
        'VariableNames', {'key', 'label', 'n_grayordinates', 'total_weight'});
end
end


% -------------------------------------------------------------------------
function w = local_vertex_areas(obj)
% Per-grayordinate weight = barycentric vertex area on the cortical mesh
% (1/3 of incident face areas); subcortical voxels get unit weight.
w = ones(size(obj.dat, 1), 1);
for mi = 1:numel(obj.brain_model.models)
    m = obj.brain_model.models{mi};
    rows = (m.start:(m.start + m.count - 1))';
    if ~strcmp(m.type, 'surf'), continue; end
    [V, F] = local_hemi_VF(obj.surface_space, m.struct);
    va = local_face_vertex_area(V, F);          % [numvert x 1]
    w(rows) = va(m.vertlist + 1);
end
end


function [V, F] = local_hemi_VF(space, structname)
% Prefer midthickness (true areas); fall back to inflated if unavailable.
try
    g = load_surface_geom(space, 'midthickness');
catch
    g = load_surface_geom(space, 'inflated');
end
if contains(upper(structname), 'LEFT'), V = g.vertices_lh; F = g.faces_lh;
else, V = g.vertices_rh; F = g.faces_rh; end
end


function va = local_face_vertex_area(V, F)
% Barycentric per-vertex area: each vertex gets 1/3 of each incident face's area.
v1 = V(F(:,1), :); v2 = V(F(:,2), :); v3 = V(F(:,3), :);
fa = 0.5 * sqrt(sum(cross(v2 - v1, v3 - v1, 2).^2, 2));   % [nFaces x 1]
va = accumarray(F(:), repmat(fa, 3, 1) / 3, [size(V,1) 1]);
end
