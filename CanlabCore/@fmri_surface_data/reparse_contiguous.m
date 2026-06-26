function [obj, ncl] = reparse_contiguous(obj, varargin)
% reparse_contiguous Label contiguous clusters of grayordinates (mesh + volume).
%
% :Usage:
% ::
%     [obj, ncl] = reparse_contiguous(obj)
%     [obj, ncl] = reparse_contiguous(obj, 'which_image', 2)
%
% Finds connected components among the "active" (nonzero, non-NaN) grayordinates
% and stores an integer cluster label per grayordinate in obj.brain_model.cluster
% (0 = not in any cluster). For cortical-surface models, contiguity is computed
% on the cortical mesh edge graph (using the bundled mesh for the object's
% surface_space and MATLAB's built-in graph/conncomp -- no external toolbox). For
% subcortical voxel models, contiguity is 26-connected in the volume (spm_clusters).
%
% This is the surface analogue of image_vector.reparse_contiguous and is the
% basis for the region() conversion (a later milestone).
%
% :Optional Inputs:
%   **'which_image':** which map (column) defines "active" grayordinates. Default 1.
%
% :Outputs:
%   **obj:** the object with obj.brain_model.cluster populated [nGrayordinates x 1].
%   **ncl:** total number of clusters found.
%
% :See also: fmri_surface_data, load_surface_geom, region, reconstruct_image

which_image = 1;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'which_image', which_image = varargin{i+1};
        otherwise, error('reparse_contiguous:badopt', 'Unknown option: %s', varargin{i});
    end
end

d = double(obj.dat(:, which_image));
active = d ~= 0 & ~isnan(d);
cluster = zeros(size(obj.dat, 1), 1);
offset = 0;

models = obj.brain_model.models;

for mi = 1:numel(models)
    m = models{mi};
    rows = (m.start:(m.start + m.count - 1))';
    act_local = active(rows);
    if ~any(act_local), continue; end

    if strcmp(m.type, 'surf')
        F = local_hemi_faces(obj.surface_space, m.struct);
        G = local_mesh_graph(F, m.numvert);
        meshverts = m.vertlist(:) + 1;                 % 1-based, aligned with rows
        active_nodes = meshverts(act_local);
        comp = conncomp(subgraph(G, active_nodes));    % component id per active node
        cluster(rows(act_local)) = comp(:) + offset;
        offset = offset + max(comp);

    else   % voxel model: 26-connected components in the volume
        ijk = m.voxlist(:, act_local);                 % 3 x nactive (0-based)
        if isempty(ijk), continue; end
        a = spm_clusters(ijk + 1);                     % 1-based voxel coords
        cluster(rows(act_local)) = a(:) + offset;
        offset = offset + max(a);
    end
end

obj.brain_model.cluster = cluster;
ncl = offset;
obj.history{end+1} = sprintf('reparse_contiguous: %d clusters (mesh + volume, image %d)', ncl, which_image);
end


% -------------------------------------------------------------------------
function F = local_hemi_faces(surface_space, structname)
geom = load_surface_geom(surface_space, 'inflated');   % faces shared across surftypes
if contains(upper(structname), 'LEFT')
    F = geom.faces_lh;
else
    F = geom.faces_rh;
end
end


% -------------------------------------------------------------------------
function G = local_mesh_graph(F, numvert)
% Build an undirected graph over numvert mesh vertices from triangle faces.
E = [F(:, [1 2]); F(:, [2 3]); F(:, [1 3])];
E = sort(E, 2);
E = unique(E, 'rows');
G = graph(E(:, 1), E(:, 2), [], numvert);
end
