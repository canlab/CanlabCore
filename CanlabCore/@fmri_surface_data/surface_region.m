function reg = surface_region(obj, varargin)
% surface_region Summarize contiguous grayordinate clusters as region structs.
%
% :Usage:
% ::
%     reg = surface_region(obj)
%     reg = surface_region(obj, 'which_image', 2)
%
% Converts the contiguous clusters of an fmri_surface_data (active = nonzero,
% non-NaN grayordinates) into a struct array of per-cluster summaries -- the
% surface analogue of region(). Clusters are found with reparse_contiguous (mesh
% edge graph for cortex, 26-connectivity for subcortex). Cortical-cluster
% centroids are computed from the bundled mesh (midthickness if available);
% subcortical-cluster centroids come from voxel mm coordinates.
%
% (A full CANlab region object is volume-centric; this returns a lightweight
% surface-aware struct. For the subcortical part you can also use
% region(to_fmri_data(obj)).)
%
% :Optional Inputs:
%   **'which_image':** which map defines "active" grayordinates. Default 1.
%
% :Outputs:
%   **reg:** struct array, one element per cluster, with fields:
%       .struct        brain structure name (e.g. 'CORTEX_LEFT')
%       .type          'surf' or 'vox'
%       .cluster_id    integer cluster label
%       .grayord_rows  row indices into .dat for this cluster
%       .vertex_indices 0-based mesh vertices (surf clusters) or []
%       .XYZmm         [3 x 1] centroid in mm
%       .numVox        number of grayordinates in the cluster
%       .val           mean value (which_image) over the cluster
%
% :See also: reparse_contiguous, region, to_fmri_data, fmri_surface_data

which_image = 1;
for i = 1:2:numel(varargin)
    if strcmpi(varargin{i}, 'which_image'), which_image = varargin{i+1}; end
end

obj = reparse_contiguous(obj, 'which_image', which_image);
cluster = obj.brain_model.cluster;
d = double(obj.dat(:, which_image));

reg = struct('struct', {}, 'type', {}, 'cluster_id', {}, 'grayord_rows', {}, ...
    'vertex_indices', {}, 'XYZmm', {}, 'numVox', {}, 'val', {});

if isempty(cluster) || all(cluster == 0), return; end

% mesh vertices for centroids (cortex)
geom = [];
% subcortical affine for centroids (vox)
hasvol = isfield(obj.brain_model, 'vol') && ~isempty(obj.brain_model.vol);

for mi = 1:numel(obj.brain_model.models)
    m = obj.brain_model.models{mi};
    rows = (m.start:(m.start + m.count - 1))';
    cl_local = cluster(rows);
    ids = unique(cl_local(cl_local > 0));

    if strcmp(m.type, 'surf') && ~isempty(ids)
        if isempty(geom)
            try, geom = load_geom_safe(obj.surface_space); catch, geom = []; end
        end
        if contains(upper(m.struct), 'LEFT') && ~isempty(geom), Vh = geom.vertices_lh;
        elseif ~isempty(geom), Vh = geom.vertices_rh; else, Vh = []; end
    end

    for c = ids(:)'
        sel = find(cl_local == c);
        gr = rows(sel);
        s = struct();
        s.struct = m.struct;
        s.type = m.type;
        s.cluster_id = c;
        s.grayord_rows = gr;
        s.numVox = numel(gr);
        s.val = mean(d(gr));
        if strcmp(m.type, 'surf')
            verts = m.vertlist(sel) + 1;          % 1-based mesh verts
            s.vertex_indices = m.vertlist(sel);   % store 0-based
            if exist('Vh', 'var') && ~isempty(Vh)
                s.XYZmm = mean(Vh(verts, :), 1)';
            else
                s.XYZmm = [NaN; NaN; NaN];
            end
        else
            s.vertex_indices = [];
            ijk = m.voxlist(:, sel) + 1;          % 1-based
            if hasvol
                mm = obj.brain_model.vol.sform * [mean(ijk - 1, 2); 1];  % sform is 0-based
                s.XYZmm = mm(1:3);
            else
                s.XYZmm = [NaN; NaN; NaN];
            end
        end
        reg(end+1) = s; %#ok<AGROW>
    end
end
end


% -------------------------------------------------------------------------
function g = load_geom_safe(space)
% Use midthickness for accurate centroids; fall back to inflated.
try
    g = call_private_geom(space, 'midthickness');
catch
    g = call_private_geom(space, 'inflated');
end
end

function g = call_private_geom(space, surftype)
% load_surface_geom is private to @fmri_surface_data; reachable from class methods.
g = load_surface_geom(space, surftype);
end
