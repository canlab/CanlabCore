function geom = load_surface_geom(surface_space, surftype)
% load_surface_geom Load bundled L/R cortical meshes for a surface space.
%
% Returns faces/vertices for both hemispheres of a standard surface, loaded from
% the meshes bundled in CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/.
% No gifti toolbox is required (.mat meshes are load()ed; .surf.gii is read with
% canlab_read_gifti).
%
% :Inputs:
%   **surface_space:** 'fsLR_32k' (32492 verts/hemi) or 'fsaverage_164k' (163842).
%   **surftype:**      'inflated' (default), 'midthickness', 'sphere', 'veryinflated'.
%
% :Outputs:
%   **geom:** struct with .vertices_lh/.faces_lh/.vertices_rh/.faces_rh (faces
%             1-based), .space, .surftype, .numvert.
%
% :See also: surface, fmri_surface_data, add_surface

if nargin < 2 || isempty(surftype), surftype = 'inflated'; end

switch surface_space
    case 'fsLR_32k'
        switch lower(surftype)
            case {'inflated','veryinflated'}
                fL = 'S12000.L.inflated_MSMAll.32k_fsl_LR.mat';
                fR = 'S12000.R.inflated_MSMAll.32k_fsl_LR.mat';
            case 'midthickness'
                fL = 'S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii';
                fR = 'S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii';
            case 'sphere'
                fL = 'S1200.L.sphere.32k_fs_LR.mat';
                fR = 'S1200.R.sphere.32k_fs_LR.mat';
            otherwise
                error('load_surface_geom:surftype', 'Unknown fsLR_32k surftype: %s', surftype);
        end
        nv = 32492;

    case 'fsaverage_164k'
        switch lower(surftype)
            case {'inflated','veryinflated'}
                fL = 'surf_freesurf_inflated_Left.mat';
                fR = 'surf_freesurf_inflated_Right.mat';
            otherwise
                error('load_surface_geom:surftype', ...
                    'fsaverage_164k currently supports surftype ''inflated'' only (got %s).', surftype);
        end
        nv = 163842;

    otherwise
        error('load_surface_geom:space', ...
            ['No bundled mesh for surface_space ''%s''. Supported: fsLR_32k, ' ...
             'fsaverage_164k. (vol2surf produces fsaverage_164k; native CIFTI is fsLR_32k.)'], ...
            surface_space);
end

[vL, faL] = local_load_mesh(fL);
[vR, faR] = local_load_mesh(fR);

geom = struct('vertices_lh', vL, 'faces_lh', faL, ...
              'vertices_rh', vR, 'faces_rh', faR, ...
              'space', surface_space, 'surftype', surftype, 'numvert', nv);
end


% -------------------------------------------------------------------------
function [vertices, faces] = local_load_mesh(fname)
p = which(fname);
if isempty(p), error('load_surface_geom:notfound', 'Mesh file not found on path: %s', fname); end
if endsWith(lower(p), '.gii')
    g = canlab_read_gifti(p);
    vertices = g.vertices;
    faces = g.faces;                 % already 1-based
else
    S = load(p);
    vertices = S.vertices;
    faces = S.faces;                 % .mat meshes store 1-based faces
end
end
