function out = resample_surface(obj, target_space, varargin)
% resample_surface Resample cortical-surface data to another surface space.
%
% :Usage:
% ::
%     out = resample_surface(obj, target_space)
%     out = resample_surface(obj, target_space, 'interp', 'nearest')
%     resample_surface(obj, 'list')      % print the available target spaces
%
% Resamples the CORTICAL grayordinates of an fmri_surface_data object from its
% current surface space to a target surface space, natively (no FreeSurfer /
% Connectome Workbench). Supported spaces and the mapping between them:
%
%   fsaverage_164k  <->  fs_LR_32k         (HCP CIFTI)   -- spherical interpolation
%   fsaverage_164k  <->  fsaverage6/5/4                  -- nested icosahedral
%   fs_LR_32k       <->  fsaverage6/5/4                  -- via fsaverage frame
%
% fsaverage<->fs_LR uses the vendored HCP/Connectome-Workbench registration
% ("deformed") sphere that expresses fsaverage vertices in the fs_LR-aligned
% frame, so both meshes live in a common spherical frame and are resampled with
% barycentric (or nearest-neighbor) interpolation on the sphere. fsaverage down-
% sampling (164k -> fsaverage6/5/4) is exact: the lower-resolution FreeSurfer
% icosahedra are nested subsets (the first N vertices) of fsaverage-164k.
%
% Interpolation is barycentric ("linear") by default -- more accurate than
% nearest-neighbor and, because the interpolation weights depend only on the mesh
% geometry, they are computed ONCE and applied to all maps (a single sparse
% matrix multiply), so multi-map objects resample quickly. Binary masks and label
% (.dlabel) images default to NEAREST to preserve their discrete values.
%
% Only the cortical surface models are resampled. Any subcortical (voxel) models
% are volumetric and space-independent, so they are carried through unchanged.
%
% :Inputs:
%   **obj:**          an fmri_surface_data object whose surface_space is one of
%                     fsaverage_164k, fsaverage6, fsaverage5, fsaverage4, fs_LR_32k.
%   **target_space:** a target-space keyword (see 'list'); case-insensitive, with
%                     aliases (e.g. 'fsaverage','fsLR','hcp','fsaverage6', ...).
%
% :Optional Inputs:
%   **'interp':** 'linear' (barycentric, default for continuous data) or
%                 'nearest' (default for binary masks / label images).
%
% :Outputs:
%   **out:** an fmri_surface_data object in the target surface_space (cortex
%            resampled; subcortex, if any, unchanged).
%
% :Examples:
% ::
%     s  = vol2surf(ttest(load_image_set('emotionreg')));   % fsaverage_164k
%     s32 = resample_surface(s, 'fsLR_32k');                % -> HCP fs_LR-32k
%     s6  = resample_surface(s, 'fsaverage6');              % -> nested fsaverage6
%
%     m  = fmri_surface_data(which('S1200.MyelinMap_MSMAll.32k_fs_LR.dscalar.nii'));
%     mfs = resample_surface(m, 'fsaverage_164k');          % fs_LR-32k -> fsaverage
%
% :References:
%   fs_LR<->fsaverage registration spheres: Van Essen DC, Glasser MF, Dierker DL,
%   Harwell J, Coalson T (2012). Parcellations and hemispheric asymmetries of
%   human cerebral cortex analyzed on surface-based atlases. Cerebral Cortex
%   22(10):2241-2262 (HCP standard_mesh_atlases / resample_fsaverage).
%
% :See also: vol2surf, surf2vol, to_fmri_data, fmri_surface_data, apply_parcellation

% ---- 'list' : print available target spaces and return ----
if (ischar(target_space) || isstring(target_space)) && any(strcmpi(target_space, {'list', 'spaces', 'help'}))
    out = local_print_spaces();
    return
end

% ---- Parse options ----
interp_opt = '';                                   % '' = auto
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'interp', interp_opt = lower(varargin{i + 1});
        otherwise, error('resample_surface:badopt', 'Unknown option: %s', varargin{i});
    end
end

% ---- A patch handle or mesh struct target: resolve its space by vertex count.
% This lets rendering code resample onto an arbitrary isosurface (e.g. an addbrain
% patch) by passing the patch itself, not a keyword.
if ~(ischar(target_space) || isstring(target_space))
    target_space = local_space_from_patch(target_space);
end

% ---- Resolve source / target spaces ----
[tgt_name, tgt_N] = local_resolve_space(target_space);
[src_name, src_N] = local_resolve_space(obj.surface_space);

if strcmp(src_name, tgt_name)
    warning('resample_surface:noop', 'Object is already in %s; returning a copy.', tgt_name);
    out = obj;
    return
end

% ---- Choose interpolation method ----
if isempty(interp_opt)
    if local_is_discrete(obj), method = 'nearest'; else, method = 'linear'; end
else
    if ~any(strcmp(interp_opt, {'linear', 'barycentric', 'nearest'}))
        error('resample_surface:interp', 'interp must be ''linear'' or ''nearest''.');
    end
    if strcmp(interp_opt, 'barycentric'), interp_opt = 'linear'; end
    method = interp_opt;
end

% ---- Ensure the interpolation engine is on the path ----
p_anchor = which('fsavg_sphere_lh.mat');
if isempty(which('spherical_icosahedral_interpolation'))
    if ~isempty(p_anchor), addpath(fullfile(fileparts(p_anchor), 'src')); end
    if isempty(which('spherical_icosahedral_interpolation'))
        error('resample_surface:engine', ['Cannot find spherical_icosahedral_interpolation ' ...
            '(bundled under canlab_canonical_brains/Canonical_brains_surfaces/src).']);
    end
end
% onavg spheres live in an 'onavg' subfolder that may not be on the path yet.
if isempty(which('onavg_sphere_fsLR_lh_41k.mat')) && ~isempty(p_anchor)
    addpath(fullfile(fileparts(p_anchor), 'onavg'));
end

% ---- Dense cortical hemisphere data (medial wall -> 0 for interpolation) ----
r = reconstruct_image(obj);
nMaps = size(obj.dat, 2);
Ld = zeros(src_N, nMaps); Rd = zeros(src_N, nMaps);
if isfield(r, 'cortex_left')  && ~isempty(r.cortex_left),  Ld = double(r.cortex_left);  end
if isfield(r, 'cortex_right') && ~isempty(r.cortex_right), Rd = double(r.cortex_right); end
if size(Ld, 1) ~= src_N || size(Rd, 1) ~= src_N
    error('resample_surface:srcverts', ...
        'Expected %d vertices/hemisphere for %s; got L=%d R=%d.', src_N, src_name, size(Ld,1), size(Rd,1));
end
Ld(~isfinite(Ld)) = 0; Rd(~isfinite(Rd)) = 0;

% ---- Resample each hemisphere ----
src_fam = local_family(src_name); tgt_fam = local_family(tgt_name);
if strcmp(src_fam, 'fsaverage') && strcmp(tgt_fam, 'fsaverage') && tgt_N < src_N
    % Exact nested icosahedral downsample: lower-res vertices are the first N.
    Lt = Ld(1:tgt_N, :); Rt = Rd(1:tgt_N, :);
else
    [VsL, FsL, VsR, FsR] = local_sphere(src_name);
    [VqL, ~,  VqR, ~]    = local_sphere(tgt_name);
    WL = local_weights(VsL, FsL, VqL, method);
    WR = local_weights(VsR, FsR, VqR, method);
    Lt = WL * Ld; Rt = WR * Rd;
end

% ---- Reassemble the target object (cortex + unchanged subcortex) ----
out = local_build(obj, single(Lt), single(Rt), tgt_name, tgt_N, method);
end


% =========================================================================
% Spherical interpolation weights (computed once from geometry, reused for all
% maps as a sparse matrix W: target_verts x source_verts, so Cq = W * C).
function W = local_weights(Vs, Fs, Vq, method)
Ns = size(Vs, 1); Nt = size(Vq, 1);
probe = zeros(Ns, 1);
if strcmpi(method, 'nearest')
    [~, ~, vec] = spherical_icosahedral_interpolation(double(Vs), double(Fs), probe, double(Vq), 'nearest');
    W = sparse((1:Nt)', double(vec(:)), 1, Nt, Ns);
else
    [~, coord, vec] = spherical_icosahedral_interpolation(double(Vs), double(Fs), probe, double(Vq));
    rows = repmat((1:Nt)', 1, 3);
    W = sparse(rows(:), double(vec(:)), coord(:), Nt, Ns);
end
end


% =========================================================================
function out = local_build(obj, Lt, Rt, tgt_name, tgt_N, method)
% Rebuild an fmri_surface_data in the target space: dense cortex + any subcortex.
mL = struct('struct', 'CORTEX_LEFT',  'type', 'surf', 'start', 1, 'count', tgt_N, ...
    'numvert', tgt_N, 'vertlist', 0:tgt_N - 1, 'voxlist', []);
mR = struct('struct', 'CORTEX_RIGHT', 'type', 'surf', 'start', tgt_N + 1, 'count', tgt_N, ...
    'numvert', tgt_N, 'vertlist', 0:tgt_N - 1, 'voxlist', []);
models = {mL, mR};
new_dat = [Lt; Rt];
start = 2 * tgt_N + 1;

% Carry subcortical (voxel) models + their data through unchanged.
has_sub = false;
for i = 1:numel(obj.brain_model.models)
    m = obj.brain_model.models{i};
    if ~strcmp(m.type, 'vox'), continue; end
    has_sub = true;
    rows = m.start:(m.start + m.count - 1);
    new_dat = [new_dat; obj.dat(rows, :)]; %#ok<AGROW>
    m.start = start; start = start + m.count;
    models{end + 1} = m; %#ok<AGROW>
end

bm = struct('type', 'dense', 'length', size(new_dat, 1), 'models', {models}, ...
    'vol', obj.brain_model.vol);
bm.cluster = [];
if has_sub
    bm.grayordinate_type = sprintf('grayordinate_%d', size(new_dat, 1));
else
    bm.grayordinate_type = 'cortex_only';
end

out = fmri_surface_data('dat', new_dat, 'brain_model', bm, ...
    'surface_space', tgt_name, 'imagetype', obj.imagetype);
if ~isempty(obj.image_names), out.image_names = obj.image_names; end
if ~isempty(obj.label_table),  out.label_table  = obj.label_table;  end
out.history = obj.history;
out.history{end + 1} = sprintf('resample_surface: %s -> %s (%s interp)', ...
    obj.surface_space, tgt_name, method);
end


% =========================================================================
% Per-hemisphere sphere vertices (and faces) in the COMMON fs_LR-aligned frame.
function [VL, FL, VR, FR] = local_sphere(space)
switch space
    case 'fs_LR_32k'
        L = local_load_sphere('S1200.L.sphere.32k_fs_LR.mat');
        R = local_load_sphere('S1200.R.sphere.32k_fs_LR.mat');
        VL = L.vertices; FL = L.faces; VR = R.vertices; FR = R.faces;
    case {'fsaverage_164k', 'fsaverage6', 'fsaverage5', 'fsaverage4'}
        % fsaverage vertices expressed in the fs_LR-aligned frame (deformed sphere).
        L = local_load_sphere('fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.mat');
        R = local_load_sphere('fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.mat');
        n = local_space_n(space);
        VL = L.vertices(1:n, :); VR = R.vertices(1:n, :);    % nested subset
        if n == size(L.vertices, 1)
            FL = L.faces; FR = R.faces;
        else
            FL = convhulln(double(VL)); FR = convhulln(double(VR));   % ic6/5/4 faces
        end
    case {'onavg_41k', 'onavg_10k'}
        % onavg vertices in the fs_LR-aligned frame (TemplateFlow tpl-onavg
        % space-fsLR registration spheres), so they resample in the common frame.
        d = local_onavg_density(space);
        L = local_load_sphere(sprintf('onavg_sphere_fsLR_lh_%s.mat', d));
        R = local_load_sphere(sprintf('onavg_sphere_fsLR_rh_%s.mat', d));
        VL = L.vertices; FL = L.faces; VR = R.vertices; FR = R.faces;
    otherwise
        error('resample_surface:space', 'No sphere geometry for space %s.', space);
end
end


function d = local_onavg_density(space)
switch space
    case 'onavg_41k', d = '41k';
    case 'onavg_10k', d = '10k';
    otherwise, error('resample_surface:onavg', 'Unknown onavg density %s.', space);
end
end


function S = local_load_sphere(fname)
p = which(fname);
if isempty(p)
    error('resample_surface:sphere', ['Cannot find %s on the path (bundled under ' ...
        'canlab_canonical_brains/Canonical_brains_surfaces).'], fname);
end
S = load(p);
S.vertices = double(S.vertices); S.faces = double(S.faces);
end


% =========================================================================
function [name, n] = local_resolve_space(kw)
kw = lower(strrep(strrep(char(kw), '-', '_'), ' ', '_'));
switch kw
    case {'fslr_32k', 'fs_lr_32k', 'fslr32k', 'fslr', 'hcp', '32k', 'grayordinate', '91k', 'fs_lr'}
        name = 'fs_LR_32k';
    case {'fsaverage', 'fsaverage_164k', 'fsaverage164k', 'fsavg', 'fsaverage7', '164k'}
        name = 'fsaverage_164k';
    case {'fsaverage6', 'fsavg6', 'ico6'}, name = 'fsaverage6';
    case {'fsaverage5', 'fsavg5', 'ico5'}, name = 'fsaverage5';
    case {'fsaverage4', 'fsavg4', 'ico4'}, name = 'fsaverage4';
    case {'onavg', 'onavg_41k', 'onavg41k', 'onavg_ico6'}, name = 'onavg_41k';
    case {'onavg_10k', 'onavg10k', 'onavg_ico5'}, name = 'onavg_10k';
    otherwise
        error('resample_surface:unknownspace', ...
            ['Unknown surface space "%s". Run resample_surface(obj, ''list'') for the ' ...
             'available spaces.'], kw);
end
n = local_space_n(name);
end


function n = local_space_n(name)
switch name
    case 'fs_LR_32k',      n = 32492;
    case 'fsaverage_164k', n = 163842;
    case 'fsaverage6',     n = 40962;
    case 'fsaverage5',     n = 10242;
    case 'fsaverage4',     n = 2562;
    case 'onavg_41k',      n = 40962;
    case 'onavg_10k',      n = 10242;
    otherwise, error('resample_surface:space', 'Unknown space %s.', name);
end
end


function name = local_space_from_patch(target)
% Resolve a surface space from an isosurface patch handle or a mesh struct, by
% its per-hemisphere vertex count (a single patch is one hemisphere).
if isstruct(target) && isfield(target, 'vertices')
    nv = size(target.vertices, 1);
elseif all(ishandle(target(:))) && strcmp(get(target(1), 'Type'), 'patch')
    nv = size(get(target(1), 'Vertices'), 1);
else
    error('resample_surface:target', ...
        ['Target must be a surface-space keyword, a patch handle, or a struct with ' ...
         '.vertices. Run resample_surface(obj, ''list'').']);
end
% Standard cortical meshes are recognized by vertex count. Ambiguous counts
% (40962/10242/2562 could be fsaverage or onavg) resolve to the FreeSurfer
% fsaverage mesh, which is what addbrain surfaces use.
switch nv
    case 32492,  name = 'fs_LR_32k';
    case 163842, name = 'fsaverage_164k';
    case 40962,  name = 'fsaverage6';
    case 10242,  name = 'fsaverage5';
    case 2562,   name = 'fsaverage4';
    otherwise
        error('resample_surface:unknownmesh', ...
            ['The target surface has %d vertices/hemisphere, which is not a recognized ' ...
             'standard cortical mesh, so the data cannot be resampled onto it. Render ' ...
             'onto a standard fsaverage or fs_LR surface, or project via a volume ' ...
             '(surf2vol + render_on_surface).'], nv);
end
end


function fam = local_family(name)
% Families gate the exact nested-subset fast path (fsaverage-only). onavg is a
% different (uniform-area) tessellation, so it always uses interpolation.
if strcmp(name, 'fs_LR_32k'), fam = 'fsLR';
elseif strncmp(name, 'onavg', 5), fam = 'onavg';
else, fam = 'fsaverage'; end
end


% =========================================================================
function tf = local_is_discrete(obj)
% Binary masks and label/atlas images resample by nearest-neighbor.
tf = false;
if ~isempty(obj.imagetype) && any(strcmpi(obj.imagetype, {'dlabel', 'label'}))
    tf = true; return
end
v = obj.dat(:); v = double(v(isfinite(v)));
u = unique(v);
if numel(u) <= 2 && all(ismember(u, [0 1])), tf = true; end     % binary mask
end


% =========================================================================
function t = local_print_spaces()
kw   = {'fsaverage_164k'; 'fsaverage6'; 'fsaverage5'; 'fsaverage4'; 'fs_LR_32k'; 'onavg_41k'; 'onavg_10k'};
nver = {163842; 40962; 10242; 2562; 32492; 40962; 10242};
desc = { ...
    'FreeSurfer fsaverage, 7th-order icosahedron (vol2surf output)'; ...
    'FreeSurfer fsaverage6 (6th-order icosahedron; nested subset of 164k)'; ...
    'FreeSurfer fsaverage5 (5th-order icosahedron; nested subset)'; ...
    'FreeSurfer fsaverage4 (4th-order icosahedron; nested subset)'; ...
    'HCP fs_LR-32k (native CIFTI grayordinate cortex; Van Essen 2012)'; ...
    'onavg equal-area template, den-41k (Feilong et al. 2024; CC0)'; ...
    'onavg equal-area template, den-10k'};
alias = { ...
    'fsaverage, fsavg, fsaverage7, 164k'; ...
    'fsavg6, ico6'; 'fsavg5, ico5'; 'fsavg4, ico4'; ...
    'fsLR, fs_LR, hcp, 32k, 91k'; ...
    'onavg, onavg41k'; 'onavg10k'};
t = table(kw, nver, desc, alias, 'VariableNames', {'space', 'verts_per_hemi', 'description', 'aliases'});
fprintf('\n==== resample_surface: available target spaces ====\n');
disp(t);
fprintf(['Usage: out = resample_surface(obj, ''<space>'');  add ''interp'',''nearest'' for ' ...
    'label/binary maps.\n\n']);
end
