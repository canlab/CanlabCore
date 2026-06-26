function out_handles = render_on_surface(obj, surface_handles, varargin)
% render_on_surface Color existing surface patches with grayordinate data.
%
% :Usage:
% ::
%     render_on_surface(surf_obj, surface_handles)
%     render_on_surface(surf_obj, surface_handles, 'clim', [-3 3], 'which_image', 2)
%
% Colors one or more existing surface patch handles using an fmri_surface_data
% object's per-vertex values. For each patch:
%   - If the patch has the same vertex count as one of the object's cortical
%     hemispheres (e.g. the matching fs_LR / fsaverage inflated, midthickness, or
%     sphere mesh from addbrain), the per-vertex values are applied DIRECTLY
%     (no resampling) -- true native-space rendering. Left/right is resolved by
%     the patch Tag ('left'/'right') when present, else by vertex count.
%   - Otherwise (an arbitrary MNI surface, e.g. addbrain('left')), the object is
%     projected to a volume (surf2vol / to_fmri_data) and the standard
%     image_vector.render_on_surface is reused to sample the volume at the patch
%     vertices.
%
% :Inputs:
%   **obj:**             an fmri_surface_data object.
%   **surface_handles:** vector of patch handles (e.g. from addbrain or surface).
%
% :Optional Inputs:
%   **'which_image':** map (column) to render. Default 1.
%   **'clim':**        [lo hi] color limits. Default symmetric from the data.
%   **'pos_colormap' / 'neg_colormap':** [n x 3] colormaps (default hot / cool).
%   Other options are passed through to image_vector.render_on_surface on the
%   volume fallback path.
%
% :Outputs:
%   **out_handles:** the surface handles (colored in place).
%
% :See also: surface, canlab_surface_vertexcolors, image_vector.render_on_surface

which_image = 1;
clim = [];
poscm = [];
negcm = [];
i = 1;
passthrough = {};
while i <= numel(varargin)
    key = varargin{i};
    if ischar(key) || isstring(key)
        switch lower(char(key))
            case 'which_image', which_image = varargin{i+1}; i = i + 2; continue
            case 'clim',        clim = varargin{i+1};        i = i + 2; continue
            case 'pos_colormap', poscm = varargin{i+1};      i = i + 2; continue
            case 'neg_colormap', negcm = varargin{i+1};      i = i + 2; continue
        end
    end
    passthrough{end+1} = varargin{i}; %#ok<AGROW>
    i = i + 1;
end

% Dense per-hemisphere data (medial wall = NaN)
r = reconstruct_image(obj);
Ldat = []; Rdat = [];
if isfield(r, 'cortex_left'),  Ldat = r.cortex_left(:, which_image);  end
if isfield(r, 'cortex_right'), Rdat = r.cortex_right(:, which_image); end
nL = numel(Ldat); nR = numel(Rdat);

% Color limits shared across hemispheres
if isempty(clim)
    allv = [Ldat; Rdat];
    allv = allv(~isnan(allv) & allv ~= 0);
    m = max(abs(allv)); if isempty(m) || m == 0, m = 1; end
    clim = [-m m];
end

for h = 1:numel(surface_handles)
    hp = surface_handles(h);
    nv = size(get(hp, 'Vertices'), 1);
    tag = '';
    try, tag = lower(get(hp, 'Tag')); catch, end %#ok<NOCOM>

    isLeft  = contains(tag, 'left')  || contains(tag, ' l ') || endsWith(tag, ' l');
    isRight = contains(tag, 'right') || endsWith(tag, ' r');

    if nv == nL && (isLeft || (~isRight && ~isempty(Ldat)))
        local_apply(hp, Ldat, clim, poscm, negcm);
    elseif nv == nR && (isRight || ~isempty(Rdat))
        local_apply(hp, Rdat, clim, poscm, negcm);
    else
        % Arbitrary surface: project to volume and reuse image_vector renderer
        vol = obj_to_volume(obj);
        render_on_surface(vol, surface_handles, passthrough{:});
        out_handles = surface_handles;
        return
    end
end

out_handles = surface_handles;
end


% -------------------------------------------------------------------------
function local_apply(hp, vals, clim, poscm, negcm)
rgb = canlab_surface_vertexcolors(vals, clim, poscm, negcm);
set(hp, 'FaceVertexCData', rgb, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 1);
end
