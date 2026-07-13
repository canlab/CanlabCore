function out_handles = render_on_surface(obj, surface_handles, varargin)
% render_on_surface Color existing surface patches with grayordinate data.
%
% :Usage:
% ::
%     render_on_surface(surf_obj, surface_handles)
%     render_on_surface(surf_obj, surface_handles, 'colormap', 'hot', 'cmaprange', [0 5])
%     render_on_surface(surf_obj, surface_handles, 'maxcolor', [1 1 0], 'mincolor', [1 0 0])
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
%     image_vector.render_on_surface is reused to sample the volume.
%
% Color options are harmonized with the volume visualization pipeline
% (render_on_surface / addblobs / set_colormap), so the same keywords color
% surface data and volume blobs.
%
% :Inputs:
%   **obj:**             an fmri_surface_data object.
%   **surface_handles:** vector of patch handles (e.g. from addbrain or surface).
%
% :Optional Inputs (color, matching the volume pipeline):
%   **'which_image':**  map (column) to render. Default 1.
%   **'clim' / 'cmaprange':** [lo hi] color limits. Default symmetric from data.
%   **'colormap' / 'colormapname':** a named MATLAB colormap ('hot','parula',...)
%        or an [n x 3] matrix -> a single (sequential) colormap over clim.
%   **'pos_colormap' / 'neg_colormap':** [n x 3] split colormaps (default hot / cool).
%   **'splitcolor':**   {neg_low, neg_high, pos_low, pos_high} colors -> split map.
%   **'maxcolor' / 'mincolor':** endpoint colors -> a single gradient over clim.
%   **'color' / 'solid':** a single solid color for all in-data vertices.
%
% :Outputs:
%   **out_handles:** the surface handles (colored in place).
%
% :See also: surface, canlab_surface_vertexcolors, image_vector.render_on_surface

which_image = 1;
clim = [];
poscm = [];
negcm = [];
cmap_single = [];        % single sequential colormap over clim
splitcolors = [];        % {neg_low neg_high pos_low pos_high}
maxcolor = [];
mincolor = [];
solidcolor = [];         % single solid color for in-data vertices
coloropts = {};          % color options to forward to the volume renderer
i = 1;
passthrough = {};
while i <= numel(varargin)
    key = varargin{i};
    if ischar(key) || isstring(key)
        switch lower(char(key))
            case 'which_image', which_image = varargin{i+1}; i = i + 2; continue
            case {'clim','cmaprange'}, clim = varargin{i+1}; coloropts = [coloropts {'clim', clim}]; i = i + 2; continue %#ok<AGROW>
            case 'pos_colormap', poscm = varargin{i+1}; coloropts = [coloropts {'pos_colormap', poscm}]; i = i + 2; continue %#ok<AGROW>
            case 'neg_colormap', negcm = varargin{i+1}; coloropts = [coloropts {'neg_colormap', negcm}]; i = i + 2; continue %#ok<AGROW>
            case {'colormap','colormapname'}, cmap_single = local_resolve_cmap(varargin{i+1}); coloropts = [coloropts {'colormap', varargin{i+1}}]; i = i + 2; continue %#ok<AGROW>
            case 'splitcolor', splitcolors = varargin{i+1}; coloropts = [coloropts {'splitcolor', splitcolors}]; i = i + 2; continue %#ok<AGROW>
            case 'maxcolor', maxcolor = varargin{i+1}; i = i + 2; continue
            case 'mincolor', mincolor = varargin{i+1}; i = i + 2; continue
            case {'color','solid','onecolor'}, solidcolor = varargin{i+1}; coloropts = [coloropts {'color', solidcolor}]; i = i + 2; continue %#ok<AGROW>
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

% Resolve the color spec from the harmonized options
spec = struct('mode', 'split', 'clim', clim, 'poscm', poscm, 'negcm', negcm, ...
    'cmap', cmap_single, 'solid', solidcolor);
if ~isempty(solidcolor)
    spec.mode = 'solid';
elseif ~isempty(maxcolor) && ~isempty(mincolor)
    spec.mode = 'single';
    spec.cmap = local_interp_cmap(mincolor, maxcolor, 256);
    coloropts = [coloropts {'colormap', spec.cmap}];
elseif ~isempty(cmap_single)
    spec.mode = 'single';
elseif ~isempty(splitcolors) && numel(splitcolors) >= 4
    spec.negcm = local_interp_cmap(splitcolors{1}, splitcolors{2}, 256);
    spec.poscm = local_interp_cmap(splitcolors{3}, splitcolors{4}, 256);
end

for h = 1:numel(surface_handles)
    hp = surface_handles(h);
    nv = size(get(hp, 'Vertices'), 1);
    tag = '';
    try, tag = lower(get(hp, 'Tag')); catch, end %#ok<NOCOM>

    isLeft  = contains(tag, 'left')  || contains(tag, ' l ') || endsWith(tag, ' l');
    isRight = contains(tag, 'right') || endsWith(tag, ' r');

    if nv == nL && (isLeft || (~isRight && ~isempty(Ldat)))
        local_apply(hp, Ldat, spec);
    elseif nv == nR && (isRight || ~isempty(Rdat))
        local_apply(hp, Rdat, spec);
    else
        % Arbitrary surface: project to a volume and reuse the image_vector
        % renderer, forwarding the harmonized color options.
        vol = obj_to_volume(obj);
        if which_image > 1 && size(vol.dat, 2) >= which_image
            vol = get_wh_image(vol, which_image);
        end
        render_on_surface(vol, surface_handles, coloropts{:}, passthrough{:});
        out_handles = surface_handles;
        return
    end
end

out_handles = surface_handles;
end


% -------------------------------------------------------------------------
function local_apply(hp, vals, spec)
rgb = local_vertex_rgb(vals, spec);
set(hp, 'FaceVertexCData', rgb, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 1);
end


function rgb = local_vertex_rgb(vals, spec)
vals = double(vals(:));
N = numel(vals);
graycolor = [.5 .5 .5];

switch spec.mode
    case 'solid'
        rgb = repmat(graycolor, N, 1);
        indata = ~isnan(vals) & vals ~= 0;
        rgb(indata, :) = repmat(spec.solid(:)', nnz(indata), 1);

    case 'single'
        cmap = spec.cmap;
        n = size(cmap, 1);
        rgb = repmat(graycolor, N, 1);
        lo = spec.clim(1); hi = spec.clim(2);
        if hi == lo, hi = lo + 1; end
        ok = ~isnan(vals);
        idx = round(1 + (vals(ok) - lo) / (hi - lo) * (n - 1));
        idx = min(max(idx, 1), n);
        rgb(ok, :) = cmap(idx, :);

    otherwise   % 'split' (default): hot for +, cool for -, gray for 0/NaN
        rgb = canlab_surface_vertexcolors(vals, spec.clim, spec.poscm, spec.negcm);
end
end


function cmap = local_resolve_cmap(c)
if isnumeric(c) && size(c, 2) == 3
    cmap = c;
elseif ischar(c) || isstring(c)
    try
        cmap = feval(char(c), 256);
    catch
        error('fmri_surface_data:render_on_surface:colormap', ...
            'Unknown colormap "%s". Use a MATLAB colormap name or an [n x 3] matrix.', char(c));
    end
else
    error('fmri_surface_data:render_on_surface:colormap', ...
        'colormap must be a name or an [n x 3] matrix.');
end
end


function cmap = local_interp_cmap(clo, chi, n)
clo = double(clo(:))';
chi = double(chi(:))';
t = linspace(0, 1, n)';
cmap = (1 - t) .* clo + t .* chi;
end
