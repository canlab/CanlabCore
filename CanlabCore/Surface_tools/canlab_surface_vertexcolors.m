function rgb = canlab_surface_vertexcolors(vals, clim, poscm, negcm, graycolor)
% canlab_surface_vertexcolors Map per-vertex values to truecolor RGB (split colormap).
%
% :Usage:
% ::
%     rgb = canlab_surface_vertexcolors(vals, clim, poscm, negcm, graycolor)
%
% Maps a vector of per-vertex values to an [N x 3] truecolor matrix using a
% split colormap: positive values use a "hot"-style map, negative values a
% "cool"-style map, and zero / NaN vertices (e.g. medial wall, out-of-mask) are
% rendered in a neutral gray. Designed for direct FaceVertexCData assignment on
% surface patches (FaceColor 'interp'), so no axis colormap/caxis juggling is
% needed and the medial wall stays gray.
%
% :Inputs:
%   **vals:** [N x 1] per-vertex values (NaN allowed for medial wall).
%
% :Optional Inputs:
%   **clim:**      [lo hi] color limits. Default: symmetric [-m m] where m is the
%                  max abs nonzero value. Mapping uses max(abs(clim)).
%   **poscm:**     [P x 3] positive colormap. Default hot(256).
%   **negcm:**     [Q x 3] negative colormap. Default cool(256).
%   **graycolor:** 1x3 color for zero/NaN vertices. Default [.5 .5 .5].
%
% :Outputs:
%   **rgb:** [N x 3] truecolor matrix in [0,1].
%
% :See also: surface, render_on_surface, fmri_surface_data

vals = double(vals(:));
N = numel(vals);

if nargin < 5 || isempty(graycolor), graycolor = [.5 .5 .5]; end
if nargin < 4 || isempty(negcm), negcm = cool(256); end
if nargin < 3 || isempty(poscm), poscm = hot(256); end
if nargin < 2 || isempty(clim)
    finite_nonzero = vals(~isnan(vals) & vals ~= 0);
    m = max(abs(finite_nonzero));
    if isempty(m) || m == 0, m = 1; end
    clim = [-m m];
end

hi = max(abs(clim));
if hi == 0, hi = 1; end

rgb = repmat(graycolor, N, 1);
pos = vals > 0 & ~isnan(vals);
neg = vals < 0 & ~isnan(vals);

idx = @(v, n) min(max(round(1 + (abs(v) / hi) * (n - 1)), 1), n);
if any(pos), rgb(pos, :) = poscm(idx(vals(pos), size(poscm, 1)), :); end
if any(neg), rgb(neg, :) = negcm(idx(vals(neg), size(negcm, 1)), :); end
end
