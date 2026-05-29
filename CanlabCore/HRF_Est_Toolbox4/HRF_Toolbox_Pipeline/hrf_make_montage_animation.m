function out_video = hrf_make_montage_animation(img_obj_or_file, out_video, varargin)
%HRF_MAKE_MONTAGE_ANIMATION Fast montage-based 3D animation over time.
% Supports raw volume, beta maps, or t maps in fmri_data format.
%
% Color limits are held constant across frames by passing 'cmaprange',
% 'pos_colormap', 'neg_colormap' directly to canlabCore's montage --
% setting CLim on the figure axes after montage() draws does NOT override
% the per-blob color mapping that addblobs builds internally, which is why
% an earlier "set CLim post-hoc" approach left one half of the bound
% drifting frame to frame.
%
% Name-value parameters
% ---------------------
%   'FrameStep'    - integer, animate every Nth volume (default 2)
%   'FPS'          - playback rate (default 8)
%   'Threshold'    - scalar; voxels with |value| < Threshold are zeroed
%                    (via statistic_image/threshold). Default [] (no thresh).
%   'UseSigMask'   - true (default) to honor obj.sig for statistic_image input
%   'TitlePrefix'  - prefix for the per-frame title
%   'ColorLimits'  - 'auto' (default) -> symmetric robust [-q q] from the
%                    98th percentile of |non-zero, finite, sig-masked|
%                    values across every frame in the animation,
%                    OR [lo hi] to set explicitly. Held constant across
%                    every frame.
%   'ColorMap'     - 'rdylbu' (default), 'hotcool' (canlab classic), 'bwr',
%                    or a {pos_cm, neg_cm} 2-cell pair of Nx3 RGB matrices.
%   'GrayBuffer'   - default false. canlabCore's HCP surface display has
%                    a gray buffer around blob boundaries that's helpful
%                    for unthresholded volumetric maps but clutters
%                    thresholded statistical map animations.
%   'MontageType'  - default 'hcp'. Pass-through layout name to canlab
%                    montage() (e.g. 'hcp', 'compact2', 'full').

p = inputParser;
p.addRequired('img_obj_or_file');
p.addRequired('out_video', @(x) ischar(x) || isstring(x));
p.addParameter('FrameStep', 2, @(x) isscalar(x) && x >= 1);
p.addParameter('FPS', 8, @(x) isscalar(x) && x > 0);
p.addParameter('Threshold', [], @(x) isempty(x) || isscalar(x));
p.addParameter('UseSigMask', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('TitlePrefix', 'Frame', @(x) ischar(x) || isstring(x));
p.addParameter('ColorLimits', 'auto', @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && numel(x) == 2));
p.addParameter('ColorMap', 'rdylbu', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('GrayBuffer', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('MontageType', 'hcp', @(x) ischar(x) || isstring(x));
p.parse(img_obj_or_file, out_video, varargin{:});
opts = p.Results;

if isa(img_obj_or_file, 'image_vector')
    dat = img_obj_or_file;
else
    dat = fmri_data(char(img_obj_or_file));
end

if ~isempty(opts.Threshold)
    dat = threshold(dat, [opts.Threshold Inf], 'raw-between');
end

% Compute the color range ONCE over the whole animation.
cmaprange = local_resolve_clim(dat, opts);

% Build pos / neg colormaps; canlab's montage/addblobs splits these for
% diverging maps. RdYlBu is the diverging ColorBrewer palette; hotcool
% reproduces canlab's classic default.
[pos_colormap, neg_colormap] = local_resolve_colormaps(opts.ColorMap);

% These options flow through montage -> addblobs (volumetric) and
% render_on_surface (surface) so the blob color mapping and gray buffer
% are identical across every frame of the animation.
montage_args = { ...
    'cmaprange', cmaprange, ...
    'pos_colormap', pos_colormap, ...
    'neg_colormap', neg_colormap, ...
    'gray_buffer', logical(opts.GrayBuffer)};

vw = VideoWriter(char(out_video), 'MPEG-4');
vw.FrameRate = opts.FPS;
open(vw);

fig = figure('Color', 'w', 'Visible', 'off');
for t = 1:opts.FrameStep:size(dat.dat, 2)
    clf(fig);
    this = get_wh_image(dat, t);
    if logical(opts.UseSigMask) && isa(this, 'statistic_image') && ~isempty(this.sig)
        this.dat(~logical(this.sig)) = 0;
    end
    montage(this, char(opts.MontageType), montage_args{:});
    title(sprintf('%s %d', char(opts.TitlePrefix), t), 'Interpreter', 'none');
    frame = getframe(fig);
    writeVideo(vw, frame);
end
close(vw);
close(fig);
end


function clim = local_resolve_clim(dat, opts)
% Pick a [lo hi] CLim that's stable across the whole animation.
% 'auto' -> symmetric range from the 98th percentile of |non-zero, finite|
% values across all frames (so outliers don't blow out the scale).
if isnumeric(opts.ColorLimits) && numel(opts.ColorLimits) == 2
    clim = double(opts.ColorLimits(:)');
    return
end

vals = double(dat.dat);
if isa(dat, 'statistic_image') && logical(opts.UseSigMask) && ~isempty(dat.sig)
    vals(~logical(dat.sig)) = 0;
end
vals = vals(isfinite(vals) & vals ~= 0);
if isempty(vals)
    clim = [-1, 1];
    return
end
q = quantile(abs(vals), 0.98);
if ~(q > 0)
    q = max(abs(vals));
end
if ~(q > 0)
    q = 1;
end
clim = [-q, q];
end


function [pos_cm, neg_cm] = local_resolve_colormaps(spec)
% Returns pos_colormap (low-positive -> high-positive) and
% neg_colormap (low-negative -> high-negative) Nx3 matrices.
% canlab's addblobs interpolates between the first and last rows.
if iscell(spec) && numel(spec) == 2
    pos_cm = spec{1};
    neg_cm = spec{2};
    return
end
name = lower(strtrim(char(spec)));
switch name
    case 'rdylbu'
        % ColorBrewer RdYlBu, split at yellow:
        %   positive half: pale yellow -> orange -> dark red
        %   negative half: pale blue   -> medium blue -> dark blue
        pos_cm = colormap_tor([1.0000 1.0000 0.7490], [0.6471 0.0000 0.1490], ...
                              [0.9922 0.6824 0.3804]);
        neg_cm = colormap_tor([0.8784 0.9529 0.9725], [0.1922 0.2118 0.5843], ...
                              [0.4549 0.6784 0.8196]);
    case 'hotcool'
        % canlab's classic default: yellow -> red for pos, cyan -> blue for neg
        pos_cm = colormap_tor([1 1 0], [1 0 0]);
        neg_cm = colormap_tor([0 1 1], [0 0 1]);
    case 'bwr'
        % blue-white-red diverging
        pos_cm = colormap_tor([1 1 1], [1 0 0]);
        neg_cm = colormap_tor([1 1 1], [0 0 1]);
    otherwise
        error('hrf_make_montage_animation:UnknownColorMap', ...
            'Unknown ColorMap: ''%s''. Use ''rdylbu'', ''hotcool'', ''bwr'', or a {pos_cm, neg_cm} pair.', ...
            name);
end
end
