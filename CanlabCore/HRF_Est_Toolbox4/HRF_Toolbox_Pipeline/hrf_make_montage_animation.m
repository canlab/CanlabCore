function out_video = hrf_make_montage_animation(img_obj_or_file, out_video, varargin)
%HRF_MAKE_MONTAGE_ANIMATION Fast montage-based 3D animation over time.
% Supports raw volume, beta maps, or t maps in fmri_data format.
%
% Color limits are held constant across frames so the animation has a
% consistent interpretation (the default per-frame montage CLim would
% rescale every frame independently and obscure amplitude changes).
%
% Name-value parameters
% ---------------------
%   'FrameStep'    - integer, animate every Nth volume (default 2)
%   'FPS'          - playback rate (default 8)
%   'Threshold'    - scalar; voxels with |value| < Threshold are zeroed
%                    (via statistic_image/threshold). Default [] (no thresh).
%   'UseSigMask'   - true (default) to honor obj.sig for statistic_image input
%   'TitlePrefix'  - prefix for the per-frame title
%   'ColorLimits'  - 'auto' (default) -> symmetric robust range from data,
%                    or [lo hi] to set explicitly. Held constant across
%                    every frame of the animation.

p = inputParser;
p.addRequired('img_obj_or_file');
p.addRequired('out_video', @(x) ischar(x) || isstring(x));
p.addParameter('FrameStep', 2, @(x) isscalar(x) && x >= 1);
p.addParameter('FPS', 8, @(x) isscalar(x) && x > 0);
p.addParameter('Threshold', [], @(x) isempty(x) || isscalar(x));
p.addParameter('UseSigMask', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('TitlePrefix', 'Frame', @(x) ischar(x) || isstring(x));
p.addParameter('ColorLimits', 'auto', @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && numel(x) == 2));
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

% Compute the color range ONCE so every frame uses the same scale.
clim = local_resolve_clim(dat, opts);

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
    % montage(this);
    montage(this, 'hcp');
    title(sprintf('%s %d', char(opts.TitlePrefix), t), 'Interpreter', 'none');
    local_apply_clim(fig, clim);
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


function local_apply_clim(fig, clim)
% Override montage()'s per-frame autoscaling. montage() may not honour CLim
% set ahead of time on the parent figure, so we set it on every axes after
% each frame is drawn.
axes_in_fig = findall(fig, 'Type', 'axes');
for a = 1:numel(axes_in_fig)
    if isempty(get(axes_in_fig(a), 'Children'))
        continue
    end
    try
        set(axes_in_fig(a), 'CLim', clim);
    catch
    end
end
end
