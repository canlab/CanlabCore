function out_video = hrf_make_montage_animation(img_obj_or_file, out_video, varargin)
%HRF_MAKE_MONTAGE_ANIMATION Fast montage-based 3D animation over time.
% Supports raw volume, beta maps, or t maps in fmri_data format.

p = inputParser;
p.addRequired('img_obj_or_file');
p.addRequired('out_video', @(x) ischar(x) || isstring(x));
p.addParameter('FrameStep', 2, @(x) isscalar(x) && x >= 1);
p.addParameter('FPS', 8, @(x) isscalar(x) && x > 0);
p.addParameter('Threshold', [], @(x) isempty(x) || isscalar(x));
p.addParameter('UseSigMask', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('TitlePrefix', 'Frame', @(x) ischar(x) || isstring(x));
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
    frame = getframe(fig);
    writeVideo(vw, frame);
end
close(vw);
close(fig);
end
