function movie_file = hrf_animate_wholebrain_stats(stats_input, varargin)
%HRF_ANIMATE_WHOLEBRAIN_STATS Make a quick montage movie from 4D HRF maps.
%
% movie_file = hrf_animate_wholebrain_stats(stats_input, ...)
%
% stats_input may be hrf_fit_wholebrain_stats output, a statistic_image
% object, or a 4D NIfTI filename. The function steps through volumes with
% CANlab montage(), which keeps it compatible with statistic_image
% thresholding.

p = inputParser;
p.addRequired('stats_input');
p.addParameter('Object', 't', @(x) ischar(x) || isstring(x));
p.addParameter('OutputFile', 'hrf_montage_animation.mp4', @(x) ischar(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x));
p.addParameter('MetadataTable', table(), @(x) isempty(x) || istable(x));
p.addParameter('UseThreshold', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('FrameRate', 4, @(x) isscalar(x) && x > 0);
p.addParameter('Clim', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
p.addParameter('MontageArgs', {}, @(x) iscell(x));
p.addParameter('FigureVisible', 'off', @(x) ischar(x) || isstring(x));
p.parse(stats_input, varargin{:});
opts = p.Results;

[obj, metadata_table] = local_get_object(stats_input, opts.Object, opts.MetadataTable);
wh = local_select_volumes(obj, metadata_table, char(opts.Condition));
movie_file = char(opts.OutputFile);

if isempty(wh)
    error('No volumes selected for animation.');
end

v = VideoWriter(movie_file, local_video_profile(movie_file));
v.FrameRate = opts.FrameRate;
open(v);
cleaner = onCleanup(@() close(v));

fig = figure('Color', 'w', 'Visible', char(opts.FigureVisible));
for i = 1:numel(wh)
    clf(fig);
    frame_obj = get_wh_image(obj, wh(i));
    if logical(opts.UseThreshold) && isa(frame_obj, 'statistic_image') && ~isempty(frame_obj.sig)
        frame_obj.dat(~frame_obj.sig) = 0;
    end

    montage(frame_obj, opts.MontageArgs{:});
    if ~isempty(opts.Clim)
        try
            set(findobj(fig, 'Type', 'axes'), 'CLim', opts.Clim);
        catch
        end
    end
    title(local_frame_title(frame_obj, metadata_table, wh(i)), 'Interpreter', 'none');
    drawnow;
    writeVideo(v, getframe(fig));
end
close(fig);
end

function [obj, metadata_table] = local_get_object(stats_input, which_obj, metadata_table)
if isstruct(stats_input) && isfield(stats_input, 'b') && isfield(stats_input, 't')
    switch lower(char(which_obj))
        case {'beta', 'b'}
            obj = stats_input.b;
        case {'t', 'tmap', 'tmaps'}
            obj = stats_input.t;
        otherwise
            error('Unknown Object: %s. Use ''beta'' or ''t''.', char(which_obj));
    end
    if isempty(metadata_table) && isfield(stats_input, 'metadata_table')
        metadata_table = stats_input.metadata_table;
    end
elseif isa(stats_input, 'image_vector')
    obj = stats_input;
elseif ischar(stats_input) || isstring(stats_input)
    obj = statistic_image(char(stats_input), 'type', 'generic');
else
    error('Unsupported stats_input type.');
end
end

function wh = local_select_volumes(obj, metadata_table, condition)
n_images = size(obj.dat, 2);
wh = 1:n_images;
if isempty(condition)
    return
end

if ~isempty(metadata_table) && any(strcmp('condition', metadata_table.Properties.VariableNames))
    wh = find(strcmp(cellstr(string(metadata_table.condition)), condition));
elseif isa(obj, 'statistic_image') && ~isempty(obj.image_labels)
    wh = find(contains(cellstr(string(obj.image_labels)), condition));
else
    error('Condition selection needs metadata_table.condition or statistic_image.image_labels.');
end
end

function ttl = local_frame_title(obj, metadata_table, idx)
ttl = sprintf('Volume %d', idx);
if ~isempty(metadata_table) && height(metadata_table) >= idx && any(strcmp('image_label', metadata_table.Properties.VariableNames))
    ttl = char(metadata_table.image_label(idx));
elseif isa(obj, 'statistic_image') && ~isempty(obj.image_labels)
    ttl = obj.image_labels{1};
end
end

function profile = local_video_profile(movie_file)
[~, ~, ext] = fileparts(movie_file);
switch lower(ext)
    case '.avi'
        profile = 'Motion JPEG AVI';
    otherwise
        profile = 'MPEG-4';
end
end
