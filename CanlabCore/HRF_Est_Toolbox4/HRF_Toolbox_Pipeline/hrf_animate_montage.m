function out = hrf_animate_montage(panels, varargin)
%HRF_ANIMATE_MONTAGE Synchronized grid of HRF animations (word-clouds + brains).
%
% Renders several per-lag animations side by side, advancing them TOGETHER so
% every tile shows the same peristimulus lag at once -- e.g. a Neurosynth-term
% word-cloud next to a signature word-cloud next to a whole-brain montage, all
% synchronized over the HRF. Each panel produces its own per-lag frames; the
% driver tiles the l-th frame of every panel into one composite movie.
%
% :Usage:
% ::
%     panels = {
%       struct('type','wordcloud', 'args', {{dirs, 'Set','neurosynth', 'Condition','heat'}}), ...
%       struct('type','wordcloud', 'args', {{dirs, 'Unit','signature', 'Set','all', 'Condition','heat'}}), ...
%       struct('type','movie', 'file', 'group_heat_average.mp4', 'title','brain') };
%     hrf_animate_montage(panels, 'OutputFile','synced.mp4', 'FrameRate',5);
%
% :Inputs:
%   **panels:** cell array of panel specs (each a struct):
%     .type = 'wordcloud' -> rendered by hrf_animate_wordcloud; pass its
%                            name-value inputs (incl. the source) in .args (a
%                            cell). Do NOT set its OutputFile.
%           = 'movie'/'brain' -> read per-lag frames from an existing movie/gif
%                            in .file (e.g. a brain montage from
%                            hrf_make_montage_animation / hrf_pooled_wholebrain_animation).
%           = 'frames' -> use .frames directly (cell of [H x W x 3] uint8).
%     .title (optional) - label drawn above the tile.
%
% :Optional Inputs:
%   **'OutputFile':** '.mp4'/'.avi'/'.gif' to write. Default ''.
%   **'FrameRate':**  fps. Default 4.
%   **'Layout':**     [rows cols] tile grid. Default a near-square auto layout.
%   **'Title':**      super-title drawn across the top. Default ''.
%   **'Verbose':**    default true.
%
% :Output:
%   **out:** struct with .nlag, .npanel, .frames (composite RGB per lag),
%            .file (written path or '').
%
% See also: hrf_animate_wordcloud, hrf_make_montage_animation,
%           hrf_pooled_wholebrain_animation.

p = inputParser;
p.addRequired('panels', @(x) iscell(x) && ~isempty(x));
p.addParameter('OutputFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('FrameRate', 4, @(x) isscalar(x) && x > 0);
p.addParameter('Layout', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
p.addParameter('Title', '', @(x) ischar(x) || isstring(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(panels, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);

% ---- render / collect each panel's per-lag frames -----------------------
np = numel(panels);
pf = cell(1, np); ptitle = cell(1, np);
for i = 1:np
    [pf{i}, ptitle{i}] = local_panel_frames(panels{i}, verbose, i);
end
nlags = cellfun(@numel, pf);
valid = nlags > 0;
if ~any(valid), error('hrf_animate_montage:NoFrames', 'No panel produced frames.'); end
if any(~valid)
    warning('hrf_animate_montage:EmptyPanel', 'Dropping %d panel(s) that produced no frames.', sum(~valid));
    pf = pf(valid); ptitle = ptitle(valid); nlags = nlags(valid); np = numel(pf);
end
nlag = min(nlags);
if numel(unique(nlags)) > 1 && verbose
    fprintf('hrf_animate_montage: panels have %s frames; using the common first %d.\n', mat2str(nlags), nlag);
end

% ---- common tile size (pad, no distortion) + grid layout ----------------
tileH = 0; tileW = 0;
for i = 1:np
    sz = size(pf{i}{1}); tileH = max(tileH, sz(1)); tileW = max(tileW, sz(2));
end
lab_h = 26;                                   % strip above each tile for its title
if isempty(opts.Layout)
    rows = floor(sqrt(np)); cols = ceil(np / max(rows, 1));
else
    rows = opts.Layout(1); cols = opts.Layout(2);
end

vw = local_open_video(opts.OutputFile, opts.FrameRate);
comp_frames = cell(1, nlag);
for l = 1:nlag
    tiles = cell(rows, cols);
    for i = 1:np
        tiles{i} = local_label_tile(local_pad_center(pf{i}{l}, tileH, tileW), ptitle{i}, lab_h);
    end
    blank = 255 * ones(tileH + lab_h, tileW, 3, 'uint8');
    for k = np + 1:rows * cols, tiles{k} = blank; end
    rowimgs = cell(rows, 1);
    for r = 1:rows, rowimgs{r} = cat(2, tiles{r, :}); end
    comp = cat(1, rowimgs{:});
    if ~isempty(char(opts.Title)), comp = local_super_title(comp, char(opts.Title)); end
    comp_frames{l} = comp;
    vw = local_write_frame(vw, comp);
end
vw = local_close_video(vw);

out = struct('nlag', nlag, 'npanel', np, 'frames', {comp_frames}, ...
    'file', local_out_path(vw, opts.OutputFile));
if verbose
    fprintf('hrf_animate_montage: %d panels x %d lags, %dx%d grid%s\n', np, nlag, rows, cols, ...
        local_note(out.file));
end
end


% =========================================================================
function [frames, ttl] = local_panel_frames(spec, verbose, idx)
frames = {}; ttl = '';
if ~isstruct(spec), error('hrf_animate_montage:Spec', 'Each panel must be a struct.'); end
if isfield(spec, 'title'), ttl = char(string(spec.title)); end
type = 'wordcloud'; if isfield(spec, 'type'), type = lower(char(spec.type)); end
switch type
    case 'wordcloud'
        args = {}; if isfield(spec, 'args'), args = spec.args; end
        if verbose, fprintf('  [panel %d] wordcloud ...\n', idx); end
        o = hrf_animate_wordcloud(args{:}, 'ReturnFrames', true, 'OutputFile', '', 'doverbose', false);
        frames = o.frames;
        if isempty(ttl), ttl = o.title; end
    case {'movie', 'brain'}
        if ~isfield(spec, 'file'), error('hrf_animate_montage:NoFile', 'movie/brain panel needs a .file.'); end
        if verbose, fprintf('  [panel %d] frames from %s ...\n', idx, char(spec.file)); end
        frames = local_read_frames(char(spec.file));
    case 'frames'
        if ~isfield(spec, 'frames'), error('hrf_animate_montage:NoFrames2', 'frames panel needs .frames.'); end
        frames = spec.frames;
    otherwise
        error('hrf_animate_montage:Type', 'Unknown panel type: %s', type);
end
frames = frames(~cellfun(@isempty, frames));
end


function img = local_pad_center(img, H, W)
img = local_u8(img);
if size(img, 3) == 1, img = repmat(img, 1, 1, 3); end
[h, w, ~] = size(img);
if h > H || w > W                              % shouldn't happen (H,W are maxima) but guard
    img = img(1:min(h, H), 1:min(w, W), :); [h, w, ~] = size(img);
end
canvas = 255 * ones(H, W, 3, 'uint8');
r0 = floor((H - h) / 2) + 1; c0 = floor((W - w) / 2) + 1;
canvas(r0:r0 + h - 1, c0:c0 + w - 1, :) = img;
img = canvas;
end


function tile = local_label_tile(img, ttl, lab_h)
strip = 255 * ones(lab_h, size(img, 2), 3, 'uint8');
tile = cat(1, strip, img);
if isempty(ttl), return; end
tile = local_puttext(tile, [size(img, 2) / 2, lab_h / 2], ttl, 16);
end


function comp = local_super_title(comp, ttl)
strip = 255 * ones(34, size(comp, 2), 3, 'uint8');
strip = local_puttext(strip, [size(comp, 2) / 2, 17], ttl, 20);
comp = cat(1, strip, comp);
end


function img = local_puttext(img, pos, str, fs)
% Draw centered black text via insertText if available; otherwise leave the
% image unchanged (labels are cosmetic, no Computer Vision Toolbox required).
if exist('insertText', 'file') ~= 2, return; end
try
    img = insertText(img, pos, str, 'AnchorPoint', 'Center', 'FontSize', fs, ...
        'BoxOpacity', 0, 'TextColor', 'black');
catch
end
end


function y = local_u8(x)
if isa(x, 'uint8'), y = x; return; end
x = double(x);
if max(x(:)) <= 1 + eps, x = x * 255; end
y = uint8(min(max(x, 0), 255));
end


function frames = local_read_frames(file)
if exist(file, 'file') ~= 2, error('hrf_animate_montage:MissingFile', 'Not found: %s', file); end
[~, ~, e] = fileparts(file);
if strcmpi(e, '.gif')
    [g, map] = imread(file, 'frames', 'all');
    n = size(g, 4); frames = cell(1, n);
    for i = 1:n, frames{i} = local_u8(ind2rgb(g(:, :, 1, i), map)); end
else
    v = VideoReader(file); frames = {};
    while hasFrame(v), frames{end + 1} = local_u8(readFrame(v)); end %#ok<AGROW>
end
end


% ---- video helpers (mirror hrf_animate_wordcloud) -----------------------
function vw = local_open_video(file, fps)
vw = struct('type', 'none', 'obj', [], 'file', char(file), 'fps', fps, 'frames', {{}});
f = char(file); if isempty(f), return; end
[~, ~, e] = fileparts(f);
if strcmpi(e, '.gif')
    vw.type = 'gif';
elseif strcmpi(e, '.avi')
    vw.obj = VideoWriter(f, 'Motion JPEG AVI'); vw.obj.FrameRate = fps; open(vw.obj); vw.type = 'avi';
else
    try, vw.obj = VideoWriter(f, 'MPEG-4'); catch, f = regexprep(f, '\.mp4$', '.avi'); vw.file = f; vw.obj = VideoWriter(f, 'Motion JPEG AVI'); end
    vw.obj.FrameRate = fps; open(vw.obj); vw.type = 'mp4';
end
end

function vw = local_write_frame(vw, img)
if strcmp(vw.type, 'none'), return; end
if strcmp(vw.type, 'gif'), vw.frames{end + 1} = img; else, writeVideo(vw.obj, local_u8(img)); end
end

function vw = local_close_video(vw)
if strcmp(vw.type, 'gif')
    F = vw.frames; if isempty(F), return; end
    [~, gmap] = rgb2ind(cat(1, F{:}), 256); dt = 1 / max(vw.fps, 1);
    for i = 1:numel(F)
        A = rgb2ind(F{i}, gmap);
        if i == 1, imwrite(A, gmap, vw.file, 'gif', 'LoopCount', Inf, 'DelayTime', dt);
        else, imwrite(A, gmap, vw.file, 'gif', 'WriteMode', 'append', 'DelayTime', dt); end
    end
elseif ~strcmp(vw.type, 'none') && ~isempty(vw.obj)
    close(vw.obj);
end
end

function pth = local_out_path(vw, file)
if isempty(char(file)), pth = ''; else, pth = vw.file; end
end

function s = local_note(f)
if isempty(f), s = ' (not saved)'; else, s = sprintf(' -> %s', f); end
end
