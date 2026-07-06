function wt_save(figh, fname, varargin)
% Capture a figure (or uifigure controller) to a tightly-cropped PNG for the
% visualization walkthrough docs.
%
% wt_save(figh, fname)                 % regular figure via exportgraphics
% wt_save(figh, fname, 'app')          % uifigure (controller) via exportapp
% wt_save(figh, fname, 'pad', 12)      % white padding (px) around content
%
% Trims surrounding white so montages/surfaces/controllers come out clean
% regardless of the empty layout space CANlab figures leave.

isapp = any(strcmpi(varargin, 'app'));
nobars = any(strcmpi(varargin, 'nobars'));      % strip colorbar legends (cleaner surfaces)
pad = 10;  wh = find(strcmpi(varargin, 'pad'), 1);          if ~isempty(wh), pad = varargin{wh+1}; end
res = 130; wh = find(strcmpi(varargin, 'resolution'), 1);   if ~isempty(wh), res = varargin{wh+1}; end

figdir = fileparts(fname);
if ~isempty(figdir) && ~exist(figdir, 'dir'), mkdir(figdir); end

if nobars
    % Delete colorbar legends and the (now-empty) axes that hosted them, so the
    % exported surface is clean. render_on_surface builds these with colorbar().
    delete(findall(figh, 'Type', 'ColorBar'));
end

tmp = [tempname '.png'];
drawnow
if isapp
    exportapp(figh, tmp);
else
    exportgraphics(figh, tmp, 'Resolution', res, 'BackgroundColor', 'white');
end

% Autocrop white borders
im = imread(tmp);
delete(tmp);
mask = min(im, [], 3) < 248;           % non-white content (white = 255 all channels)
if ~any(mask(:)), imwrite(im, fname); return, end
% Require a minimum number of content pixels per kept row/col, so stray faint
% anti-aliasing dots don't defeat the crop (montage figures leave big margins).
minpix = 4;
rows = find(sum(mask, 2) >= minpix);  cols = find(sum(mask, 1) >= minpix);
if isempty(rows) || isempty(cols)     % fall back to any-content if too aggressive
    rows = find(any(mask, 2));  cols = find(any(mask, 1));
end
r1 = max(1, rows(1) - pad);   r2 = min(size(im, 1), rows(end) + pad);
c1 = max(1, cols(1) - pad);   c2 = min(size(im, 2), cols(end) + pad);
imwrite(im(r1:r2, c1:c2, :), fname);

info = imfinfo(fname);
fprintf('  saved %s  (%dx%d)\n', fname, info.Width, info.Height);
end
