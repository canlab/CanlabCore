function [whsl, plate] = display_slices(dat, varargin)
% Quick montage view of a single volume in an image_vector / fmri_data object.
%
% - By default, displays an axial montage. The slice spacing is chosen
%   automatically so that roughly 24 non-empty slices are shown (fewer if
%   the volume contains fewer than 24 non-empty slices in that direction).
% - Pass 'saggital'/'sagittal' or 'coronal' to switch to one of those
%   views; the same auto slice-count logic is applied.
% - Pass an explicit 'spacing' value to override the auto-spacing.
% - Pass 'three_views' to get the legacy combined axial + coronal +
%   sagittal figure that used to be the default.
% - Pass 'multi_image' for objects with multiple images in `.dat` to get
%   one slice per image, tiled as subplots. For statistic_image inputs a
%   canonical T1 underlay is loaded and significant voxels are drawn in
%   hot/cool split colors on top; for plain fmri_data / image_vector the
%   slice is rendered directly with the fast imagesc style used by the
%   single-volume montage.
% - Also easy to select slices: Enter 'startslice' and 'endslice', each followed by values, to display
%   a range of slices. Default is entry in mm coordinate, but 'voxelspace' allows you to specify them in
%   voxel coordinates, i.e., in slices in the actual data image.
%
% - This function is good for getting a quick view of the contents of an
%   image, e.g., for structural images.
% - The montage.m method and fmridisplay are good for displaying blobs on
%   standard brains, and are much slower.
%
% :Usage:
% ::
%
%    [whsl, data_matrix] = display_slices(dat, [myview], ['spacing', slicespacing], ['vertical'])
%
% :Inputs:
%
% **'axial', 'saggital', 'sagittal', 'coronal':**
%   Keywords to control view (default: 'axial')
%
% **'three_views'**
%   Draw the legacy 3-panel figure (axial + coronal + sagittal) in a
%   dedicated figure window. Equivalent to the pre-2026 default behavior.
%
% **'multi_image'**
%   Tile one slice per image in `.dat` instead of a slice-montage of one
%   image. Combine with 'axial' / 'coronal' / 'saggital' to choose
%   orientation (default axial), and 'slice', followed by an mm
%   coordinate, to choose the slice location (default 0). For
%   statistic_image inputs, a canonical T1 (via canlab_get_underlay_image)
%   is rendered as a gray underlay with significant voxels drawn in
%   hot/cool split colors on top.
%
% **'slice'**
%   Followed by mm coordinate of the slice to show in 'multi_image' mode
%   (default 0). Ignored in single-volume / three_views modes.
%
% **'nimages'**
%   Followed by maximum number of images to render in 'multi_image' mode.
%   Default is all columns of `.dat`, capped at 64.
%
% **'names'**
%   Followed by a cell array of strings used as subplot titles in
%   'multi_image' mode.
%
% **'spacing'**
%   Followed by slice spacing in mm (or voxels with 'voxelspace').
%   If omitted, spacing is chosen automatically to give ~24 slices.
%
% **'vertical'**
%   Vertical stack instead of horizontal
%
% **'slices_per_row'**
%   Followed by slices per row in mm
%
% **'startslice'**
%   Followed by start slice in mm
%
% **'endslice'**
%   Followed by ending slice in mm
%
% **'mm'**
%   Enter input coordinates in mm
%
% **'voxelspace'**
%   Enter input coords in voxels
%
% **'clim'**
%   Followed by color limit for plot
%
% :Outputs:
%
%   **whsl:**
%        Vector of which slices have been chosen
%
%   **data_matrix:**
%        Concatenated data matrix to be imaged
%
% :Examples:
% ::
%
%     % --- single-volume views (use a single image, e.g., the mean) ---
%     imgs = load_image_set('emotionreg');           % sample multi-image fmri_data
%     m    = mean(imgs);                              % single-image volume
%     display_slices(m);                              % axial montage, auto slice count (~24)
%     display_slices(m, 'coronal');                   % coronal montage, auto slice count
%     display_slices(m, 'saggital');                  % sagittal montage, auto slice count
%     display_slices(m, 'three_views');               % legacy 3-panel combined figure
%     display_slices(m, 'axial', 'spacing', 8);                % override: 8 mm spacing
%     display_slices(m, 'saggital', 'spacing', 10, 'vertical'); % vertical stack
%     display_slices(m, 'axial', 'slices_per_row', 20, 'spacing', 4, 'startslice', -10, 'endslice', 10);
%
%     % --- multi-image mode: one slice per image ---
%     display_slices(imgs, 'multi_image');                          % axial, z = 0 mm
%     display_slices(imgs, 'multi_image', 'slice', -10);            % axial, z = -10 mm
%     display_slices(imgs, 'multi_image', 'coronal', 'slice', 0);   % coronal, y = 0
%     display_slices(imgs, 'multi_image', 'saggital', 'slice', -4); % sagittal, x = -4
%     display_slices(imgs, 'multi_image', 'nimages', 6);            % only first 6
%
%     % --- multi-image, statistic_image (1 image): T1 underlay + blobs ---
%     t = ttest(imgs); t = threshold(t, .01, 'unc');                % statistic_image
%     display_slices(t, 'multi_image');                             % axial, z = 0
%     display_slices(t, 'multi_image', 'coronal', 'slice', -20);
%
%     % --- multi-image, statistic_image with several images in .dat ---
%     % Build a 2-column statistic_image by t-testing two subgroups, so
%     % t_multi.dat(:, 2) holds the second contrast.
%     half1 = ttest(get_wh_image(imgs, 1:floor(end/2)));
%     half2 = ttest(get_wh_image(imgs, floor(end/2)+1:size(imgs.dat,2)));
%     half1 = threshold(half1, .01, 'unc');
%     half2 = threshold(half2, .01, 'unc');
%     t_multi             = half1;
%     t_multi.dat         = [half1.dat,  half2.dat];   % t_multi.dat(:,2) = 2nd map
%     t_multi.p           = [half1.p,    half2.p];
%     t_multi.sig         = [half1.sig,  half2.sig];
%     t_multi.image_labels = {'Half 1', 'Half 2'};
%     display_slices(t_multi, 'multi_image', ...
%         'names', t_multi.image_labels);                            % 2 underlays + blobs
%     display_slices(t_multi, 'multi_image', 'coronal', 'slice', -10);
%
% Make a montage of a subset of slices and make it fill the figure:
% create_figure('sagg')
% display_slices(m, 'saggital', 'slices_per_row', 4, 'spacing', 3, 'startslice', -20, 'endslice', 20)
% set(gca, 'Position', [.05 .05 .9 1]);
% colormap gray
%
% ..
%    Tor Wager, Aug 2012. Update: Jan 2018, May 2026
% ..

do3views = false; % legacy 3-view composite (opt-in via 'three_views')
domulti_image = false; % multi-image mode (one slice per image)
myview = 'axial'; % default for single view
spacing = 1;
user_set_spacing = false;  % track whether caller passed an explicit 'spacing'
target_nslices = 24;       % target slice count when auto-spacing
dovertical = false;
stackorient = 2;  % 2 for horiz, 1 for vertical
dotight = true;   % tight: show only non-zero slices
s = 12;           % slices per row
% startslice = [];  % first slice to show; set later
% endslice = [];    % last slice to show (empty = all); set later
entered_mm_coords = true;  % interpret entries in mm or voxel space coordinates
clim = [];
slice_mm = 0;        % mm-coord for multi_image mode
nimages_req = [];    % requested image count for multi_image mode (empty = all)
image_names = {};    % subplot titles in multi_image mode
do_plot = true;      % set false via 'noplot' to suppress imagesc (build plate only)

% Process inputs - single case
% in current axes
% ------------------------------------------------------------------------
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case {'axial', 'saggital', 'sagittal', 'coronal'}, myview = varargin{i};
            case {'three_views', 'threeviews', '3views', 'multiview'}, do3views = true;
            case {'multi_image', 'multiimage'}, domulti_image = true;
            case 'slice', slice_mm = varargin{i + 1};
            case 'nimages', nimages_req = varargin{i + 1};
            case 'names', image_names = varargin{i + 1};
            case 'spacing', spacing = varargin{i + 1}; user_set_spacing = true;
            case 'vertical', stackorient = 1;
            case 'slices_per_row', s = varargin{i + 1};
            case 'startslice', startslice = varargin{i + 1};
            case 'endslice', endslice = varargin{i + 1};
            case 'mm', entered_mm_coords = true;            % we have entered coords in mm [default]
            case 'voxelspace', entered_mm_coords = false;   % entered values in voxel space

            case 'clim', clim = varargin{i + 1};
            case {'noplot', 'noimage'}, do_plot = false;

            otherwise
                disp('Warning: Unknown input.')
        end
    end
end

% wh = strcmp(varargin, 'voxelspace'); % entered in voxel space
% if any(wh), entered_mm_coords = false; end

if ~domulti_image && size(dat.dat, 2) > 1
    error('Use display_slices on objects containing only one image, or pass ''multi_image'' for a per-image slice montage.');
end

% ------------------------------------------------------------------------
% Three-views mode: render axial, coronal, sagittal as three horizontal
% montages and stack them vertically into a single full-figure axes. Each
% view uses auto slice spacing (~24 slices) and the default 12 slices per
% row, so each "stack" is a compact 2-row strip. Padding goes on the
% right of narrower rows so the composite is rectangular.
%
if do3views
    [~, ax_plate]  = display_slices(dat, 'axial',    'noplot');
    [~, cor_plate] = display_slices(dat, 'coronal',  'noplot');
    [~, sag_plate] = display_slices(dat, 'saggital', 'noplot');

    plate = compose_three_views(ax_plate, cor_plate, sag_plate);

    fh = create_figure('slice_display');
    set(fh, 'Color', 'w');
    ax_h = axes('Parent', fh, 'Position', [0 0 1 1], 'Units', 'normalized');
    imagesc(ax_h, plate);
    axis(ax_h, 'image');
    set(ax_h, 'YDir', 'Reverse');
    axis(ax_h, 'off');

    apply_split_colormap_and_clim(ax_h, plate, clim);
    fill_axes_in_figure(ax_h, plate);

    whsl = [];  % three_views produces no single slice list
    return
end

% ------------------------------------------------------------------------
% Multi-image mode: one slice per image in .dat
% ------------------------------------------------------------------------
if domulti_image
    [whsl, plate] = display_slices_multi_image(dat, myview, slice_mm, ...
        nimages_req, image_names, clim);
    return
end

% ------------------------------------------------------------------------
% Main function for single montage
%

if stackorient == 2, platestack = 1; else platestack = 1; end

vdat = reconstruct_image(dat);

% flip x dim if x voxel size is negative. other display/write methods also
% do this.
if sign(dat.volInfo.mat(1)) < 0
    vdat = flipdim(vdat, 1);
end

% Find out which dimension to vary slices along
% ------------------------------------------------------------------------
switch myview
    case 'axial'
        wh_col = 3;
    case {'saggital', 'sagittal'}
        wh_col = 1;
    case 'coronal'
        wh_col = 2;
    otherwise error('unknown slice view.')
end

% Get vector of slices in mm or voxel space, depending on input args
% ------------------------------------------------------------------------
% re-set default slices if we have not entered them
% depends on voxel or mm. if mm, convert to voxels.
if ~any(strcmp(varargin, 'startslice'))
    
    if entered_mm_coords  % Figure out mm coords of first slice
        
        xyzvoxel = voxel2mm([1 1 1]', dat.volInfo.mat);
        startmm = xyzvoxel;
        startslice = xyzvoxel(wh_col); 
    else
        startslice = 1;
    end
end

if ~any(strcmp(varargin, 'endslice'))
    if entered_mm_coords % Figure out mm coords of last slice
        
        xyzvoxel = voxel2mm(size(vdat)', dat.volInfo.mat);
        endmm = xyzvoxel;
        endslice = xyzvoxel(wh_col);
    else
        endslice = size(vdat, wh_col);
    end
end

% reverse slices if radiological orientation (-x vox)
if startslice > endslice

    tmp = endslice;
    endslice = startslice;
    startslice = tmp;

end

% Auto-pick spacing if the caller didn't provide one, so the montage shows
% at most ~target_nslices non-empty slices through the volume.
% ------------------------------------------------------------------------
if ~user_set_spacing
    switch wh_col
        case 1
            sumabsval = squeeze(sum(sum(abs(vdat), 2), 3));
        case 2
            sumabsval = squeeze(sum(sum(abs(vdat), 1), 3));
        case 3
            sumabsval = squeeze(sum(sum(abs(vdat), 1), 2));
    end

    nonempty = find(sumabsval >= 100 * eps);

    if isempty(nonempty)
        spacing_voxels = 1;
    else
        extent_voxels = nonempty(end) - nonempty(1) + 1;
        spacing_voxels = max(1, ceil(extent_voxels / target_nslices));
    end

    if entered_mm_coords
        voxsize = abs(dat.volInfo.mat(wh_col, wh_col));
        if voxsize == 0, voxsize = 1; end
        spacing = spacing_voxels * voxsize;
    else
        spacing = spacing_voxels;
    end
end

whsl = startslice:spacing:endslice;  % voxel or mm


% mm to voxel conversion
% ------------------------------------------------------------------------
if entered_mm_coords
    % Convert from mm 2 voxel, if we have entered mm coordinates
    
    xyzmm = 10 .* ones(length(whsl), 3);
    
    xyzmm(:, wh_col) = whsl;
    
    xyzvoxel = mm2voxel(xyzmm', dat.volInfo.mat, 1);
    whsl = xyzvoxel(:, wh_col)';

    % ISSUE: in some cases, xyzvoxel is < 1
    if any(whsl < 1)
        disp('Warning! Interpolation issue: Some slices < 1')
        whsl(whsl < 1) = 1;
        
        sz = size(vdat);
        whsl(whsl > sz(wh_col)) = sz(wh_col);
        whsl = unique(whsl);
        
    end
    
end

% Now everything is in voxel coords from here on out.

if dotight
   % figure out which slices to skip because they are empty.
   switch myview
       case 'axial'
           subvol = abs(vdat(:, :, whsl));
           sumabsval = sum(sum(subvol, 1), 2);
           
       case {'saggital', 'sagittal'}
           subvol = abs(vdat(whsl, :, :));
           sumabsval = sum(sum(subvol, 2), 3);
           
       case 'coronal'
           subvol = abs(vdat(:, whsl, :));
           sumabsval = sum(sum(subvol, 1), 3);
           
       otherwise
           error('unknown slice view.')
   end

   wh_omit = sumabsval < 100 * eps;
   whsl(wh_omit) = []; 

end

% Eliminate empty areas of image (skip dim we're varying along)
% ------------------------------------------------------------------------
if dotight
    vdat = eliminate_empty_areas(vdat, myview);
end

% stack images to display together
% ------------------------------------------------------------------------

% starting and ending values for slices
% indices into whsl
nrows = ceil(length(whsl) ./ s);

for j = 1:nrows
    rowslices{j} = whsl((j-1)*s+1 : min(j*s, length(whsl)));
end

plate = cell(nrows, 1);

for j = 1:nrows
    stacked = [];
    
    for i = rowslices{j}  %st(j):en(j)
        
        switch myview
            case 'axial'
                stacked{i} = vdat(:, :, i);
                stacked{i} = rot90(stacked{i});
            case {'saggital', 'sagittal'}
                stacked{i} = squeeze(vdat(i, :, :));
                stacked{i} = rot90(stacked{i});
            case 'coronal'
                stacked{i} = squeeze(vdat(:, i, :));
                stacked{i} = rot90(stacked{i});
            otherwise error('unknown slice view.')
        end
        
    end
    
    stacked = cat(stackorient, stacked{:});
    
    plate{j} = stacked;
    
end % rows/plates

maxrows = max(cellfun(@(x) size(x, 1), plate, 'UniformOutput', true));
maxcols = max(cellfun(@(x) size(x, 2), plate, 'UniformOutput', true));
padval = 0; % mean(plate{j}(:)); % value to use for pad

% Pad if needed
for j = 1:nrows
    if size(plate{j}, 2) < maxcols
        % pad
        plate{j} = [plate{j} padval .* ones(maxrows, maxcols - size(plate{j}, 2))];
    end
    
end % rows

plate = cat(platestack, plate{:});

% create image
% ------------------------------------------------------------------------
if do_plot
    ax_h = gca;
    han = imagesc(plate);

    axis image
    set(ax_h, 'YDir', 'Reverse')
    axis off

    apply_split_colormap_and_clim(ax_h, plate, clim);

    % Only take over the figure when we own it (single axes in the figure).
    fh_here = ancestor(ax_h, 'figure');
    if numel(findall(fh_here, 'Type', 'axes')) <= 1
        fill_axes_in_figure(ax_h, plate);
    end
end

end % Main function

% function vdat = eliminate_empty_areas(vdat)
% % Eliminate empty areas of image
% % ------------------------------------------------------------------------
% % this will mess up conversion to voxel coords from mm if empty slices!
% % turn off in direction of variation
% 
%     nullvox = vdat == 0 | isnan(vdat);
%     bottom = all(nullvox, 3);
%     
%         nullx = squeeze(all(bottom, 2));
%         vdat(nullx, :, :) = [];
%     
%         nully = squeeze(all(bottom, 1));
%         vdat(:, nully, :) = [];
%     
%         side = squeeze(all(nullvox, 1));
%         nullz = squeeze(all(side, 1));
%         vdat(:, :, nullz) = [];
%         
% end

function plate = compose_three_views(ax_plate, cor_plate, sag_plate)
% Stack three horizontal montage strips vertically:
%
%   +----------------------------+
%   |       Axial montage        |
%   +----------------------------+
%   |       Coronal montage      |
%   +----------------------------+
%   |       Sagittal montage     |
%   +----------------------------+
%
% Pad each strip on the right with zeros so they share a common width.
% Zero maps to white in canlab_hot_cool_colormap, so the padding is
% invisible against the white figure background.

target_w  = max([size(ax_plate, 2), size(cor_plate, 2), size(sag_plate, 2)]);
ax_plate  = pad_right(ax_plate,  target_w);
cor_plate = pad_right(cor_plate, target_w);
sag_plate = pad_right(sag_plate, target_w);

plate = [ax_plate; cor_plate; sag_plate];
end % compose_three_views


function M = pad_right(M, target_w)
[~, w] = size(M);
if w < target_w
    M = [M, zeros(size(M, 1), target_w - w)];
end
end % pad_right


function M = pad_bottom(M, target_h)
[h, ~] = size(M);
if h < target_h
    M = [M; zeros(target_h - h, size(M, 2))];
end
end % pad_bottom


function apply_split_colormap_and_clim(ax_h, data, clim_override)
% Apply canlab_hot_cool_colormap to ax_h with a symmetric CLim so that
% the value 0 maps to the white center of the colormap.

if nargin < 3, clim_override = []; end

if ~isempty(clim_override) && numel(clim_override) == 2
    clim_use = clim_override;
else
    v = data(:);
    v = v(~isnan(v) & v ~= 0);
    if isempty(v)
        m = 1;
    else
        m = max(abs(v));
        if m == 0, m = 1; end
    end
    clim_use = [-m m];
end

set(ax_h, 'CLim', clim_use);
colormap(ax_h, canlab_hot_cool_colormap(256));
set(ancestor(ax_h, 'figure'), 'Color', 'w');
end % apply_split_colormap_and_clim


function fill_axes_in_figure(ax_h, plate)
% Resize axes to fill the figure and (when the figure is newly created)
% match the figure to the plate aspect ratio so axis image leaves
% essentially no whitespace.

set(ax_h, 'Position', [0 0 1 1], 'Units', 'normalized');

fh = ancestor(ax_h, 'figure');
plate_h = size(plate, 1);
plate_w = size(plate, 2);
if plate_h <= 0 || plate_w <= 0, return; end
aspect = plate_w / plate_h;

set(fh, 'Units', 'pixels');
pos = get(fh, 'Position');
new_h = pos(4);
new_w = round(new_h * aspect);
screen = get(0, 'ScreenSize');
max_w = round(screen(3) * 0.9);
max_h = round(screen(4) * 0.9);
if new_w > max_w
    new_w = max_w;
    new_h = round(new_w / aspect);
end
if new_h > max_h
    new_h = max_h;
    new_w = round(new_h * aspect);
end
set(fh, 'Position', [pos(1) pos(2) new_w new_h]);
end % fill_axes_in_figure


function [wh_slice, slice_cell] = display_slices_multi_image(dat, myview, slice_mm, ...
    nimages_req, image_names, clim)
% Render one slice per image in dat as a tiled subplot figure.
% For statistic_image, draw a canonical T1 underlay with hot/cool blobs
% over the significant voxels; otherwise imagesc the slice directly.

nimgs_total = size(dat.dat, 2);
if isempty(nimages_req)
    nimgs = nimgs_total;
else
    nimgs = min(nimages_req, nimgs_total);
end
if nimgs > 64
    warning('display_slices:tooManyImages', ...
        'multi_image: showing only the first 64 of %d images', nimgs);
    nimgs = 64;
end

if isempty(image_names)
    image_names = arrayfun(@(k) sprintf('Img %d', k), 1:nimgs, 'UniformOutput', false);
end

switch myview
    case 'axial', wh_col = 3;
    case {'saggital', 'sagittal'}, wh_col = 1;
    case 'coronal', wh_col = 2;
    otherwise, error('display_slices:badView', 'unknown slice view: %s', myview);
end

xyzmm = zeros(1, 3); xyzmm(wh_col) = slice_mm;
xyzvox = mm2voxel(xyzmm', dat.volInfo.mat, 1);
wh_slice = xyzvox(wh_col);

do_flip = sign(dat.volInfo.mat(1)) < 0;

sz = dat.volInfo.dim(1:3);
% mm2voxel returns the voxel index in stored order; flipdim along dim 1
% remaps voxel k -> sz(1)-k+1, so adjust the sagittal slice accordingly.
if do_flip && wh_col == 1
    wh_slice = sz(1) - wh_slice + 1;
end
wh_slice = max(1, min(sz(wh_col), wh_slice));

is_stat = isa(dat, 'statistic_image');

u_slice = [];
if is_stat
    underlay_file = canlab_get_underlay_image;
    if isempty(underlay_file) || ~exist(underlay_file, 'file')
        warning('display_slices:noUnderlay', ...
            'Canonical underlay not found; falling back to raw slice display.');
        is_stat = false;
    else
        u_obj = fmri_data(underlay_file, 'noverbose');
        u_obj = resample_space(u_obj, dat);
        u_vol = reconstruct_image(u_obj);
        if do_flip, u_vol = flipdim(u_vol, 1); end
        u_slice = extract_oriented_slice(u_vol, myview, wh_slice);
    end
end

slice_cell = cell(nimgs, 1);
sig_cell   = cell(nimgs, 1);
all_max    = 0;

for i = 1:nimgs
    sng = get_wh_image(dat, i);
    vi  = reconstruct_image(sng);
    if do_flip, vi = flipdim(vi, 1); end
    sl  = extract_oriented_slice(vi, myview, wh_slice);

    sig_sl = true(size(sl));
    if is_stat && isprop(sng, 'sig') && ~isempty(sng.sig)
        sig_full = sng;
        sig_full.dat = double(sng.sig(:, 1));
        sm = reconstruct_image(sig_full);
        if do_flip, sm = flipdim(sm, 1); end
        sig_sl = extract_oriented_slice(sm, myview, wh_slice) ~= 0;
    end

    slice_cell{i} = sl;
    sig_cell{i}   = sig_sl;
    vals = sl(sig_sl);
    if ~isempty(vals), all_max = max(all_max, max(abs(vals))); end
end

if isempty(clim) || numel(clim) ~= 2
    if all_max == 0, all_max = 1; end
    clim_used = [-all_max all_max];
else
    clim_used = clim;
end

ncol = ceil(sqrt(nimgs));
nrow = ceil(nimgs / ncol);

fh = create_figure('display_slices_multi');
set(fh, 'Color', 'w');
clf(fh);

% Delete any leftover axes from create_figure / prior renders.
delete(findall(fh, 'Type', 'axes'));

pad     = 0.005;                 % gap between cells (figure-normalized)
title_h = 0.030;                 % fraction reserved for each subplot title
cell_w  = (1 - (ncol + 1) * pad) / ncol;
axes_h  = (1 - (nrow + 1) * pad - nrow * title_h) / nrow;
if axes_h <= 0
    axes_h  = (1 - (nrow + 1) * pad) / nrow;
    title_h = 0;
end

cmap_for_blobs = canlab_hot_cool_colormap(256);

for i = 1:nimgs
    r = ceil(i / ncol);
    c = i - (r - 1) * ncol;
    left          = pad + (c - 1) * (cell_w + pad);
    v_top_r       = 1 - (r - 1) * (title_h + axes_h + pad) - pad;
    axes_bottom_r = v_top_r - title_h - axes_h;

    ax = axes('Parent', fh, 'Units', 'normalized', ...
        'Position', [left, axes_bottom_r, cell_w, axes_h]);

    sl = slice_cell{i};

    if is_stat
        blob = sl;
        blob(~sig_cell{i}) = NaN;
        render_underlay_with_blobs(ax, u_slice, blob, clim_used, cmap_for_blobs);
    else
        imagesc(ax, sl);
        axis(ax, 'image'); set(ax, 'YDir', 'Reverse'); axis(ax, 'off');
        set(ax, 'CLim', clim_used);
        colormap(ax, cmap_for_blobs);
    end

    if title_h > 0
        title(ax, image_names{i}, 'FontSize', 9, 'Interpreter', 'none');
    end
end

% Exactly nimgs axes were created -> no leftover empty subplot frames.

end % display_slices_multi_image


function sl = extract_oriented_slice(vol, myview, wh_slice)
sz = size(vol);
switch myview
    case 'axial'
        wh_slice = max(1, min(sz(3), wh_slice));
        sl = rot90(squeeze(vol(:, :, wh_slice)));
    case {'saggital', 'sagittal'}
        wh_slice = max(1, min(sz(1), wh_slice));
        sl = rot90(squeeze(vol(wh_slice, :, :)));
    case 'coronal'
        wh_slice = max(1, min(sz(2), wh_slice));
        sl = rot90(squeeze(vol(:, wh_slice, :)));
    otherwise
        error('display_slices:badView', 'unknown slice view: %s', myview);
end
end % extract_oriented_slice


function render_underlay_with_blobs(ax_h, u_slice, blob_slice, clim, cmap)
% Render a grayscale T1 underlay with the blob slice composited on top in
% truecolor RGB. Blob colors are sampled from cmap (canlab_hot_cool_colormap)
% using the symmetric CLim so the same value-to-color mapping applies
% across all subplots.

if nargin < 5 || isempty(cmap)
    cmap = canlab_hot_cool_colormap(256);
end

% Preserve the caller's Position across cla('reset')
saved_units    = get(ax_h, 'Units');
saved_position = get(ax_h, 'Position');
cla(ax_h, 'reset');
set(ax_h, 'Units', saved_units, 'Position', saved_position);

u = u_slice;
u(isnan(u)) = 0;
lo = min(u(:)); hi = max(u(:));
if hi <= lo, hi = lo + 1; end
u_rgb = repmat((u - lo) / (hi - lo), [1 1 3]);
image(ax_h, u_rgb);
axis(ax_h, 'image'); set(ax_h, 'YDir', 'reverse'); axis(ax_h, 'off');
hold(ax_h, 'on');

mask = ~isnan(blob_slice) & blob_slice ~= 0;
if ~any(mask(:)), hold(ax_h, 'off'); return; end

if isempty(clim) || numel(clim) ~= 2
    m = max(abs(blob_slice(mask)));
    if m == 0, m = 1; end
    clim = [-m m];
end

ncolors = size(cmap, 1);
m_abs = max(abs(clim));
if m_abs == 0, m_abs = 1; end

[H, W] = size(blob_slice);
blob_rgb = zeros(H, W, 3);

v = blob_slice(mask);
% Map v in [-m_abs, m_abs] to colormap row 1..ncolors
idx = round((v + m_abs) ./ (2 * m_abs) .* (ncolors - 1)) + 1;
idx = max(1, min(ncolors, idx));

for ch = 1:3
    tmp = blob_rgb(:, :, ch);
    tmp(mask) = cmap(idx, ch);
    blob_rgb(:, :, ch) = tmp;
end

h_b = image(ax_h, blob_rgb);
set(h_b, 'AlphaData', double(mask));
hold(ax_h, 'off');
end % render_underlay_with_blobs


function vdat = eliminate_empty_areas(vdat, myview)
% Eliminate empty areas of image
% ------------------------------------------------------------------------
% this will mess up conversion to voxel coords from mm if empty slices!
% turn off in direction of variation

    nullvox = vdat == 0 | isnan(vdat);
    bottom = all(nullvox, 3);
    
    if strcmp(myview, 'saggital') || strcmp(myview, 'sagittal')
        % skip
    else
        nullx = squeeze(all(bottom, 2));
        vdat(nullx, :, :) = [];
    end
    
    if strcmp(myview, 'coronal')
        % skip
    else
        nully = squeeze(all(bottom, 1));
        vdat(:, nully, :) = [];
    end
    
    if strcmp(myview, 'axial')
        % skip
    else
        side = squeeze(all(nullvox, 1));
        nullz = squeeze(all(side, 1));
        vdat(:, :, nullz) = [];
    end
    
end
