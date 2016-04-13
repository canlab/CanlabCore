function draw_anatomical_roi_2008(meth, varargin)
% Quick start guide:
% type draw_anatomical_roi_2008
%
% :Usage:
% ::
%
%     draw_anatomical_roi_2008('init');
%     draw_anatomical_roi_2008('load', 'ROI_midbrain.img');
%     draw_anatomical_roi_2008('init', 'overlay',  'remi_mean_T2.img');
%
% :Inputs:
%
%   no arguments: init
%
%   **'init':**
%        initialize gui and orthviews and remove previous ROI
%
%   **'load':**
%        load an ROI from a mask
%
%   **'free':**
%        draw an ROI freehand (click on one of the slices in the Slices window 1st)
%
%   **'poly':**
%        don't run this yet
%
%   **'add':**
%        add a region you've drawn on a slice to your ROI 
%
%   **'remove':**
%        remove a region you've drawn from your ROI
%
%   **'smooth':**
%        3-d smoothing of ROI
%
%   **'write':**
%        write mask image of ROI and return clusters to workspace
%
%   **'exit':**
%        exit. ROI data is stored in the Slices figure, so you can
%        continue to edit, etc. after exiting.
%
% :Note:
% You can use cluster_orthviews to image multiple blobs, and then
% draw relative to those.
% this function saves it's data in the Slices figure, so you can draw,
% re-initialize the orthviews, and keep drawing before you save.
%
% :Examples:
% ::
%
%    cluster_orthviews(red, {[1 0 0]}, 'overlay', 'remi_mean_T2.img');
%    cluster_orthviews(stn, {[0 1 0]}, 'add');
%    set(findobj('Tag','Graphics'), 'WindowButtonUpFcn', 'draw_anatomical_roi_2008(''moveslice'');');
%    % Use the spm_orthviews menu to ZOOM IN...and keep drawing!
%
% Example of brainstem ROI drawing:
% ::
%
%    draw_anatomical_roi_2008('init', 'overlay',  'remi_mean_T2.img');
%    set(findobj('Tag','Graphics'), 'WindowButtonUpFcn', '');
%    cluster_orthviews(red, {[1 0 0]}, 'overlay', 'remi_mean_T2.img');
%    stn = mask2clusters('ROI_STN.img');
%    cluster_orthviews(stn, {[0 1 1]}, 'add');
%    % Now zoom in to the midbrain in SPM orthviews and draw new ROIs
%
% ..
%    Tor Wager, Dec 2008
% ..


if nargin < 1, meth = 'init'; end

while ~strcmp(meth, 'exit')

%     switch meth
%         % for drawing tools only, get tag based on which axis you clicked on
%         % update orientation data.orient
%         case {'free', 'poly'}
%             get_current_orientation;
%     end


    switch meth

        case 'choose_next'
            meth = input('Enter command(init/load/free/poly/add/remove/smooth/write/exit): ', 's'); % 'free';

        case 'init'
            init_slice_figure(varargin{:});
            move_slice;
            meth = 'exit';

        case 'load'
            load_vol(varargin{:});
            meth = 'exit';
  
        case 'moveslice', move_slice;
            % This is primarily used as the mouse-up function in
            % spm_orthviews, and is run when you click the orthviews
            meth = 'exit';

        case 'free'
            disp('Draw freehand region on the axial slice in Slices fig.');
            draw_freehand;
            meth = 'exit';
            
        case 'freesagg'
            disp('Draw freehand region on the saggital slice in Slices fig.');
            data = get_gui_data;
            axis_handle = data.saggh;
            draw_freehand(axis_handle);
            meth = 'exit';

        case 'zoomin', zoom_in;
            meth = 'exit';
            
        case 'zoomout', zoom_out;
            meth = 'exit';
            
        case 'poly'
            draw_poly;
            meth = 'exit';

        case 'add'
            add_region_to_roi_mask;
            meth = 'exit';

        case 'remove'
            add_region_to_roi_mask('remove');
            meth = 'exit';
            
        case 'smooth'
            smooth_vol(varargin{:}); 
            meth = 'exit';
            
        case 'write'
            write_roi_mask;
            meth = 'exit';

        case 'exit'
            % quit out
            return

        otherwise
            disp('unknown method.')
            meth = 'choose_next';

    end


end % while loop

end

% --------------------------------
% Get data from GUI figure
% --------------------------------

function [data, f] = get_gui_data
f = findobj('Type', 'Figure', 'Tag', 'Slices');
data = guidata(f);
end

% --------------------------------
% Update orientation of figure in just-clicked-on axis
% --------------------------------

function get_current_orientation

[data, f] = get_gui_data;


mytag = get(gca, 'Tag');
switch mytag
    case 'ax', orient = 'ax';
    case 'sagg', orient = 'sagg';
    otherwise
        disp('Click on one of the slices in the Slices window.');
        orient = [];
end

data.orient = orient;
guidata(f, data);

end

% --------------------------------
% Initialize orthviews and slice figure (which stores data)
% --------------------------------

function init_slice_figure(varargin)
disp('Initializing orthviews and slices.');
disp('This clears any previous ROI you''ve been building.')
disp('Use free or poly commands to draw regions on slices,')
disp('and then ''add'' to add them to your ROI if you are');
disp('satisfied with them.')

ovl = which('scalped_avg152T1.img');

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'ovl', 'overlay'}, ovl = varargin{i+1}; varargin{i+1} = [];
            otherwise, warning('drawRoi:unknownoption', ['Unknown input string option:' varargin{i}]);
        end
    end
end

V = spm_vol(ovl);
dat = spm_read_vols(V);
dim = V.dim;

spm_image('init', ovl);
orthfig = findobj('Type', 'Figure', 'Tag', 'Graphics');
set(orthfig, 'WindowButtonUpFcn', 'draw_anatomical_roi_2008(''moveslice'');');

f = create_figure('Slices', 2, 2);
%data = guidata(f);
data = [];
data.ovl = ovl;
data.V = V;
data.dat = dat;
data.dim = dim;

subplot(2, 2, 1); hold off;
data.axh = gca;
set(gca, 'Tag', 'ax'); % Tag is used for getting which axis you clicked; must be re-set every time hold is off
subplot(2, 2, 2); hold off;
data.saggh = gca;
set(gca, 'Tag', 'sagg');

subplot(2, 2, 3); hold off;
% data.corh = gca;
% set(gca, 'Tag', 'cor');

create_uibuttons;

create_figure('Volume');
data.surfh = gca;

Rvol = zeros(data.dim);
data.Rvol = Rvol;
guidata(f, data);
end

% --------------------------------
% Load from mask .img file
% --------------------------------
function load_vol(varargin)

if nargin < 1, img = spm_select(1);
else img = varargin{1};
end

[data, f] = get_gui_data;

dat = scn_map_image(img, data.V.fname);
dat = ~isnan(dat) & dat ~= 0;
data.Rvol = dat;

guidata(f, data);

move_slice;
draw_volume_render;
end

% --------------------------------
% Redraw slices based on current orthviews position
% --------------------------------

function move_slice
%disp('Showing slices based on current coordinate position.');
coord = spm_orthviews('Pos');

[data, f] = get_gui_data;
data.coord = coord;

vox = mm2voxel(coord, data.V.mat);
data.wh_slice_ax = round(vox(3));
data.wh_slice_sagg = round(vox(1));

figure(f)

if isempty(data.axh) || isempty(data.saggh) || ~ishandle(data.axh) || ~ishandle(data.saggh)
    disp('Slices window was closed. Re-''init''');
    return
end

% draw slices
axslice = data.dat(:, :, data.wh_slice_ax);
saggslice = rot90(squeeze(data.dat(data.wh_slice_sagg, :, :)));

axes(data.axh); hold off;
imagesc(axslice);
axis tight, axis image, hold on

axes(data.saggh); hold off;
imagesc(saggslice);
axis tight, axis image, hold on

% set axis limits if we have them from before
% (we don't have them if we're initializing)
if isfield(data, 'axislimits') && ~isempty(data.axislimits)
    set(data.axh, 'XLim', data.axislimits.axx);
    set(data.axh, 'YLim', data.axislimits.axy);
    set(data.saggh, 'XLim', data.axislimits.saggx);
    set(data.saggh, 'YLim', data.axislimits.saggy);
    
    axslice = axslice(ceil(data.axislimits.axy(1)) : floor(data.axislimits.axy(2)), ...
        ceil(data.axislimits.axx(1)) : floor(data.axislimits.axx(2)));

    saggslice = saggslice(ceil(data.axislimits.saggy(1)) : floor(data.axislimits.saggy(2)), ...
        ceil(data.axislimits.saggx(1)) : floor(data.axislimits.saggx(2)));
end

% set color limits (for good contrast)
stdax = std(axslice(:));
stdsagg = std(saggslice(:));
meanax = mean(axslice(:));
meansagg =  mean(saggslice(:));
axclim = [meanax - 2*stdax meanax + 2*stdax];
saggclim = [meansagg - 2*stdsagg meansagg + 2*stdsagg];
%axclim = [prctile(axslice(:), 5) prctile(axslice(:), 95)];
%saggclim = [prctile(saggslice(:), 5) prctile(saggslice(:), 95)];
set(data.axh, 'CLim', axclim)
set(data.saggh, 'CLim', saggclim)



set(data.axh, 'Tag', 'ax'); % Tag is used for getting which axis you clicked; must be re-set every time hold is off
set(data.saggh, 'Tag', 'sagg');

% Add previous regions
myslice = zeros(data.dim(1:2));
myslice(:, :, 2) = data.Rvol(:, :, data.wh_slice_ax);

% because of bug when zooming in, we want to make the space tight around
% regions

FVC = isocaps(myslice, 0, 'zmax');
axes(data.axh);
patch(FVC, 'EdgeColor', 'g', 'FaceColor', 'g');
%set(patchh, 'FaceAlpha', .5)

clear myslice, myslice(:, :, 1) = zeros(data.dim([3 2]));
myslice(:, :, 2) = rot90(squeeze(data.Rvol(data.wh_slice_sagg, :, :)));
FVC = isocaps(myslice, 0, 'zmax');
axes(data.saggh);
patch(FVC, 'EdgeColor', 'g', 'FaceColor', 'g'); %, %'FaceAlpha', .5);

axis tight, axis image, hold on, colormap gray

guidata(f, data);
end

% --------------------------------
% Region drawing tools
% --------------------------------

function draw_freehand(varargin)

data = get_gui_data;

if nargin > 0, axis_handle = varargin{1};
else axis_handle = data.axh;
end
axes(axis_handle);

get_current_orientation;
[data, f] = get_gui_data;

if isfield(data, 'fh') && ishandle(data.fh), delete(data.fh); end

try
    h = imfreehand(axis_handle);
catch
    disp('Error with freehand draw due to unknown axis-handling bug.'); 
    disp('You should probably close the Slices figure')
    disp('and re-initialize (''write'' first if you want).');
    return
end

api = iptgetapi(h);
position = api.getPosition();
xi = position(:, 1);
yi = position(:, 2);

data.fh = fill(xi, yi, 'r', 'FaceAlpha', .4);

if isempty(data.orient);  % bad axis
    disp('Click on one of the slices in the ''Slices'' figure and re-enter free command.')
    return
end

switch data.orient

    case 'ax'
        R = poly2mask(position(:, 1),position(:, 2),data.dim(1),data.dim(2));

    case 'sagg'
        R = poly2mask(position(:, 1),position(:, 2),data.dim(3),data.dim(2));
        R = rot90(R, -1);
end

data.R = R;

guidata(f, data);

% re-set window button
%set(f, 'WindowButtonUpFcn', 'draw_anatomical_roi_2008(''free'')');

end

% --------------------------------
% Zoom in and set explicit axis limits
% --------------------------------
function zoom_in

disp('Click on the top left and then bottom right')
disp('of the desired zoom area in the Slices window.');
pos = ginput(2);
get_current_orientation;

[data, f] = get_gui_data;

set(gca, 'XLim', sort(pos(1:2, 1)));
set(gca, 'YLim', sort(pos(1:2, 2)));

% get axis limits.
% We will use these to constrain view area
data.axislimits.axx = get(data.axh, 'XLim');
data.axislimits.axy = get(data.axh, 'YLim');
data.axislimits.saggx = get(data.saggh, 'XLim');
data.axislimits.saggy = get(data.saggh, 'YLim');

guidata(f, data);

end

% --------------------------------
% Zoom out: remove explicit axis limits
% --------------------------------
function zoom_out
[data, f] = get_gui_data;
if isfield(data, 'axislimits');
    data = rmfield(data, 'axislimits');
    guidata(f, data);
end

move_slice;
end

% --------------------------------
% Region drawing tools: poly
% --------------------------------
function draw_poly
disp('Drawing movable polygon region.');
disp('NOT COMPLETE OPTION YET.');

% [data, f] = get_gui_data;
% if isfield(data, 'fh') && ishandle(data.fh), delete(data.fh); end
% 
% if isempty(data.orient);  % bad axis
%     disp('Click on one of the slices in the ''Slices'' figure and re-enter free command.')
%     return
% end
% 
% % polygon
% [data.R, xi, yi] = roipoly;
% 
% data.fh = fill(xi, yi, 'r', 'FaceAlpha', .4);
% 
% guidata(f, data);
end

% -----------------------------------------
% Add (or remove) region to the ROI volume
% -----------------------------------------

function add_region_to_roi_mask(varargin)

doremove = 0; % remove region instead of adding it
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'remove', doremove = 1;
            otherwise, warning('drawRoi:unknownoption', ['Unknown input string option:' varargin{i}]);
        end
    end
end

if doremove
    disp('Removing region from ROI mask.');
    addstr = 'noadd';
else
    disp('Adding region to ROI mask.');
    addstr = 'add';
end

[data, f] = get_gui_data;
thickness = 3;
%data.R = repmat(data.R, [1 1 thickness]);

switch data.orient
    case 'ax'
        slice_vals = data.wh_slice_ax - (thickness - 1)/2 : data.wh_slice_ax + (thickness - 1)/2;
        Rslices = data.Rvol(:, :, slice_vals);
        if doremove
            oldR = data.R(:, :, slice_vals);
            Rslices(oldR == 1) = 0; 
        else
            Rslices(data.R == 1) = 1;
        end
        data.Rvol(:, :, slice_vals) = Rslices;
        %data.Rvol(:, :, slice_vals) = double(data.Rvol(:, :, slice_vals)) + data.R;

    case 'sagg'
        slice_vals = data.wh_slice_sagg - (thickness - 1)/2 : data.wh_slice_sagg + (thickness - 1)/2;
        
        % if new slice, permute and replicate to thickness.  if not, don't
        % if we click "add" on an empty slice, it will have old data in
        % it...
        if length(size(data.R)) ~= 3
            data.R = repmat(data.R, [1 1 thickness]);
            data.R = permute(data.R, [3 1 2]);      % Slices with data to replace
        end
        
        if doremove
            oldR = data.Rvol(slice_vals, :, :);
            oldR(logical(data.R)) = 0;
            data.Rvol(slice_vals, :, :) = oldR;
        
        else % add
            data.Rvol(slice_vals, :, :) = double(data.Rvol(slice_vals, :, :) | data.R);
        end
        
%         data.R = permute(data.R, [3 1 2]);      % Slices with data to replace
%         Rslices = data.Rvol(slice_vals, :, :);  % previous data to manipulate, from whole volume
%         if doremove
%             oldR = data.R(slice_vals, :, :);
%             Rslices(oldR == 1) = 0; 
%         else % add
%             oldR = permute(data.R, [3 1 2]);
%             Rslices(data.R == 1) = 1;           % preserve previous 1 values
%         end
%         data.Rvol(slice_vals, :, :) = Rslices;
%         %data.Rvol(slice_vals, :, :) =double(data.Rvol(slice_vals, :, :)) + data.R;

end

guidata(f, data);

% show it
move_slice;
draw_volume_render(addstr);

end

% --------------------------------
% Write ROI mask
% --------------------------------

function write_roi_mask
f = findobj('Type', 'Figure', 'Tag', 'Slices');
data = guidata(f);

data.Vout.fname = input('Enter file name (e.g., ROI.img): ', 's');
data.Vout.mat = data.V.mat;
data.Vout.dim = data.V.dim;
data.Vout.n = [1 1];
data.Vout = spm_create_vol(data.Vout);
spm_write_vol(data.Vout, data.Rvol);
disp(['Written: ' data.Vout.fname]);

cl_roi = mask2clusters(data.Vout.fname);

cluster_orthviews(cl_roi, {[1 0 0]}, 'overlay', data.ovl);

assignin('base', 'cl_roi', cl_roi)
disp('Assigned cluster cl_roi in base workspace.');

% set callback so we can keep drawing
orthfig = findobj('Type', 'Figure', 'Tag', 'Graphics');
set(orthfig, 'WindowButtonUpFcn', 'draw_anatomical_roi_2008(''moveslice'');');

end

% --------------------------------
% orthviews and vol rendering
% --------------------------------

function draw_volume_render(addstr)

if nargin < 1, addstr = 'add'; end

data = get_gui_data;

cl = mask2clusters(data.Rvol, data.V.mat);
if isempty(cl)
    spm_image('init', data.ovl);
else
    cluster_orthviews(cl, {[0 1 0]}, addstr, 'overlay', data.ovl);
end

% do this again, just in case it gets messed up
orthfig = findobj('Type', 'Figure', 'Tag', 'Graphics');
set(orthfig, 'WindowButtonUpFcn', 'draw_anatomical_roi_2008(''moveslice'');');

axes(data.surfh); hold on;
if isfield(data, 'surfobjh') && ishandle(data.surfobjh)
    delete(data.surfobjh)
else
    % draw new surface
    %                 brainh = addbrain; set(brainh, 'FaceAlpha', .2);
    view(132, 30);
    lighting gouraud; lightFollowView; material dull;
    axis image, axis tight
    rotate3d off
end

data.cl = cl;
if ~isempty(cl)
    data.surfobjh = imageCluster('cluster', cl, 'color', 'g');
    axis image; axis tight; lightRestoreSingle; lighting gouraud
    rotate3d off
else
    data.surfobjh = [];
end

end

function smooth_vol(varargin)

disp('Smoothing ROI region.');

[data, f] = get_gui_data;

data.Rvol = smooth3(data.Rvol, 'g', 5);
data.Rvol = data.Rvol > .5;

guidata(f, data);

move_slice;
draw_volume_render('noadd');

end


function create_uibuttons

[data, f] = get_gui_data;

x1start = 30;
x2start = 200;

str = 'draw_anatomical_roi_2008(''free'');';
str = expand_callback_str(str);

uicontrol(f,'String','Draw (axial)',...
'Position',[x1start 200-35*0 150 30],...
'CallBack', str,...
'Interruptible','on',...
'ForegroundColor','k','FontWeight','b');

str = 'draw_anatomical_roi_2008(''freesagg'');';
str = expand_callback_str(str);

uicontrol(f,'String','Draw (saggital)',...
'Position',[x2start 200-35*0 150 30],...
'CallBack', str,...
'Interruptible','on',...
'ForegroundColor','k','FontWeight','b');

% ------
str = 'draw_anatomical_roi_2008(''add'');';
str = expand_callback_str(str);

uicontrol(f,'String','Add region',...
'Position',[x1start 200-35*1 150 30],...
'CallBack', str,...
'Interruptible','on',...
'ForegroundColor','k','FontWeight','b');

% ------
str = 'draw_anatomical_roi_2008(''remove'');';
str = expand_callback_str(str);

uicontrol(f,'String','Remove region',...
'Position',[x1start 200-35*2 150 30],...
'CallBack', str,...
'Interruptible','on',...
'ForegroundColor','k','FontWeight','b');

% ------
str = 'draw_anatomical_roi_2008(''write'');';
str = expand_callback_str(str);

uicontrol(f,'String','Write ROI mask',...
'Position',[x1start 200-35*3 150 30],...
'CallBack', str,...
'Interruptible','on',...
'ForegroundColor','k','FontWeight','b');

% ------
str = 'draw_anatomical_roi_2008(''smooth'');';
str = expand_callback_str(str);

uicontrol(f,'String','Smooth 3D ROI',...
'Position',[x1start 200-35*4 150 30],...
'CallBack', str,...
'Interruptible','on',...
'ForegroundColor','k','FontWeight','b');

% ------
str = 'draw_anatomical_roi_2008(''zoomin'');';
str = expand_callback_str(str);

uicontrol(f,'String','Zoom in',...
'Position',[x2start 200-35*1 150 30],...
'CallBack', str,...
'Interruptible','on',...
'ForegroundColor','k','FontWeight','b');
axis off

str = 'draw_anatomical_roi_2008(''zoomout'');';
str = expand_callback_str(str);
uicontrol(f,'String','Zoom out',...
'Position',[x2start 200-35*2 150 30],...
'CallBack', str,...
'Interruptible','on',...
'ForegroundColor','k','FontWeight','b');

axis off
end

function str = expand_callback_str(str)

str2 = [];
for i = 1:length(str)
    if str(i) == '''', str2(end+1) = ''''; str2(end+1) = '''';
    else str2(end+1) = str(i);
    end
end

str = ['disp(''' char(str2) '''), ' str ];     % display then execute

end
