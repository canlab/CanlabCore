function [slices_fig_h, slice_mm_coords, slice_vox_coords, newax] = cluster_orthviews_montage(spacing, myview, varargin)
% :Usage:
% ::
%
%    [slices_fig_h, slice_mm_coords, slice_vox_coords, axis_handles] = cluster_orthviews_montage(spacing, myview, [overlay], [other optional args])
%
% Runs on top of spm_orthviews, creates montages from current orthviews
% display, whatever it is
%
% :Examples:
% ::
%
%    cluster_orthviews_montage(6, 'coronal');   % 6 mm spacing
%    cluster_orthviews_montage(10, 'sagittal', 'range', [-10 10]); % 10 mm spacing sag view with only parasagittal slices
%    cluster_orthviews_montage(12, 'axial');    % 12 mm spacing, axial view
%
% additional options: enter AFTER overlay:
%   - 'whichorth', whichorth = varargin{i+1}; varargin{i:i+1} = [];
%   - 'onerow', doonerow = 1; varargin{i} = [];
%   - 'range', followed by [min max] mm coords for slices
%   - 'xhairs', xhairs = 1; turn on cross-hairs on slice plot
%
% used in cluster_orthviews_classes
%
% ..
%    tor wager, aug 2006; updated (minor) April 2011. Update: Aug 2012 -
%    changed default slices to match canlab_results_fmridisplay
% ..


overlay = [];
xhairs = 0;
doonerow = 0;
whichorth = [];

if length(varargin) > 0
    overlay = varargin{1}; varargin{1} = [];
end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'whichorth', whichorth = varargin{i+1}; varargin{i:i+1} = [];
            case 'onerow', doonerow = 1; varargin{i} = [];
                
            case 'range'
                srange = varargin{i+1}; 
                if length(srange) ~= 2, error('Range must be [min max] mm coords'); end
                
            case 'xhairs', xhairs = 1;
                
            otherwise, warning('scn_tools:badInput', ['Unknown input string option:' varargin{i}]);
        end
    end
end

if isempty(overlay), overlay = which('SPM8_colin27T1_seg.img'); end %which('scalped_single_subj_T1.img'); end

% set up orthviews to be OK
try % will crash if no blobs
    scn_export_spm_window('setup', overlay);
catch
    disp('Error setting up orthviews window. Running with no blobs?');
end
%   spm_orthviews_white_background

%     set(get(findobj('Tag','Graphics'), 'Children'), 'Color', 'white');
set(findobj('Type', 'axes', 'Parent', findobj('Tag','Graphics')), 'Color', 'black');
if xhairs, spm_orthviews('Xhairs', 'on'); end

% Set window to be equal/constant
[volInfo, dat] = iimg_read_img(overlay, 2);
dat = dat(volInfo.wh_inmask);
spm_orthviews('Window', whichorth, [prctile(dat, 2) prctile(dat, 98)]);

myviews = {'axial' 'coronal' 'sagittal'};   % for selecting SPM window
whview = find(strcmp(myviews, myview));
if isempty(whview), error('myview must be axial, coronal, or sagittal.'); end


switch whview
    case 1 % axial
        if ~exist('srange', 'var'), srange = [-40 50]; end
        
        cen = [srange(1):spacing:srange(2)]';
        slice_mm_coords = cen;
        cen = [zeros(length(cen), 2) cen];
        xyzvox = mm2voxel(cen, volInfo.mat);
        slice_vox_coords = xyzvox(:, 3);
    case 2
        if ~exist('srange', 'var'), srange = [-90 55]; end
        
        cen = [srange(1):spacing:srange(2)]';
        slice_mm_coords = cen;
        cen = [zeros(length(cen),1) cen zeros(length(cen),1)];
        xyzvox = mm2voxel(cen, volInfo.mat);
        slice_vox_coords = xyzvox(:, 2);
    case 3
        if ~exist('srange', 'var'), srange = [-20 20]; end
        
        cen = [srange(1):spacing:srange(2)]';
        slice_mm_coords = cen;
        cen = [ cen zeros(length(cen), 2)];
        xyzvox = mm2voxel(cen, volInfo.mat);
        slice_vox_coords = xyzvox(:, 1);
end




myviews2 = {'sagittal' 'coronal' 'axial' };  % for selecting coord
whcoord = strmatch(myview, myviews2) ;

% get text string base
mystr = {'x = ' 'y = ' 'z = '};
textbase = mystr{whcoord};

% get optimal number of axes
num_axes = size(cen, 1);
rc = ceil(sqrt(num_axes));

if isempty(whichorth)
    axh = get_orth_axishandles;
else
    axh = get_orth_axishandles_whichorth(whichorth);
end

axh = axh(whview);

% Set up figure

slices_fig_h = create_figure(['montage_' myview]);
set(slices_fig_h, 'Color', 'k');
if not(feature('ShowFigureWindows'));
  set(slices_fig_h, 'Visible', 'off');
end

% copy existing colormap
fh = findobj('Tag', 'Graphics');
if strcmp(get(fh, 'Type'), 'figure') && ishandle(fh)
    set(slices_fig_h, 'colormap', get(fh, 'colormap'))
end

				% Resize figure based on view
if feature('ShowFigureWindows');
  ss = get(0, 'ScreenSize');
  if doonerow
    switch myview
      case {'axial'}
        set(slices_fig_h, 'Position', [round(ss(3)/12) round(ss(4)*.9) round(ss(3)*.9) round(ss(4)/7) ])
      case {'coronal'}
        set(slices_fig_h, 'Position', [round(ss(3)/12) round(ss(4)*.5) round(ss(3)*.9) round(ss(4)/7) ])
      case 'sagittal'
        set(slices_fig_h, 'Position', [round(ss(3)/12) round(ss(4)*.7) round(ss(3)*.6) round(ss(4)/5.5) ])
    end
  else
    switch myview
      case {'axial'}
		%how far right, how far up, how big across, how big up
        set(slices_fig_h, 'Position', [round(ss(3)/12) round(ss(4)/12) round(ss(3)*.7) round(ss(4)*.7) ])
      case {'coronal'}
        set(slices_fig_h, 'Position', [round(ss(3)/12) round(ss(4)/12) round(ss(3)*.7) round(ss(4)*.7) ])
      case {'sagittal'}
        set(slices_fig_h, 'Position', [round(ss(3)/12) round(ss(4)/12) round(ss(3)*.7) round(ss(4)*.7) ])
    end
  end
end


for i = 1:num_axes
    
    if doonerow
        newax(i) = subplot(1, num_axes, i);
    else
        newax(i) = subplot(rc, rc, i);
    end
    
    axis off;
end


for i = 1:num_axes
    spm_orthviews('Reposition', cen(i, :));
    
    spm_orthviews_showposition;
    
    if ~ishandle(axh)
        disp('SPM figure was deleted or does not exist, possibly because a copied montage figure was closed? Skipping');
        continue
    end
    
    copyobj(get(axh, 'Children'), newax(i));
    axes(newax(i));
    axis image
    
    ehan = get(newax(i), 'Children');

    %DeleteFcn: This prevents deleting windows from clearing SPM orthviews
    %windows
    set(ehan, 'DeleteFcn', []);
    
    % Transparent background
    wh = strcmp('image', get(ehan, 'Type')); % | strcmp('line', get(ehan, 'Type'));
    ehan = ehan(wh);
    
  % set transparent value for clear axes
    myAlphaData = ~(all(get(ehan, 'CData') == 0, 3) | (all(get(ehan, 'CData') == 1, 3)));

        % If we set alphadata to clear for BG and axis color to none, we get clear
        % axes
   set(ehan, 'AlphaDataMapping', 'scaled', 'AlphaData', myAlphaData)
               
    
    % get rid of green text
    h = findobj(newax(i), 'Type', 'text');
    delete(h);
    
    h = findobj(newax(i), 'Type', 'text');
    if length(h) > 1, h= h(1); end % kludgy fix for multiple matches ***
    
    % try to set a reasonable font size
    pos = get(newax(i), 'Position');
    height = pos(3);
    fs = max(10, round(height * 70));
    set(h, 'FontSize', fs)
    % set position of text (move down and right)
    pos = get(h, 'Position');
    pos(2) = 0; %pos(2) - fs./2;
    pos(1) = 0;
    set(h, 'Position', pos)
    % set text string based on cen
    %set(h, 'String', [textbase num2str(cen(i))]);
    
    %         elseif ishandle(h)
    %             delete(h);
    %         end
end

% adjust axes
h = findobj(slices_fig_h, 'Type', 'axes');
pos = get(h, 'Position');
pos = cat(1, pos{:});
% pos2 = get(h, 'OuterPosition');
% pos2 = cat(1, pos2{:});

% % left bottom justify
%pos(:, 1) = pos(:, 1) - min(pos(:, 1)) + .03;


if doonerow
    % respace based on figure: X start
    pos(:, 1) = linspace((.98-(1./length(h))), .02, length(h));
    % X width
    pos(:, 3) = 1./length(h);
    % % left bottom justify
    pos(:, 2) = pos(:, 2) - min(pos(:, 2)) + .05;
else
   newrows = [1; 1+find(diff(pos(:,1))>0)]; 
   spacing = newrows(2)-newrows(1);
   for i = 1:numel(newrows)
       pos(newrows(i):spacing+newrows(i)-1,1) = linspace((.98-(1./spacing)), .02, spacing);
       pos(newrows(i):spacing+newrows(i)-1,2) = .05+(i-1)*(.75/(numel(newrows)-1))*ones(spacing,1);
   end
   pos(:,3) = 1./spacing;
   switch myview
       case {'coronal'}
           pos(:,4) = pos(:,4)*1.2;
   end
end

for i = 1:length(h)
    set(h(i), 'OuterPosition', pos(i, :));
    set(h(i), 'Position', pos(i, :));
end

% If we set alphadata to clear for BG and axis color to none, we get clear axes
set(h, 'Color', 'none')
               
% try to set a reasonable enlargement factor
% n = 1 + .2 * log(num_axes ./ 2);
% n = max(n, 1); n = min(n, 3);
% enlarge_axes(gcf, 1, n)

% set background color to print as whatever it is on screen
set(gcf,'InvertHardcopy', 'off');

h = findobj(slices_fig_h, 'Type', 'text');
set(h, 'FontSize', 24)

end




function axish = get_orth_axishandles

% Get figure and axis handles
fh = findobj('Tag', 'Graphics');
ch = get(fh, 'Children');
for i= 1:length(ch)
    mytype = get(ch(i), 'Type');
    wh(i) = strcmp(mytype, 'axes');
end
axish = ch(find(wh));

if isempty(axish)
    disp('SPM figure orthviews do not exist'); 
    disp('You must clear figures for orthviews with previous copies of orthviews objects')
    disp('Before setting up spm orthviews.')
    disp('Clearing existing figures...try re-running now.')
    
    % clear existing windows
    figh = findobj('Tag', 'montage_axial');
    if ishandle(figh), clf(figh); end
    
    figh = findobj('Tag', 'montage_coronal');
    if ishandle(figh), clf(figh); end
    
    figh = findobj('Tag', 'montage_sagittal');
    if ishandle(figh), clf(figh); end
    
end

% get which axis is which
for i = 1:length(axish)
    poss(i, :) = get(axish(i), 'Position');
end

% get rid of extra axes we may have created in the 4th quadrant
other_axes = find(any(poss(:, 1) > .45 & poss(:, 2) < .2, 2));
axish(other_axes) = [];
poss(other_axes, :) = [];

% sort into order:  axial, coronal, saggital
ssum = sum(poss(:, 1:2), 2);
[ssum, ind] = sort(ssum);
axish = axish(ind);

end

function axish = get_orth_axishandles_whichorth(whichorth)
% Get the axis handles for the current orthviews
global st
for i = 1:length(st.vols), wh(i) = ~isempty(st.vols{i}); end
wh = find(wh); wh = wh(whichorth);
axish = cat(1, st.vols{wh}.ax{:});
axish = sort(cat(1, axish(:).ax));

end

