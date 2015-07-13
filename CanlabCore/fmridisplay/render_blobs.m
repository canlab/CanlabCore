function [blobhan, cmaprange, mincolor, maxcolor] = render_blobs(currentmap, mymontage, SPACE, varargin)

% This is a helper function for fmridisplay objects, called by the addblobs method
%
% Usage:
% -------------------------------------------------------------------------
% [blobhan, cmaprange, mincolor, maxcolor] = render_blobs(currentmap, mymontage, SPACE, varargin)
% See fmridisplay.m and addblobs.m method in fmridisplay for more details and options.
%
% Author and copyright information:
% -------------------------------------------------------------------------
%     Copyright (C) 2015 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Inputs:
% -------------------------------------------------------------------------
% currentmap : see addblobs method. Montage within fmridisplay object.
% mymontage : ditto
% SPACE: space of map to sample to (object display SPACE in fmridisplay object)
% 
% Optional inputs:
% There are many optional inputs that control display features of blobs.
% COLOR
% 'color',      followed by color vector, e.g., [0 1 1]
% 'maxcolor'    followed by color vector for max color range, e.g., [0 1 1]
% 'mincolor'    followed by color vector for min color range, e.g., [0 0 1]
% 'onecolor'    force solid-color blobs
% 'splitcolor'  Positive and negative values are mapped to different
%               colormaps. Default is +=hot, -=cool colors.  Followed
%               optionally by cell array with 
%               vectors of 4 colors defining max/min for +/- range, e.g., {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]}
%
% OUTLINING
% 'outline' 
% 'linewidth',  followed by width value, e.g., 1
%
% COLOR RANGE
% 'cmaprange',  followed by range of values, e.g., [0 40], [-3 3]. Used in
%               color and transparency setting under some circumstances.
% 
% TRANSPARENCY
% {'trans', 'transparent','scaledtransparency', 'constanttrans', [val], 'transvalue', [val]}
% 'trans'               Transparent blobs; with no other input, transparency = 0.75 (1 is opaque, 0 is transparent/invisible)
% 'scaledtransparency'  Transparency is a function of voxel value, lower values are more transparent
% 'transvalue'          Followed by width value, e.g., 1. also 'constanttrans'
%
% OTHER OPTIONS
%
% 'smooth'      Smooth blobs
% 'contour'
% Orientation: 'sagittal', 'coronal', 'axial'
%
% Outputs:
% -------------------------------------------------------------------------
% [blobhan, cmaprange, mincolor, maxcolor] 
% All used in addblobs.m
%
% Examples:
% -------------------------------------------------------------------------
%
% Use addblobs; do not run this function directly unless you are
% programming with it.
%
% See also:
% fmridisplay/addblobs, fmridisplay, fmridisplay/multi_threshold

% Programmers' notes:
% Matlab's graphics behavior changed in 2014, causing some problems.
% This version was tested by Tor on R2015a.


myview = mymontage.orientation;

dosmooth = 0;
color = [1 0 0];

% contour and outline options
docontour = 0;
dotrans = 0;
transvalue = []; % default: map with clim
contourmin = 100*eps;
outline = 0;
mylinewidth = 2;

% color-mapped blobs options
docolormap = 1;
mincolor = [0 0 1];
maxcolor = [1 1 0];
dosplitcolor = 0;
%cmaprange = [min(currentmap.mapdata(:)) max(currentmap.mapdata(:))];
mapd = currentmap.mapdata(:); mapd = mapd(mapd ~= 0 & ~isnan(mapd));
cmaprange = double([prctile(mapd, 10) prctile(mapd, 90)]);

% adjust defaults
if any(strcmp(varargin, 'splitcolor')) && ~any(strcmp(varargin, 'cmaprange'))
    cmaprange = double([prctile(mapd(mapd < 0), 25) prctile(mapd(mapd < 0), 75) prctile(mapd(mapd > 0), 25) prctile(mapd(mapd > 0), 75) ]);
end

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'smooth', dosmooth = 1;
                
            case 'color', docolormap = 0; color = varargin{i + 1}; % single color, turn off color mapping
                
            case {'maxcolor'}, color = varargin{i + 1}; % color for single-color solid or value-mapped
            case 'mincolor', mincolor = varargin{i + 1}; % minimum color for value-mapped colors
                
            case 'onecolor', docolormap = 0; % solid-color blobs
            case 'splitcolor'
                docolormap = 1; dosplitcolor = 1;
                splitcolors = varargin{i + 1};
                if ~iscell(splitcolors) || length(splitcolors) ~= 4
                    error('Enter splitcolor followed by a 4-element cell vector of 4 colors\n{Min_neg} {Max_neg} {Min_pos} {Max_pos}');
                end
                color = splitcolors{4}; % max pos
                mincolor = splitcolors{3}; % min pos
                minnegcolor = splitcolors{1}; % min neg
                maxnegcolor = splitcolors{2}; % max neg
                
                
                % do not use other entries for colors
                %                 varargin{strcmp(varargin, 'color')} = deal(0);
                %                 varargin{strcmp(varargin, 'maxcolor')} = deal(0);
                %                 varargin{strcmp(varargin, 'mincolor')} = deal(0);
                
            case 'cmaprange', cmaprange = double(varargin{i + 1}); % enter specific values mapped to min and max colors
            case {'trans', 'transparent'}
                dotrans = 1; % transparent blobs
                transvalue = 0.75; % Default behavior, superseded by later arguments
                
            case {'constanttrans', 'transvalue'}
                dotrans = 1; transvalue = double(varargin{i + 1}); 
                %  (default: map transparency to cmaprange, unless you
                %  enter contanttrans followed by a transparency value
                
            case 'scaledtransparency' % Transparency is a function of voxel value, lower values are more transparent
                transvalue = [];  % Empty invokes the scaled mapping later
                
                % contour options
            case 'contour', docontour = 1;
            case 'outline', docontour = 1; outline = 1;
            case 'linewidth', mylinewidth = varargin{i + 1};
                
                % orientation options
            case 'sagittal', myview = 'sagittal'; %disp('Warning! NOT implemented correctly yet!!!'), %pause(5)
            case 'coronal', myview = 'coronal'; %disp('Warning! NOT implemented correctly yet!!!'), pause(5)
            case 'axial', myview = 'axial';
                
            case 'wh_montages'
                % not functional, avoid warning
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end



blobhan = [];

isvalid = ~isempty(mymontage) && isfield(mymontage, 'axis_handles') && all(ishandle(mymontage.axis_handles));

if ~isvalid, return, end


handles = mymontage.axis_handles;
n = length(handles);

% for each slice...

% Summarize voxels shown and not shown

% find closest slice to each one in original 'blob' image to map
for j = 1:n
    
    switch myview % may be same as mymontage.orientation...
        case 'axial'
            
            my_slice_coord = mymontage.slice_mm_coords(j);
            
            [dummy, wh_slice(j)] = min(abs(currentmap.SPACE.zcoords - my_slice_coord));
            
            % Wani: The previous line gets the closest slice from the montages 
            % you're displaying by calculating the distance between them. 
            % However, if the distance measure is larger than a half of the 
            % voxel size (i.e., vox_size/2), it means there is no proper montage 
            % that can display the data. So I needed to add the following
            % line for each view. Actually, I tried to use V.mat(:,3) for
            % the axial view (because the voxel sizes for x,y,z could be
            % different), but for some reasons, it doesn't work well. So
            % I'm using V.mat(:,1) for all views (and take absolute values). 
            % I thought this could be a possible bug in the future, so I made this comment. 9/9/12
            
            % if dummy > max(abs(currentmap.V.mat(:,1)))/2, wh_slice(j) = 0; end % Wani added this line 8/11/12, then see the next line.
            if (dummy - (range(currentmap.SPACE.zcoords)/(length(currentmap.SPACE.zcoords)-1))/2) > .01, wh_slice(j) = 0; end % wani modified this line to fix a weird bug 3/5/13
            numvox = cat(1, squeeze(sum(sum(abs(currentmap.mapdata) > 0))));
                        
        case 'sagittal'
            my_slice_coord = mymontage.slice_mm_coords(j);
            
            [dummy, wh_slice(j)] = min(abs(currentmap.SPACE.xcoords - my_slice_coord));
            
            if (dummy - (range(currentmap.SPACE.xcoords)/(length(currentmap.SPACE.xcoords)-1))/2) > .01, wh_slice(j) = 0; end % Wani modified this line 3/5/13
            
            numvox = cat(1, squeeze(sum(sum(abs(currentmap.mapdata) > 0, 2), 3)));
            
        case 'coronal'
            my_slice_coord = mymontage.slice_mm_coords(j);
            
            [dummy, wh_slice(j)] = min(abs(currentmap.SPACE.ycoords - my_slice_coord));
            
            if (dummy - (range(currentmap.SPACE.ycoords)/(length(currentmap.SPACE.ycoords)-1))/2) > .01, wh_slice(j) = 0; end % Wani modified this line 3/5/13
            
            numvox = cat(1, squeeze(sum(sum(abs(currentmap.mapdata) > 0, 1), 3)));
            
    end
end

% Wani decommented this again because the previous part already fixed the
% problem. 9/9/12
% k = []; % added and modified by Wani 8/11/12 (from here)
% for i = 1:length(wh_slice)-1
%     if wh_slice(i) == wh_slice(i+1)
%         k(end+1) = i+1; 
%     end
% end
% wh_slice(k) = 0; 

k = unique(wh_slice); k(k==0) = [];
voxshown = sum(numvox(k)); %voxshown = sum(numvox(unique(wh_slice))); % (to here)

fprintf('%s montage: %3.0f voxels displayed, %3.0f not displayed on these slices\n', myview, voxshown, sum(numvox) - voxshown);

% SETUP smoothing, contours
% -------------------------------------------------------
if dosmooth
    % for smoothing
    PSF = fspecial('gaussian',10, 4);
    contourmin = .2;
end

if ~docontour
    % surface
    sz = size(SPACE.Z); % dims
    
    switch myview
        case 'axial'
            sz = sz([1 2]); % contour slice size in UPSAMPLED space
        case 'sagittal'
            sz = sz([2 3]);
        case 'coronal'
            sz = sz([1 3]); % contour slice size
    end
    
    cdat = define_cdat(sz, color);
    
    if docolormap % colors are linear mixture of color and mincolor
        
        cdat2 = define_cdat(sz, mincolor);
        
        if dosplitcolor
            % cdat and cdat2 are for positive values
            % cdatminneg and maxneg are for negative values
            cdatminneg = define_cdat(sz, minnegcolor);
            cdatmaxneg = define_cdat(sz, maxnegcolor);
        end
        
    end
    
end % end if docontour

% -----------------------------------------------------------
% Loop through slices to plot
% -----------------------------------------------------------

for j = 1:length(wh_slice) % for j = 1:n - modified by Wani 7/28/12
    if wh_slice(j) ~= 0 % Wani added this if and end 8/11/12
        switch myview
            case 'axial'
                slicedat = currentmap.mapdata(:, :, wh_slice(j));
            case 'sagittal'
                slicedat = squeeze(currentmap.mapdata(wh_slice(j), :, :));
            case 'coronal'
                slicedat = squeeze(currentmap.mapdata(:, wh_slice(j), :));
        end
        
        if any(slicedat(:))
            
            if docontour
                % slow, but do for contours; otherwise render directly
                axes(mymontage.axis_handles(j));
            end
            
            % Z is upsampled slice data; any orientation, not just axial
            switch myview
                case 'axial' % X x Y
                    myx = currentmap.SPACE.Xmm(:, :, wh_slice(j));
                    myy = currentmap.SPACE.Ymm(:, :, wh_slice(j));
                    mynewx = SPACE.Xmm(:, :, 1);
                    mynewy = SPACE.Ymm(:, :, 1);
                    Z = interp2(myx, myy, slicedat, mynewx, mynewy);
                    
                case 'sagittal' % Y x Z; Xmm is for some reason Y mm coords
                    % myx should have all rows the same, myy should have all
                    % cols the same, always always for interp2
                    myy = squeeze(currentmap.SPACE.Xmm(wh_slice(j), :, :));
                    myx = squeeze(currentmap.SPACE.Zmm(wh_slice(j), :, :));
                    mynewy = squeeze(SPACE.Xmm(wh_slice(j), :, :));
                    mynewx = squeeze(SPACE.Zmm(wh_slice(j), :, :));
                    Z = interp2(myx, myy, slicedat, mynewx, mynewy);
                    
                case 'coronal' % X x Z
                    myy = squeeze(currentmap.SPACE.Ymm(:, wh_slice(j), :));
                    myx = squeeze(currentmap.SPACE.Zmm(:, wh_slice(j), :));
                    mynewy = squeeze(SPACE.Ymm(:, wh_slice(j), :));
                    mynewx = squeeze(SPACE.Zmm(:, wh_slice(j), :));
                    Z = interp2(myx, myy, slicedat, mynewx, mynewy); % Wani modified this line. 08/11/12
                   
            end
            
            if dosmooth
                % SMOOTH
                Z = imfilter(Z, PSF, 'replicate', 'conv');
            end
            
            %if j == 8, keyboard, end
            
            if docontour
                % contour outline or plot
                % -----------------------------------------------------------
                
                [c, h] = contourf(mynewy, mynewx, abs(Z), [contourmin, contourmin]);
                
                ch = get(h, 'Children');
                
                whiteh = findobj(ch, 'FaceColor', [1 1 1]);
                %set(whiteh, 'FaceAlpha', 0);
                set(whiteh, 'FaceColor', [.8 .8 .8]);
                
                % white for bg (center-fills), 'flat' for colored
                colorh = findobj(ch, 'FaceColor', 'flat');
                set(colorh, 'FaceColor', color); %, 'FaceAlpha', .8);
                
                if dotrans
                    if ~isempty(transvalue)
                        set(colorh, 'FaceAlpha', transvalue);
                    else
                    set(colorh, 'FaceAlpha', .7); % backward compatible...
                    end
                end
                
                if outline
                    set(colorh, 'FaceAlpha', 0, 'LineWidth', mylinewidth, 'EdgeColor', color);
                end
                
            else
                % surface map method
                % -----------------------------------------------------------
                
                if ~docolormap
                    % single-color map
                    slicecdat = cdat .* repmat(double(abs(Z) > 0), [1 1 3]);
                    
                elseif ~dosplitcolor
                    % color-mapped
                    Zscaled = Z;
                    Zscaled(Zscaled ~= 0 & Zscaled > max(cmaprange)) = max(cmaprange);
                    Zscaled(Zscaled ~= 0 & Zscaled < min(cmaprange)) = min(cmaprange);
                    Zscaled = (Zscaled - min(cmaprange)) ./ (max(cmaprange) - min(cmaprange));
                    
                    w = repmat(Zscaled, [1 1 3]);
                    
                    slicecdat = (w .* cdat) + (1 - w) .* cdat2;
                    
                elseif dosplitcolor
                    % split colormap around zero
                    
                    if max(cmaprange) < 0, [dummy, wh] = max(cmaprange); cmaprange(wh) = abs(min(cmaprange)); end
                    if min(cmaprange) > 0, [dummy, wh] = min(cmaprange); cmaprange(wh) = -(max(cmaprange)); end
                    
                    % make into 4-element: min neg, max neg, min pos, max pos
                    if length(cmaprange) == 2
                        cmaprange = [min(cmaprange) 0 0 max(cmaprange)]; % just like before = all the way to 0
                    end
                    
                    Zscaled = double(Z); % cast as double to avoid weird bug
                    
                    % Zscaled must be transformed to weights from -1 to 1
                    
                    %                 Zscaled(Zscaled > 0 & Zscaled > max(cmaprange)) = max(cmaprange);
                    %                 Zscaled(Zscaled < 0 & Zscaled < min(cmaprange)) = min(cmaprange);
                    
                    %Zscaled = Zscaled ./ max(abs(cmaprange)); % keep scale equal
                    %Zscaled(Zscaled > 0) = Zscaled(Zscaled > 0) ./ max(cmaprange);
                    %Zscaled(Zscaled < 0) = Zscaled(Zscaled < 0) ./ abs(min(cmaprange));
                    
                    % linear scaling into pos and neg range, respectively
                    % allows for 4-element threshold input; e.g., [-6 -3 3 6] to display between 3 and 6
                    
                    Zscaled(Zscaled > 0 & Zscaled < cmaprange(3)) = cmaprange(3) + 100*eps;
                    Zscaled(Zscaled < 0 & Zscaled > cmaprange(2)) = cmaprange(2) - 100*eps;
                    
                    Zscaled(Zscaled > 0) = (Zscaled(Zscaled > 0) - cmaprange(3)) ./ (cmaprange(4) - cmaprange(3));
                    Zscaled(Zscaled < 0) = (Zscaled(Zscaled < 0) - cmaprange(2)) ./ abs(cmaprange(1) - cmaprange(2));
                    
                    
                    
                    % pos only
                    w = repmat(Zscaled, [1 1 3]);
                    
                    slicecdat = (w .* cdat) + (1 - w) .* cdat2;
                    
                    to_keep = double(repmat(Zscaled > 0, [1 1 3]));
                    slicecdat = slicecdat .* to_keep;
                    
                    % now do neg part, then add them
                    %w = repmat(Zscaled, [1 1 3]);
                    w = abs(w); % now values we care about are neg; care about magnitude
                    slicecdat2 = (w .* cdatminneg) + (1 - w) .* cdatmaxneg;
                    
                    to_keep = double(repmat(Zscaled < 0, [1 1 3]));
                    slicecdat2 = slicecdat2 .* to_keep;
                    
                    slicecdat = slicecdat + slicecdat2;
                end
                
                %             switch myview
                %             case 'axial'
                %             h = surf(mymontage.axis_handles(j), SPACE.Ymm(:, :, wh_slice(j)), SPACE.Xmm(:, :, wh_slice(j)), Z, 'FaceColor', 'interp', 'edgecolor', 'none');
                %
                %case 'sagittal'
                
                if ~isa(slicecdat, 'double')
                    keyboard
                end
                
                % h = surf(mymontage.axis_handles(j), mynewy, mynewx, Z, 'FaceColor', 'interp', 'edgecolor', 'none');%, 'FaceAlpha', 'interp');
                % Wani: -ones(size(Z)) is helpful for boundaries for some reasons. - doesn't work with MATLAB 2014b. 

                %h = surf(mymontage.axis_handles(j), mynewy, mynewx, -ones(size(Z)), 'FaceColor', 'interp', 'edgecolor', 'none');
                
                % This works by manually creating a "layer", creating a
                % surface image with the blobs at a Z-value of 1.  The
                % underlay anatomical has a Z-value of 0, so the blobs
                % appear on top. The surface has ones wherever there are
                % blobs to render and zeros elsewhere.
                %
                % Then, 'AlphaDataMapping' is set so that non-zero values are 1
                % (opaque, or transparency-scaled) and zero values are 0
                % (transparent)
                % Finally, 'CData' is set so the color data reflects what is
                % in slicecdat. slicecdat is built manually by taking
                % Z-scores (or whatever the input values are) and scaling
                % them as desired, and creating colors based on a colormap
                % of your choosing (or default one).
                
                h = surf(mymontage.axis_handles(j), mynewy, mynewx, ones(size(Z)), 'FaceColor', 'interp', 'edgecolor', 'none');
                

                % wani: The following options doesn't work with MATLAB 2014b ('FaceAlpha', 'interp')
                %case 'coronal'
                
                %end
                
                if dotrans
                    
                    if ~docolormap, Zscaled = abs(Z); end
                    Zscaled(isnan(Zscaled)) = 0;
                    Zscaled = abs(Zscaled);
                    
                    if ~isempty(transvalue)
                        % constant transparency value
                        set(h, 'AlphaDataMapping', 'scaled', 'AlphaData', transvalue .* double(abs(Z) > 0), 'FaceAlpha', 'interp')
                        
                    else
                        % map transparency with colormap
                        % Note: changed by Tor, 2015. Zscaled values > 1
                        % were messing up transparency scale for underlay.
                        % Max should be 1.
                        set(h, 'AlphaDataMapping', 'scaled', 'AlphaData', Zscaled ./ max(abs(Zscaled(:))), 'FaceAlpha', 'interp');
                    end
                    
                else % Default: No transparency for blobs, transparent outside of blobs
                    set(h, 'AlphaDataMapping', 'scaled', 'AlphaData', double(abs(Z) > 0), 'FaceAlpha', 'interp')
                    
                end
                
                set(h, 'CData', slicecdat)
                
            end % do contour
            
            blobhan{j} = h;
            
        end
    end
end

if ~isempty(blobhan)
    blobhan = cat(1, blobhan{:});
end

% for legend
if dosplitcolor
    mincolor = [mincolor; minnegcolor];
    maxcolor = [color; maxnegcolor];
end


end  % main function




function cdat = define_cdat(sz, color)
% needs a size and a 3-element color vector

cdat = ones([sz 3]);
for ii = 1:3
    cdat(:, :, ii) = color(ii) .* cdat(:, :, ii);
end

end