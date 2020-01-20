function render_on_surface(obj, surface_handles, varargin)
% Map voxels in obj to isosurface (patch) handles in han and change surface colors according to values in obj
%
% :Usage:
% ::
%
%     render_on_surface(obj, han, [optional inputs])
%
% - Object can be thresholded or unthresholded statistic_image, or other image_vector object
% - Uses only the first image in obj
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2020 Tor Wager
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
% ..
%
% :Inputs:
%
%   **obj:**
%       - image_vector object or subclasses
%       - e.g., fmri_data, thresholded or unthresholded statistic_image obj
%       - Uses only the first image in obj
%
%   **han:**
%       - A vector of one or more handles to patch objects (from surfaces)
%       - e.g., from patch(isosurface()), patch(isocaps()),
%       fmri_data.surface(), fmri_data.isosurface()
%       - Try: surface_handles = addbrain('left_cutaway');
%       - Try: surface_handles = addbrain('coronal_slabs_4');
%
% :Optional Inputs:
%   **'clim':**
%        Followed by [lowerlimit upperlimit] values to define color limits
%
%   **'colormap':**
%        Followed by the name of a Matlab colormap string
%        e.g., 'hot' (default), 'summer', 'winter'
%        - Creates a split colormap with gray values where object values
%        are zero, and colors where object values are nonzero
%
% :Outputs:
%
%   **renders colors on surfaces:**
%   Note that it is not always possible to reset the surfaces to their
%   original colors completely; you may have to redraw them if you change
%   the threshold/voxels that should be colored.
%
% :Examples:
% ::
%
% % Load a standard dataset (emotion regulation test data)
% % Perform a t-test, and render the results on a series of coronal slabs
%
% imgs = load_image_set('emotionreg', 'noverbose');
% t = ttest(imgs, .01, 'unc');
% figure; han = addbrain('coronal_slabs_4');
%
% % When the t map has positive and negative values, creates a special
% % bicolor split colormap that has warm colors for positive vals, and cool colors
% % for negative vals:
% render_on_surface(t, han);
%
% % You can set the color limits (here, in t-values)
% render_on_surface(t, han, 'clim', [-4 4]);
%
% % You can set the colormap to any Matlab colormap:
% render_on_surface(t, han, 'colormap', 'summer');
%
% % If your image is positive-valued only (or negative-valued only), a
% % bicolor split map will not be created:
% 
% t = threshold(t, [2 Inf], 'raw-between');
% render_on_surface(t, han, 'colormap', 'winter', 'clim', [2 6]);


% Programmers' Notes:
% Tor Wager, Jan 2020
% desirable features:
% - Speed
% - control of colormap, can update
% - consistent value->color mapping across surfaces
% - need to know meaning of values, as they often represent statistic values
% - can re-run without changing colors
% - need to plot on isosurface or isocaps objects that use colormap to
% render gray-scale images. 
% This last part makes it a bit more complicated.
% Rendering only on solid-color patch objects from isosurfaces is easier:
%     c = isocolors(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, mesh_struct.voldata, han(i));
%     set(han(i), 'FaceColor', 'interp')
%     c_orig = get(han(i), 'FaceVertexCData');
%     han(i).FaceVertexCData = c;
%
% But to render on either solid-color or gray scale colormapped patches, we want to use
% a split colormap.
%
% Different surface handles (to patches) will have FaceColor is an rgb
% triplet (for solid-color patches) or 'interp' for colormapped patches
% (like isocaps). This will be set to interp for all patches, but we have
% to preserve the existing colormapped grayscale values for interp patches
% by using the split colormap.

if any(~ishandle(surface_handles))
    error('Some surface_handles are not valid handles');
end


% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

nvals = 256;                % number of values for each segment of colormap
colormapname = 'hot';       % default colormap name (all Matlab defined colormaps OK)
custom_colormap = false;
clim = [];                  
allowable_inputs = {'clim' 'colormap' 'colormapname'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case 'colormap'
                
                colormapname = varargin{i+1}; varargin{i+1} = [];
                custom_colormap = true;
                
            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isempty(clim)
    
    % Get range of values in object
    % -----------------------------------------------------------------------
    % These will be used to map color limits
    % we will map values into the range of [1 nvals] where nvals = colormap length
    % Deal with possibility of multiple images in object
    % Exclude zeros
    if isa(obj, 'statistic_image')
        
        sig = logical(obj.sig);
        dat = obj.dat(sig);
        
    else
        
        wh = obj.dat ~= 0 & ~isnan(obj.dat);
        dat = obj.dat(wh);
        
    end
    
    clim = [min(dat(:)) max(dat(:))];
    
end

% Build the colormap
% -----------------------------------------------------------------------
% See notes below

if custom_colormap
    
    cm = split_colormap(nvals, colormapname);

elseif diff(sign(clim))
    
    % Default colormap for objects with - and + values
    cm = hotcool_split_colormap(nvals, clim);
    
else 
    cm = split_colormap(nvals, colormapname);
end

% Map object data values to indices in split colormap
% -----------------------------------------------------------------------
% Define mapping function. We will apply this later to vertex color data
% Defining this here preserves the same mapping across different surfaces
map_function = @(c) nvals + (c - clim(1)) ./ range(clim) .* nvals;


% Reconstruct volume data
% -----------------------------------------------------------------------
% Needed to map to vertices later
[~, ~, mesh_struct] = reconstruct_image(obj);                % get volume data for slices


for i = 1:length(surface_handles)
    % For each surface handle entered
    % -----------------------------------------------------------------------
     
    % Get vertex colors for image (obj) to map
    % -----------------------------------------------------------------------
    % Use isocolors to get colormapped values
    % Convert to index color in split colormap by adding nvals to put in
    % colored range
    % isocolors returns nans sometimes even when no NaNs in data 
    c = isocolors(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, mesh_struct.voldata, surface_handles(i));
    wh = c == 0 | isnan(c);                    % save these to replace with gray-scale later
    
    c_colored = map_function(c);    % Map to colormap indices (nvals = starting range, nvals elements)
    
    % Get the original grayscale values
    % -----------------------------------------------------------------------
    % Exclude colored values
    
    c_gray = get(surface_handles(i), 'FaceVertexCData');
    
    if isempty(c_gray)
        % solid
        c_gray = repmat(round(nvals ./ 2), size(get(surface_handles(i), 'Vertices'), 1), 1);
    else
        c_gray(~wh) = (c_gray(~wh) - min(c_gray(~wh))) ./ range(c_gray(~wh)) .* nvals;
    end
    
    % Now merge graycsale and object-mapped colors 
    % -----------------------------------------------------------------------
    % Everything that was originally zero in isocolors c
    % gets replaced with its original grayscale values in c_gray

    c_colored(wh) = c_gray(wh);
    
    % Set object handle properties 
    % -----------------------------------------------------------------------
    % Change from 'scaled' to 'direct' mode

    set(surface_handles(i), 'FaceVertexCData', c_colored);  % dot indexing sometimes works, sometimes doesn't...depends on handle type
    set(surface_handles(i), 'FaceColor', 'interp')
    set(surface_handles(i), 'CDataMapping', 'direct')
    set(surface_handles(i), 'EdgeColor', 'none');
    
end

colorbar_han = colorbar;
%set(colorbar_han, 'YLim', [nvals+1 2*nvals]);
set(colorbar_han, 'YLim', [nvals+1 2*nvals], 'YTick', [nvals+1 2*nvals], 'YTickLabel', round([clim(1) clim(2)] *100)/100);
set(gca, 'FontSize', 18)

% if isa(obj, 'statistic_image')
%     hotcool_colormap
% else 
%     colormap hot
% 
% end

end

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% Colormap subfunctions
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------

% - This will include nvals (generally 256) gray-scale and nvals colored
% values. If nvals = 256, the first 256 colors will be grayscale, and the
% next 256 color values.
% - If the object has positive and negative values, assume we want a split
% colormap showing negative values in cool colors, positive values in hot
% colors. The split colormap will also have different ranges for negative and
% positive values
% - If object is positive-valued or negative-valued only, assume we want a unipolar colormap 


function cm = hotcool_split_colormap(nvals, clim, varargin)
%
% cm = hotcool_split_colormap(nvals, clim, [lowhot hihot lowcool hicool]), each is [r g b] triplet
%
% Get hotcool split colormap (around 0)
% note: colormaps with these colors match spm_orthviews_hotcool_colormap used in orthviews()
% [.4 .4 .3], [1 1 0] (hot/pos) and [0 0 1], [.4 .3 .4] (cool/neg)

lowhot = [.4 .4 .3];
hihot = [1 1 0];
lowcool = [0 0 1];
hicool = [.4 .3 .4];

if length(varargin) > 0
    if length(varargin) < 4
        error('Enter no optional arguments or 4 [r g b] color triplets');
    end
    lowhot = varargin{1};
    hihot = varargin{2};
    lowcool = varargin{3};
    hicool = varargin{4};
end

cmgray = gray(nvals);

thr = 0;    % Threshold to split colormap into two colors

vec = linspace(clim(1), clim(2), nvals);
pos = find(vec > thr);
np = length(pos);

[hotcm, coolcm] = deal([]);

if np > 0
    hotcm = colormap_tor(lowhot, hihot, 'n', np);
end
hotcm(1, :) = [.5 .5 .5]; % first element is gray

neg = find(vec < -thr);
nn = length(neg);

if nn > 0
    coolcm = colormap_tor(lowcool, hicool, 'n', nn);
end
coolcm(end, :) = [.5 .5 .5]; % last element is gray

cm = [coolcm; hotcm];

cm = [cmgray; cm];
colormap(cm);

end


function cm = split_colormap(nvals, colormapname)
% Get split colormap with gray then colored values
% e.g., cm = split_colormap(256, [2 8], 'hot')
cmgray = gray(nvals);

cmcolor = eval(sprintf('%s(%d)', colormapname, nvals));

cm = [cmgray; cmcolor];
colormap(cm)

end
