function obj = addblobs(obj, cl, varargin)
% This is a method for fmridisplay objects that adds blobs to one or more montages and other surface plot(s).
%
% :Usage:
% ::
%
%     obj = addblobs(obj, cl, varargin)
%
% See addthreshblobs and multi_threshold methods for a multiple thresholds version
% render_blobs does most of the hard work.
%
% :Inputs:
%
%   **obj:**
%        an fmridisplay object
%
%   **cl:**
%        a region object. If you're using an fmri_data object pass in region(fmri_data_obj)
%
% :Optional inputs:
%
% There are many optional inputs that control display features of blobs.
% These are determined by render_blobs
%
%   **COLOR:**
%
%   **'color':**
%        followed by color vector, e.g., [0 1 1], for solid-color blobs
%
%   **'maxcolor':**
%        followed by color vector for max color range, e.g., [0 1 1]
%
%   **'mincolor':**
%        followed by color vector for min color range, e.g., [0 0 1]
%
%   **'onecolor':**
%        force solid-color blobs
%
%   **'splitcolor':**
%        Positive and negative values are mapped to different
%        colormaps. Default is +=hot, -=cool colors. Followed
%        optionally by cell array with
%        vectors of 4 colors defining max/min for +/- range,
%        {minnegcolor maxnegcolor minposcolor maxposcolor}
%        e.g., {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]}
%
%   **'OUTLINING:**
%
%   **''outline'
%
%   **'linewidth':**
%        followed by width value, e.g., 1
%
%        Note: add 'no_surface' to stop updating the existing surface blobs
%
%
%   **'COLOR RANGE:**
%
%   **'cmaprange':**
%        followed by range of values, e.g., [0 40], [-3 3]. Used in
%        color and transparency setting under some circumstances.
%
%   **'TRANSPARENCY:**
%        {'trans', 'transparent','scaledtransparency', 'constanttrans', [val], 'transvalue', [val]}
%
%   **'trans':**
%        Transparent blobs; with no other input, transparency = 0.75 (1 is opaque, 0 is transparent/invisible)
%
%   **'scaledtransparency':**
%        Transparency is a function of voxel value, lower values are more transparent
%
%   **'transvalue':**
%        Followed by width value, e.g., 1. also 'constanttrans'
%
%
% :Other Options:
%
%   **'smooth':**
%        Smooth blobs
%
%   **'contour':**
%
%   **'interp':**
%      Options to pass to interp2 when rendering blobs on standard
%      brain. See help interp2 for valid interpolation methods
%
%   **'no_surface':**
%        Do not add blobs to surface handles, if they exist
%
%   **'noverbose'** : turn off verbose text reporting
%
%   **'depth'** : followed by number (mm) changes search depth (default=2)
%
% CONTROL OF WHICH MONTAGE
%
%   **'wh_montages':**
%        followed by vector of montage numbers as they appear in
%        the list of registered montages in the fmridisplay object
%
% CONTROL OF WHICH SURFACE
%
%   **'wh_surfaces':**
%        followed by vector of surface numbers as they appear in
%                   the list of registered surfaces in the fmridisplay object
%
% :Examples:
% ::
%
%    obj = addblobs(obj, cl, 'color', [0 1 1]);
%    obj = addblobs(obj, cl, 'color', [0 0 1], 'outline');
%    obj = addblobs(obj, cl, 'color', [0 1 0], 'outline', 'linewidth', 1, 'smooth');
%    obj = addblobs(obj, cl, 'color', [1 0 0], 'smooth', 'cmaprange', [0 40]);
%    obj = addblobs(obj, cl, ... 'wh_montages', 1);
%    obj = addblobs(obj, cl, 'splitcolor', ... 'cmaprange', ... 'trans');
%
% Add only to montage 2 in vector of montages in obj.montage
% ::
%
%    obj = addblobs(obj, cl, 'which_montages', 2);
%
% Map values to the colormap red->yellow.
%
% This uses the default percentile-based mapping so that 20% of voxels
% will have the low color and 20% will have the high color, and the rest will be in between:
% ::
%
%    obj = addblobs(obj, cl, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);
%
% Same, but now Map a specific range of values in image ([0 to .05])
% ::
%
%    obj = addblobs(obj, cl, 'maxcolor', [1 1 0], 'mincolor', [1 0 0], 'cmaprange', [0 .05]);
%
% Separate positive and negative activations and map to a split colormap;
%
% See render_blobs
% ::
%
%    o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]}, 'wh_montages', 1);
%
% It is possible to transparency-map values in a statistic image so you
% can show 'unthresholded' statistic values.  e.g.:
% ::
%
%    o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, 'cmaprange', [-2 2], 'trans', 'scaledtransparency');
%
% ..
%    Copyright Tor Wager, 2011
% ..

% ..
%    Add the volume to activation maps in fmridisplay object
% ..

% Check and convert to region
% -------------------------------------------------------------------
if isstruct(cl) && ~isa(cl, 'region'), cl = cluster2region(cl); end
    
if ~isa(cl, 'region')
    error('cl input must be a region object. Try region() constructor method.');
end

if isempty(cl), return, end

% Map to volume space and interpolate
% -------------------------------------------------------------------

[~, mask] = clusters2mask2011(cl); % turn clusters into mask and volume info

if sum(mask(:)) == 0, warning('No voxels in cl! Empty/zero cl.Z field? Bad cl? No results?'); return, end

XYZ = cat(2, cl(:).XYZ);

dim = max(XYZ, [], 2);

V = struct('mat', cl(1).M, 'dim', dim');

SPACE = map_to_world_space(V);

obj.activation_maps{end + 1} = struct('mapdata', mask, 'V', V, 'SPACE', SPACE, 'blobhandles', [], 'cmaprange', [], ...
    'mincolor', [0 0 1], 'maxcolor', [1 0 0], 'color', []);

% Montage selection and options
% -------------------------------------------------------------------

% select which montages; default = all
wh_montage = 1:length(obj.montage);

whm = strcmp(varargin, 'wh_montages') | strcmp(varargin, 'wh_montage') | strcmp(varargin, 'which_montages') | strcmp(varargin, 'which montages');
if any(whm)
    whm = find(whm);
    wh_montage = varargin{whm(1) + 1};
end

% select which surfaces; default = all
wh_surface = 1:length(obj.surface);
addsurfaceblobs = 1;

whs = strcmp(varargin, 'wh_surfaces') | strcmp(varargin, 'wh_surface') | strcmp(varargin, 'which_surfaces') | strcmp(varargin, 'which surfaces');
if any(whs)
    whs = find(whs);
    wh_surface = varargin{whs(1) + 1};
end

% Default color values
% -------------------------------------------------------------------

% default values
dosplitcolor = 1;
doonecolor = 0;
domaxcolor = 0;
domincolor = 0;

maxposcolor = [1 1 0]; % max pos, most extreme values
minposcolor = [1 .4 .5]; % [.8 .3 0]; % min pos
maxnegcolor = [0 .8 .8]; % [.3 .6 .9]; % max neg
minnegcolor = [0 0 1]; % min neg, most extreme values

% More purply orange and deeper blue:  {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]}
% more straight orange to yellow, cyan: {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}

%add_splitcolor_to_varargin = 0; % internal control, not input option
% pos_colormap = colormap_tor([1 0 .5], [1 1 0], [.9 .6 .1]);  %reddish-purple to orange to yellow
% neg_colormap = colormap_tor([0 0 1], [0 1 1], [.5 0 1]);  % cyan to purple to dark blue

depth=2;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'color', 'onecolor', 'solid'} % many are passed into render_blobs - don't need to do anything
                dosplitcolor = 0;
                 doonecolor = 1;
                
            case 'maxcolor'
                dosplitcolor = 0;
                 domaxcolor = 1;
                 maxcolor = varargin{i + 1};
                
            case 'mincolor'
                dosplitcolor = 0;
                 domincolor = 1;
                 mincolor = varargin{i + 1};
                
            case 'splitcolor'

                if length(varargin) > i && iscell(varargin{i + 1})
                    % we have entered colors
                    
                    splitcolors = varargin{i + 1};
                    if ~iscell(splitcolors) || length(splitcolors) ~= 4
                        error('Enter splitcolor followed by a 4-element cell vector of 4 colors\n{Min_neg} {Max_neg} {Min_pos} {Max_pos}');
                    end
                    maxposcolor = splitcolors{4}; % max pos
                    minposcolor = splitcolors{3}; % min pos
                    maxnegcolor = splitcolors{2}; % max neg
                    minnegcolor = splitcolors{1}; % min neg
                    
                else
                    % use defaults - but add the default colors to varargin
                    % because these are passed on to render_blobs
                    % this will be done below
                    
                end
                
             case 'depth'
                depth = varargin{i + 1};
                
            case 'no_surface'
                addsurfaceblobs = 0;
                
            otherwise %suppress warning because other options passed on.  warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Add relevant colors and args to varargin, because we may have passed in
% keywords without following args, intending to use defaults specified above
% These are passed to render_blobs and also used to update legend registry
% in fmridisplay obj

if dosplitcolor  
    
    mysplitcolors = {minnegcolor maxnegcolor minposcolor maxposcolor};
    
    wh = strcmp(varargin, 'splitcolor');
    varargin(wh) = [];  % remove, add to end
    varargin{end + 1} = 'splitcolor';
    varargin{end + 1} = mysplitcolors;
    
elseif doonecolor
    % do not edit var args, just pass through 
    
else  % add maxcolor/mincolor to varargin
    
    if domaxcolor
        
        wh = strcmp(varargin, 'maxcolor');
        varargin(wh) = [];  % remove, add to end
        varargin{end + 1} = 'maxcolor';
        varargin{end + 1} = maxcolor;
        
    end
    
    if domincolor
        
        wh = strcmp(varargin, 'mincolor');
        varargin(wh) = [];  % remove, add to end
        varargin{end + 1} = 'mincolor';
        varargin{end + 1} = mincolor;
        
    end
    
end
    

% Resampling whole map seems to be too slow: do this in render slice...
% fprintf('Resampling map data to underlay space.');
%
% [resampled_dat, SPACEto] = resample_space(mask, V, obj.SPACE);
%
% fprintf(' Done.\n');
%
% obj.activation_maps{end + 1} = struct('mapdata', resampled_dat, 'V_original', V, 'blobhandles', []);

wh_to_display = length(obj.activation_maps);


% Find valid handles and render blobs on them
% -------------------------------------------------------------------------

currentmap = obj.activation_maps{wh_to_display};

if isempty(currentmap) || ~isfield(currentmap, 'mapdata') || isempty(currentmap.mapdata)
    
    error('Map to display is not a valid activation map on fmridisplay object.');
    
end

% Montages
% -------------------------------------------------------------------------

for i = wh_montage
    
    if length(obj.montage) < i
        error('Requested montage does not exist! Check input montage indices.');
    end
    
    % render, and return color ranges (splitcolor will change these)
    % Note: mincolor, maxcolor used in setting fields in fmridisplay for legend 
    
    [blobhan, cmaprange, mincolor, maxcolor] = render_blobs(currentmap, obj.montage{i}, obj.SPACE, varargin{:});
    
    % Register blob handles and colormap-defining range in object
    if ~isempty(blobhan)
        obj.activation_maps{wh_to_display}.blobhandles = [obj.activation_maps{wh_to_display}.blobhandles; blobhan];
        
        obj.activation_maps{wh_to_display}.cmaprange = cmaprange;
        
        % update color info, for legend (or use defaults..set in
        % render_blobs)_
        for Arg = {'color'}
            
            whmax = find(strcmp(varargin, Arg{1}));
            if ~isempty(whmax)
                obj.activation_maps{wh_to_display}.(Arg{1}) = [obj.activation_maps{wh_to_display}.(Arg{1}) varargin{whmax(1) + 1}];
            end
            
        end
        
        for Arg = {'mincolor' 'maxcolor'}
            
            str = (['obj.activation_maps{wh_to_display}.(Arg{1}) = ' Arg{1} ';']);
            eval(str);
            
        end
        
    end
end


% Surfaces
% -------------------------------------------------------------------------

if addsurfaceblobs
    for i = wh_surface
        
        if length(obj.surface) < i
            warning('Requested surface does not exist! Check input surface indices.');
            continue
        end
        
        % Set color maps for + / - values
        if dosplitcolor || doonecolor || domaxcolor || domincolor
            pos_colormap = colormap_tor(minposcolor, maxposcolor);
            neg_colormap = colormap_tor(minnegcolor, maxnegcolor);
        end
        cluster_surf(cl, depth, 'heatmap', 'colormaps', pos_colormap, neg_colormap, obj.surface{i}.object_handle);
        
    end
end

end  % main function

