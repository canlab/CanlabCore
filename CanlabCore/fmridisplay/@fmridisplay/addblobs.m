function obj = addblobs(obj, cl, varargin)
% obj = addblobs(obj, cl, varargin)
%
% Add blobs to montage and other surface plot(s)
%
% See render_blobs for options
%
% Examples:
%
% obj = addblobs(obj, cl, 'color', [0 1 1]);
% obj = addblobs(obj, cl, 'color', [0 0 1], 'outline');
% obj = addblobs(obj, cl, 'color', [0 1 0], 'outline', 'linewidth', 1, 'smooth');
% obj = addblobs(obj, cl, 'color', [1 0 0], 'smooth', 'cmaprange', [0 40]);
% obj = addblobs(obj, cl, ... 'wh_montages', 1);
% obj = addblobs(obj, cl, 'splitcolor', ... 'cmaprange', ... 'trans');
%
% add only to montage 2 in vector of montages in obj.montage 
% obj = addblobs(obj, cl, 'which_montages', 2); 
% 
% % Map values to the colormap red->yellow.
% % This uses the default percentile-based mapping so that 20% of voxels
% will have the low color and 20% will have the high color, and the rest will be in between:
% obj = addblobs(obj, cl, 'maxcolor', [1 1 0], 'mincolor', [1 0 0]);
%
% % Same, but now Map a specific range of values in image ([0 to .05])
% obj = addblobs(obj, cl, 'maxcolor', [1 1 0], 'mincolor', [1 0 0], 'cmaprange', [0 .05]);
%
% Separate positive and negative activations and map to a split colormap;
% see render_blobs
% o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]}, 'wh_montages', 1);
%
% It is possible to transparency-map values in a statistic image so you
% can show 'unthresholded' statistic values.  e.g.:
% o2 = addblobs(o2, cl, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, 'cmaprange', [-2 2], 'trans');
%
% Copyright Tor Wager, 2011

% Add the volume to activation maps in fmridisplay object
% -------------------------------------------------------------------------

% turn clusters into mask and volume info
[dummy, mask] = clusters2mask2011(cl);

if sum(mask(:)) == 0, warning('No voxels in cl! Empty/zero cl.Z field? Bad cl?'); end

XYZ = cat(2, cl(:).XYZ);

dim = max(XYZ, [], 2);

V = struct('mat', cl(1).M, 'dim', dim');

SPACE = map_to_world_space(V);

obj.activation_maps{end + 1} = struct('mapdata', mask, 'V', V, 'SPACE', SPACE, 'blobhandles', [], 'cmaprange', [], ...
    'mincolor', [0 0 1], 'maxcolor', [1 0 0], 'color', []);

% select which montages; default = all
wh_montage = 1:length(obj.montage);

whm = strcmp(varargin, 'wh_montages') | strcmp(varargin, 'wh_montage') | strcmp(varargin, 'which_montages') | strcmp(varargin, 'which montages');
if any(whm)
    whm = find(whm);
    wh_montage = varargin{whm(1) + 1};
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


% Could add surfaces, etc. here
% -------------------------------------------------------------------------


end  % main function

