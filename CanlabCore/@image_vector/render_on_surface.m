function [cm, colorbar1_han, colorbar2_han] = render_on_surface(obj, surface_handles, varargin)
% Map voxels in obj to isosurface (patch) handles in han and change surface colors according to values in obj
%
% :Usage:
% ::
%
%     cm = render_on_surface(obj, han, [optional inputs])
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
%   **splitcolor':**
%        Positive and negative values are mapped to different
%        colormaps. Default is +=hot, -=cool colors.  Followed
%        optionally by cell array with vectors of 4 colors defining
%        max/min for +/- range, e.g., {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1
%        0]}. This is ignored if pos_colormap or neg_colormap are
%        specified.
%
%   **'indexmap':**
%       Interpret data as indices into colormap. You must specify
%       'colormap' argument as well, otherwise this does nothing and a
%       warning is thrown. If specified, 'nearest' neighbor interpolation
%       is used automatically for volume to surface projection and 'interp'
%       option is ignored.
%
%   **'scaledtransparency':**
%       Transparency is a function of voxel value, lower values are more transparent. This will look very
%       strange without an underlay. To create an underlay render a two duplicate surfaces and invoke
%       render_on_surface only on the second duplicate.
%   
%    **'transcontrast':** 
%       If scaledtransparency is used transparency alpha values are a linear function of the data. This 
%       option makes them a sigmoidal function of the data with higher values disproportionately opaque and
%       lower values disproportionately transparent. This argument is followed by a scalar value which 
%       specifies how steep the sigmoidal part of the contrast curve should be. It works like a contrast 
%       curve in photoshop, but only affects alpha transparency, not colormapping, and therefore does not 
%       alter the data being shown.
%
%    **'interp':**
%       Interpolation method to use when projecting from volume to surface. Options: 'nearest', 'linear'. 
%       If 'indexmap' is specified this option is ignored, and 'nearest' is used automatically.
%
%    **'sourcespace':**
%       If specified together with a targetsurface then nonlinear mapping between the source volume and the 
%       target surface is performed according to the MNIsurf procedure described in Wu, Ngo, Greve et al. (2018)
%       Neuroimage. This is something you absolutely want to do if your sourcespace and targetsurface are
%       supported.
%       Supported sourcespaces = {'MNI152NLin2009cAsym','MNI152NLin6Asym','colin27'}
%
%    **'srcdepth':**
% 	Requires 'sourcespace' specification, and allows further specify desired surface depth to sample volume 
%       at. Typical options might be 'pial', 'white', or 'midthickness'. Requires that a file named 
%       '<sourename>_<srcdepth>_lh.mat' and '<sourename>_<srcdepth>_rh.mat' be in your matlab path. The pial surface 
%       and even midthickness will undersample from deep sulci, but the white surface may conversely suffer from 
%       partial volume effects if you've masked white matter out of your volumetric data. srcdepth can also be a cell 
%       array, in which case multiple depths are sampled and averaged (for linear interpolation) or the mode is taken 
%       (for nearest neighbor interpolation). If there is a modal tie, the default is to use whichever belongs first 
%       in your srcdepth list. Default: {'midthickness','pial','white'}, i.e. midthickness > pial > white for modal 
%       tie breaks.
%
%    **'targetsurface':**
%       If specified together with a targetsurface then nonlinear mapping between the source volume and the 
%       target surface is performed according to the MNIsurf procedure described in Wu, Ngo, Greve et al. (2018)
%       Neuroimage. This is something you absolutely want to do if your sourcespace and targetsurface are
%       supported.
%       Supported targetsurface = {'fsLR_32k', 'fsaverage_164k'}
%
%   
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
%
% Note:
% To erase rendered blobs, use:
% sh = addbrain('eraseblobs', sh);

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
%
% Sept 2020: add solid-color 'color' option
% revise surface() method to use this

if any(~ishandle(surface_handles))
    error('Some surface_handles are not valid handles');
end


% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% cm = [];
nvals = 256;                % number of values for each segment of colormap
colormapname = 'hot';       % default colormap name (all Matlab defined colormaps OK)
custom_colormap = false;
pos_colormap = [];
neg_colormap = [];
clim = [];
axis_handle = get(surface_handles, 'Parent');          % axis handle to apply colormap to; can be altered with varargin
dolegend = true;
doindexmap = false;
splitcolors = {};
doscaledtrans = 0;
enhance_contrast = @(x1)(x1);
sourcespace = [];
srcdepth = {'midthickness','pial','white'}; % default depth to sample a surface from
targetsurface = [];
interp = 'linear';

allowable_sourcespace = {'colin27','MNI152NLin6Asym','MNI152NLin2009cAsym'};
allowable_targetsurface = {'fsaverage_164k','fsLR_32k'};

allowable_keyword_value_pairs = {'clim' 'color' 'colormap' 'colormapname' 'axis_handle' 'pos_colormap' 'neg_colormap', 'sourcespace', 'targetsurface', 'srcdepth','interp'};


% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'color'
                % Special instructions for solid colors
                mycolor = varargin{i+1};
                validateattributes(mycolor, {'numeric'}, {'>=', 0, '<=', 1, 'size', [1 3]})
                
                nvals = 256;
                pos_colormap = repmat(mycolor, nvals, 1);
                
            case 'colormap'
                
                colormapname = varargin{i+1}; varargin{i+1} = [];
                custom_colormap = true;
                
            case 'cmaprange'
                % option to match region.montage method
                
                clim =  varargin{i+1}; varargin{i+1} = [];
                if length(clim)==4
                    clim = clim([1,end]);
                end
                
            case 'splitcolor'                
                splitcolors = varargin{i + 1};
                
                if ~iscell(splitcolors) || length(splitcolors) ~= 4
                    error('Enter splitcolor followed by a 4-element cell vector of 4 colors\n{Min_neg} {Max_neg} {Min_pos} {Max_pos}');
                end
                    
                
            case 'nolegend'
                
                dolegend = false;
                
            case 'indexmap'
                if ~contains('colormap',varargin(cellfun(@ischar,varargin))),...
                        warning('render_on_surface() given ''indexmap'' argument, but no ''colormap'' argument. ''indexmap'' doesn''t do anything without a ''colormap'' argument.');
                end
                
                doindexmap = true;
                interp = 'nearest';
                
            case 'scaledtransparency'
                doscaledtrans = 1;

            case 'transcontrast'
                k = varargin{i+1};
                enhance_contrast = @(x1)(2*((1./(1+exp(-k.*x1)))-0.5));

            case allowable_keyword_value_pairs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
                % eliminate this here because we may have passed many irrelevant values in in surface() call
                %otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if ~iscell(srcdepth)
    srcdepth = {srcdepth};
end

% There are 2 ways to enter custom colormaps.
% 1. By name
% 'colormap', [name of Matlab colormap]
%
% 2. Enter colormaps generated externally
% 'pos_colormap', followed by [n x 3] color vectors and/or
% 'neg_colormap', followed by [n x 3] color vectors
%
% If you enter 'pos_colormap' / 'neg_colormap', these will be stacked and
% added to a gray colormap


% Get range of values in object
% -----------------------------------------------------------------------
% These will be used to map color limits and create legend
% we will map values into the range of [1 nvals] where nvals = colormap length
% Deal with possibility of multiple images in object
% Exclude zeros

[datvec, clim] = get_data_range(obj, clim);

% -------------------------------------------------------------------------
% Define colormap
% -------------------------------------------------------------------------

if ~isempty(pos_colormap) ||  ~isempty(neg_colormap)
    
    if any(strcmp(varargin, 'colormap'))
        error('Entering ''colormap'' argument and name cannot be done when entering custom colormaps in ''pos_colormap''/''neg_colormap''');
    end
    
    if ~isempty(pos_colormap) && ~isempty(neg_colormap) && size(pos_colormap, 1) ~= size(neg_colormap, 1)
        error('Lengths of positive-color and negative-color colormaps must match.')
    end
    
    if ~isempty(pos_colormap)
        validateattributes(pos_colormap, {'numeric'}, {'>=', 0, '<=', 1, 'size', [NaN 3]})
        nvals = size(pos_colormap, 1);
    end
    
    if ~isempty(neg_colormap)
        validateattributes(neg_colormap, {'numeric'}, {'>=', 0, '<=', 1, 'size', [NaN 3]})
        nvals = size(neg_colormap, 1);
    end
    
    if isempty(pos_colormap)
        pos_colormap = hot(nvals);
    end
    
    if isempty(neg_colormap)
        neg_colormap = cool(nvals);
    end
    
    custom_colormap = true;
    
    % Build the colormap
    % -----------------------------------------------------------------------
    % Skip colormap generator, already
    % found the range of colors for pos/neg values to split around 0 here.
    % this sets the colormap for axis_handle
    [cm, kpos, kneg] = hotcool_split_colormap(nvals, clim, axis_handle, pos_colormap(1, :), pos_colormap(end, :), neg_colormap(1, :), neg_colormap(end, :));
    
    % colormapname = [neg_colormap; pos_colormap];   % colormapname is either name or [nvals x 3] matrix
    % nvals = size(colormapname, 1);                 % needs to match for color and gray maps to work right
    
else
    
    % Build the colormap
    % -----------------------------------------------------------------------
    % See notes below
    
    if doindexmap
        
        nvals = min([256, length(unique(obj.dat))]);
        cm = colormapname;
        cm = [0.5,0.5,0.5; cm]; % add gray
        if iscell(axis_handle)
            cellfun(@(x1)colormap(x1,cm), axis_handle)
        else
            colormap(axis_handle,cm);
        end
        
    elseif custom_colormap
        
        [cm, kpos, kneg] = split_colormap(nvals, colormapname, axis_handle); % colormapname is either name or [nvals x 3] matrix
        
    elseif diff(sign(clim))
        
        % Default colormap for objects with - and + values
        [cm, kpos, kneg] = hotcool_split_colormap(nvals, clim, axis_handle, splitcolors{:});
        
    else
        [cm, kpos, kneg] = split_colormap(nvals, colormapname, axis_handle);
    end
    
end % custom posneg or other colormap

% -------------------------------------------------------------------------
% Change colors
% -------------------------------------------------------------------------

% Map object data values to indices in split colormap
% -----------------------------------------------------------------------
% Define mapping function from image values -> colormap indices (1-512).
% We will apply this later to vertex color data (only colored)
% Defining this here preserves the same mapping across different surfaces
% clim(1) mapped to lowest color, clim(2) to highest color
% map_function = @(c) 1 + nvals + (c - clim(1)) ./ range(clim) .* nvals;

map_function = @(c, x1, x2, y1, y2)  y1 + (c - x1) * (y2 - y1) ./ (x2 - x1);


% problem with above is that border gets mapped to low color when interpolating, which is
% blue for split colormap, creating blue outline around orange blobs. SO
% we need to map separately for pos and neg colored blobs
% ***

%% add surface projection information if available        
if ~isempty(sourcespace) && ~isempty(targetsurface)
    if ~ismember(sourcespace, allowable_sourcespace)
        warning('Could not find sourcespace %s in allowed sourcespace list. Will not perform source to target projection.', sourcespace);
        sourcespace = [];
        targetsurface = [];
    elseif ~ismember(targetsurface, allowable_targetsurface)
        warning('Could not find targetsurface %s in allowed sourcespace list. Will not perform source to target projection.', targetsurface);
        sourcespace = [];
        targetsurface = [];
    else
        for i = 1:length(surface_handles)
            % our dependence on surface handle's tag property is a hacky
            % solution here. We really should have a more robust soln
            assert(any(contains(surface_handles(i).Tag,{'right','Right','RIGHT','left','Left','LEFT'})), ...
                'All surface objects must have ''left'' or ''right'' in their Tag property for correct surface projection')
        end
    
        try
            src_sp_lh = cellfun(@(x1)(load(which(sprintf('%s_%s_lh.mat', sourcespace, x1)),'faces','vertices')), srcdepth, 'UniformOutput', false);
            src_sp_rh = cellfun(@(x1)(load(which(sprintf('%s_%s_rh.mat',sourcespace, x1)),'faces','vertices')), srcdepth, 'UniformOutput', false);
        catch e
            fprintf('ERROR: Could not find source surfaces for sourcespace %s', sourcespace);
            rethrow(e);
        end

        assert(length(unique(cellfun(@length, src_sp_lh))) == 1, ...
            'Left hemisphere vertex values differ across templates. Cannot average volumetric values across these srcdepth options');
        assert(length(unique(cellfun(@length, src_sp_rh))) == 1, ...
            'Right hemisphere vertex values differ across templates. Cannot average volumetric values across these srcdepth options');

        % import registration files to transform from source surface (matches volume) to target
        % surface (does not match volume and may not even match the volume's surface model and
        % require interpolation). These were precomputed by code which can be found in
        % CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/src
        if doindexmap
            switch targetsurface
                case 'fsaverage_164k'
                    reg = importdata(which(sprintf('resample_from_%s_to_fsavg_164k_nearestneighbor.mat',sourcespace)));
                case 'fsLR_32k'
                    reg = importdata(which(sprintf('resample_from_%s_to_fsLR_32k_nearestneighbor.mat',sourcespace)));
                otherwise
                    error('Unsupported target surface %s. This error should have been precluded by earlier checks. Strange that you''ve made it here', targetsurface);
            end
        else
            switch targetsurface
                case 'fsaverage_164k'
                    reg = importdata(which(sprintf('resample_from_%s_to_fsavg_164k.mat',sourcespace)));
                case 'fsLR_32k'
                    reg = importdata(which(sprintf('resample_from_%s_to_fsLR_32k.mat',sourcespace)));
                otherwise
                    error('Unsupported target surface %s. This error should have been precluded by earlier checks. Strange that you''ve made it here', targetsurface);
            end
        end
    end
elseif ~isempty(sourcespace) & isempty(targetsurface)
    warning('Source space specified without target surface. Please specify a targetsurface argument to perform accurate surface projection. You have these options,');
    warning(sprintf('targetsurface = %s\n', allowable_targetsurface{:}));
elseif isempty(sourcespace) & ~isempty(targetsurface)    
    warning('target surface specified without source space. Please specify a sourcespace argument to perform accurate surface projection. You have these options:,')
    warning(sprintf('sourcespace = %s\n', allowable_sourcespace{:}));
end


% problem with below is that vertices near empty (zero) voxels get interpolated down to zero
%map_function = @(c) 1 + nvals + (c - 0) ./ range(clim) .* nvals;

% Reconstruct volume data
% -----------------------------------------------------------------------
% Needed to map to vertices later
[~, ~, mesh_struct] = reconstruct_image(obj);                % get volume data for slices

% fixes colormap out of range bug
% BP: sure, but it also asigns colors to empty parts of the cortex. It's
% better to address the colormap out of range bug directly, which is due to
% c_gray being initialized to a 3-vector rather than a 2D matrix. We can
% fix that below instead of using this hacky solution here.
%mesh_struct.voldata(mesh_struct.voldata < min(clim)) = min(clim);

% Deal with edge interpolation effects
% problem is that vertices near empty (zero) voxels get interpolated down to zero
% solution: map all non-zero voxels to lowest limit value
% Need to preserve signs as well.
% Then replace empty voxels with gray
% lower lim to interp to 10% of the way towards zero,
% preserving sign
% lowerclimit = linspace(clim(1), 0, 10);
% lowerclimit = lowerclimit(2);
% mesh_struct.voldata(mesh_struct.voldata == 0) = lowerclimit;
% lowerclimit = 0; ... wh = c == 0 | isnan(c);

% Kludgy fix for interpolation issues around border of blobs
% Changing c or c_colored won't do it because of how gray and colored
% colormaps are stacked. With one color, we interpolate into the gray zone,
% into about the top 256-150 values.
% for i = 256:-1:150
%     cm(i, :) = [.5 .5 .5];
% end
% colormap(cm)


for i = 1:length(surface_handles)
    %% sampling from vol2surf and inter-surface projections
    % The method we use here is what Wu et al. (2018) Neuroimage refer to 
    % as MNIsurf, but we're just applying the mapping here. They were
    % computed ahead of time for faster execution.

    % Get vertex colors for image (obj) to map
    % -----------------------------------------------------------------------
    % Use isocolors to get colormapped values
    % Convert to index color in split colormap by adding nvals to put in
    % colored range
    % isocolors returns nans sometimes even when no NaNs in data

    % interpolate from mesh grid to the surface vertices intersecting
    % the grid
    if isempty(sourcespace) | isempty(targetsurface)
        c = interp3(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, mesh_struct.voldata, ...
            surface_handles(i).Vertices(:,1), surface_handles(i).Vertices(:,2), surface_handles(i).Vertices(:,3), interp);
    else
        % figure out what surface we're dealing with and grab the appropriate
        % source surfaces to sample with
        if contains(get(surface_handles(i),'Tag'),{'left','Left','LEFT'})
            % this is used for mesh interpolation from the volume
            src_sp = src_sp_lh;
    
            % this is used to pull the correct mappings to the target surface 
            % from our reg struct
            vertices = 'vertices_lh';
            weights = 'weights_rh';
        elseif contains(get(surface_handles(i),'Tag'),{'right','Right','RIGHT'})
            % this is used for mesh interpolation from the volume
            src_sp = src_sp_rh;
    
            % this is used to pull the correct mappings to the target surface 
            % from our reg struct
            vertices = 'vertices_rh';
            weights = 'weights_rh';
        else
            error('Could not determine surface laterality. This should have been precluded by earlier checks when loading source space surface, so it''s quite strange you managed to throw this error.')
        end

        % project from source space to target surface by 
        % 1) interpolating to source space specific surface
        % 2) applying precomputed barycentric linear interpolation of source
        % surface to target surface
        if size(surface_handles(i).Vertices,1) ~= size(reg.(vertices),1)
            % if surface vertices dont match the target space
            % vertices default to grid interpolation over MNIsurf
            %
            % this condition will occur with user error but also 
            % when you have multiple surfaces, e.g. a freesurfer 
            % surface + subcortical structures, and then this code 
            % handles the auxillary surfaces.
            warning('Size of surface %d does not match size of target surface. Falling back to naive mesh interpolation.',i);
            c = interp3(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, mesh_struct.voldata, ...
                surface_handles(i).Vertices(:,1), surface_handles(i).Vertices(:,2), surface_handles(i).Vertices(:,3), interp);
        else
            % if target vertices and surface vertices match proceed
            % with MNIsurf vol2surf method
            c_mesh = zeros(size(src_sp{1}.vertices,1), length(src_sp));
            for j = 1:length(src_sp)
                c_mesh(:,j) = interp3(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, mesh_struct.voldata, ...
                    src_sp{j}.vertices(:,1), src_sp{j}.vertices(:,2), src_sp{j}.vertices(:,3), interp);
            end
            if strcmp(interp,'linear')
                c_mesh = mean(c_mesh,2);
            elseif strcmp(interp,'nearest') | doindexmap
                c_mesh = ordered_mode(c_mesh);
            else
                warning('Unknown interpolation option. Defaulting to mean across srcdepths depths.')
                c_mesh = mean(c_mesh,2);
            end
            c_mesh = c_mesh(:); % ensure it's a column vector to satisfy assumption of matmul below

            c = zeros(size(reg.(vertices),1),1);
            for v = 1:size(reg.(vertices),1)
                c(v) = reg.(weights)(v,:)*c_mesh(reg.(vertices)(v,:));
            end
        end
    end
    if doindexmap   
        c(mod(c,1) > 0) = 0;
        c_colored = c;
    else
        %% conversion from raw-values to colormap values
        %if doindexmap, c = round(c); end
        
        % doesn't work to fix interpolation error.
        %     border_percent = prctile(abs(c(c ~= 0)), 80);
        %     whlow = abs(c) < border_percent | isnan(c);
        
        c_colored = c;
    
    
        whpos = c > 0;
        %     kpos = 61;   % which block of 256 colors; depends on colormap
        cpos = map_function(c(whpos), 0, clim(2), (kpos-1)*nvals+1, kpos*nvals); % map into indices in hot cm range of colormap
        c_colored(whpos) = cpos;
        
        whneg = c < 0;
        %     kneg = 55;   % which block of 256 colors
        cneg = map_function(c(whneg), clim(1), 0, (kneg-1)*nvals+1, kneg*nvals); % map into indices in cool cm range of colormap
        c_colored(whneg) = cneg;
    end
            
    wh = c == 0 | isnan(c);                    % save these to replace with gray-scale later
    
    %c_colored = map_function(c);    % Map to colormap indices (nvals = starting range, nvals elements)
    
    % FIX: rescale range so we don't map into gray range at edges due to
    % interpolation
    %     xx = c_colored(~wh); % in-blob vertices
    %     xx = 1 + nvals + max(xx) .* (xx - min(xx)) ./ (range(xx) + nvals);
    %     c_colored(~wh) = xx;
    %cm(256, :) = [1 0 1];
    
    
    %% fix empty values so they show up as gray (some value in the colormap ) and not black (0)
    % Get the original grayscale values
    % -----------------------------------------------------------------------
    % Exclude colored values
    
    c_gray = get(surface_handles(i), 'FaceVertexCData');
    
    if isempty(c_gray) || size(c_gray,1) == 1
        
        % solid
        % c_gray = repmat(round(nvals ./ 2), size(get(surface_handles(i), 'Vertices'), 1), 1);
        if doindexmap
            c_gray = zeros(size(get(surface_handles(i), 'Vertices'), 1), 1);
        else
            c_gray = repmat(0 .* nvals + .5 .* round(nvals), size(get(surface_handles(i), 'Vertices'), 1), 1);
        end
        
    else

        %orig: c_gray(~wh) = (c_gray(~wh) - min(c_gray(~wh))) ./ range(c_gray(~wh)) .* nvals;
        if range(c_gray(wh)) == 0
            % Solid surface color
            % do nothing - skip rescaling
            % -----------------------------------------------------------------------
            %c_colored(wh) = c_gray(wh);
            
        else
            % rescale to gray part of colormap
            c_gray(wh) = (c_gray(wh) - min(c_gray(wh))) ./ range(c_gray(wh)) .* nvals;
            % ***
            
            % Now merge graycsale and object-mapped colors
            % -----------------------------------------------------------------------
            % Everything that was originally zero in isocolors c
            % gets replaced with its original grayscale values in c_gray
            
            %c_colored(wh) = c_gray(wh);
            
        end
    end
    
    c_colored(wh) = c_gray(wh);
    
    %% assign the colors constructed above to the surface object's properties
    % Set object handle properties
    % -----------------------------------------------------------------------
    % Change from 'scaled' to 'direct' mode
    
    % Save FaceVertexCData f or later
    prevdata = get(surface_handles(i), 'FaceVertexCData');
    set(surface_handles(i), 'UserData', prevdata);
    
    if doindexmap
        % increment c_colored by 1 to account for gray
        set(surface_handles(i), 'FaceVertexCData', c_colored+1);
        set(surface_handles(i), 'FaceColor', 'flat');
    else
        set(surface_handles(i), 'FaceVertexCData', c_colored);  % dot indexing sometimes works, sometimes doesn't...depends on handle type
        switch interp
            case 'linear'
                set(surface_handles(i), 'FaceColor', 'interp');
            case 'nearest'
                set(surface_handles(i), 'FaceColor', 'flat');
            otherwise
                warning('Did not understand interpolation option. Using default surface.FaceColor=''inter''');
        end
    end
    
    set(surface_handles(i), 'CDataMapping', 'direct')
    set(surface_handles(i), 'EdgeColor', 'none');
    
    if doscaledtrans
        z = c/max(abs(datvec));
        z = enhance_contrast(z);
        set(surface_handles(i), 'FaceColor', 'interp');
        set(surface_handles(i), 'FaceVertexAlphaData',abs(z));
        set(surface_handles(i), 'FaceAlpha', 'interp');
    end
    
end

% Colorbars - legend
% -----------------------------------------------------------------------

[colorbar1_han, colorbar2_han] = deal([]);
if ~dolegend, return, end

if any(datvec > 0)
    
    bar1axis = axes('Position', [.55 .55 .38 .4]);
    if doindexmap
        colormap(bar1axis, cm(2:end,:));
    else
        colormap(bar1axis, cm(1+(kpos-1)*nvals:kpos*nvals, :));
    end
    colorbar1_han = colorbar(bar1axis);
    set(bar1axis, 'Visible', 'off');
    
    if doindexmap
        set(colorbar1_han, 'YTick', [0 1], 'YTickLabel', [], 'FontSize', 18);
    else
        minpos = min(datvec(datvec > 0));
        ticklabels = [minpos clim(2)];
        ticklabels = arrayfun(@(x1)(x1), ticklabels, 'UniformOutput', false); % make cell array
        for i = 1:length(ticklabels)
            if abs(ticklabels{i}) < 0.01
                % scientific notation
                ticklabels{i} = sprintf('%0.1e',ticklabels{i});
            else
                % standard notation
                ticklabels{i} = sprintf('%0.2f',ticklabels{i});
            end
        end
        set(colorbar1_han, 'YTick', [0 1], 'YTickLabel', ticklabels, 'FontSize', 12);
    end
    
end

if any(datvec < 0)
    
    bar2axis = axes('Position', [.55 .1 .38 .4]);
    if doindexmap
        colormap(bar1axis, cm(2:end,:));
    else
        colormap(bar2axis, cm(1+(kneg-1)*nvals:kneg*nvals, :));
    end
    colorbar2_han = colorbar(bar2axis);
    set(bar2axis, 'Visible', 'off');
    
    if doindexmap
        set(colorbar2_han, 'YTick', [0 1], 'YTickLabel', [], 'FontSize', 18);
    else
        maxneg = max(datvec(datvec < 0));
        ticklabels = [clim(1) maxneg];
        ticklabels = arrayfun(@(x1)(x1), ticklabels, 'UniformOutput', false); % make cell array
        for i = 1:length(ticklabels)
            if abs(ticklabels{i}) < 0.01
                % scientific notation
                ticklabels{i} = sprintf('%0.1e',ticklabels{i});
            else
                % standard notation
                ticklabels{i} = sprintf('%0.2f',ticklabels{i});
            end
        end
        set(colorbar2_han, 'YTick', [0 1], 'YTickLabel', ticklabels, 'FontSize', 12);
    end
    
end

end % main function


% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% Subfunctions
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------


function [datvec, clim] = get_data_range(obj, clim)

if isa(obj, 'statistic_image')
    
    sig = logical(obj.sig);
    dat = obj.dat(sig);
    
else
    
    wh = obj.dat ~= 0 & ~isnan(obj.dat);
    dat = obj.dat(wh);
    
end

datvec = dat(:);

if any(isinf(datvec))
    warning('Some image values are Inf. Expect erratic behavior/errors.');
    whinf = isinf(datvec);
    datvec(whinf) = sign(datvec(whinf)) .* max(abs(datvec(~whinf)));
end

if isempty(clim)
    
    %clim = [min(datvec) max(datvec)];  % clim: data values for min and max, should become min/max colors
    
    if any(datvec < 0) && any(datvec > 0) % split colormap
        
        clim = double([prctile(datvec(datvec < 0), 10) prctile(datvec(datvec > 0), 90)]); % Match defaults for montage in render_blobs
        
    else
        % may need adjustment
        clim = [min(datvec) max(datvec)];
        
    end
    
    % if no variance, we have constant data values - special case.
    % reducing lower clim(1) will effectively map all data to max color value
    if abs(diff(clim)) < 100 * eps
        clim(1) = clim(2) - 1;
    end
    
end

end % function


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


function [cm, kpos, kneg] = hotcool_split_colormap(nvals, clim, axis_handle, varargin)
%
% cm = hotcool_split_colormap(nvals, clim, [lowhot hihot lowcool hicool]), each is [r g b] triplet
%
% Get hotcool split colormap (around 0)
% note: colormaps with these colors match spm_orthviews_hotcool_colormap used in orthviews()
% [.4 .4 .3], [1 1 0] (hot/pos) and [0 0 1], [.4 .3 .4] (cool/neg)

% Split blue/yellow
% lowhot = [.4 .4 .3];
% hihot = [1 1 0];
% lowcool = [0 0 1];
% hicool = [.4 .3 .4];

% Default addblobs - orange/pink
hihot = [1 1 0]; % max pos, most extreme values
lowhot = [1 .4 .5]; % [.8 .3 0]; % min pos
hicool = [0 .8 .8]; % [.3 .6 .9]; % max neg
lowcool = [0 0 1]; % min neg, most extreme values

graybuffer = 20;  % for border - fade to gray in lowest k values

if ~isempty(varargin)
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

% vec = linspace(clim(1), clim(2), nvals);
% pos = find(vec > thr);
% np = length(pos);

[hotcm, coolcm] = deal([]);

if graybuffer
    
    hotcm = [colormap_tor([.5 .5 .5], lowhot, 'n', graybuffer); ...
        colormap_tor(lowhot, hihot, 'n', nvals-graybuffer)];
    
else
    
    hotcm = colormap_tor(lowhot, hihot, 'n', nvals); %np);
    
end

%hotcm(1, :) = [.5 .5 .5]; % first element is gray
% end

% neg = find(vec < -thr);
% nn = length(neg);

% if nn > 0
%coolcm(end, :) = [.5 .5 .5]; % last element is gray
% end

if graybuffer
    
    coolcm = [colormap_tor(lowcool, hicool, 'n', nvals-graybuffer); ...
        colormap_tor(hicool, [.5 .5 .5], 'n', graybuffer)];
    
else
    
    coolcm = colormap_tor(lowcool, hicool, 'n', nvals); % nn);
    
end


% for interpolation issues
% 1:256 is gray-scale anat; 257:512 is buffer; 513:768 is negcm; 769:1024
% is buffer; 1025:1280 is poscm
buffercm = repmat([.5 .5 .5], nvals, 1);

% which block of nvals (usually 256) has colored values
kpos = 61;
kneg = 55;


%cm = [coolcm; buffercm; buffercm; cmgray; hotcm];
cm = [cmgray; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; coolcm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    hotcm];

% Kludgy fix for interpolation issues
% cm(128+3*256:4*256, :) = .5;
% cm(1+3*256:50+3*256, :) = .5;

% cm(128+0*256:1*256, :) = .5;
% cm(1+0*256:50+0*256, :) = .5;

% cm = [cmgray; cm];

% Change colormap(s)
if isempty(axis_handle)
    colormap(cm)
    
elseif iscell(axis_handle)
    
    for i = 1:length(axis_handle)
        colormap(axis_handle{i}, cm);
    end
    
else
    colormap(axis_handle, cm);
end

end


function [cm, kpos, kneg] = split_colormap(nvals, cmcolor, axis_handle)
% Get split colormap with gray then colored values
% e.g., cm = split_colormap(256, [2 8], 'hot')
%
% colormapname is either name or [nvals x 3] matrix

%validateattributes(cmcolor, {'numeric' 'char'})

% which block of nvals (usually 256) has colored values
kpos = 12;
kneg = 12;

cmgray = gray(nvals);

if ischar(cmcolor)
    % Named colormap function
    cmcolor = eval(sprintf('%s(%d)', cmcolor, nvals));
end

validateattributes(cmcolor, {'numeric'}, {'>=', 0, '<=', 1, 'size', [NaN 3]})

buffercm = repmat([.5 .5 .5], nvals, 1);

cm = [cmgray; repmat(buffercm, 10, 1); cmcolor];

% Change colormap(s)
if isempty(axis_handle)
    colormap(cm)
    
elseif iscell(axis_handle)
    
    for i = 1:length(axis_handle)
        colormap(axis_handle{i}, cm);
    end
    
else
    colormap(axis_handle, cm);
end

end

function M = ordered_mode(dat)
    % matlab's mode will return the smallest value when there's ambiguity
    % in the modal value. We don't want that kind of bias in our maps.
    % Instead we prioritize things based on the order of the original
    % values so the user can decide which to prioritize.
    % 
    % by default the mode is computed within row of dat

    [M,F,C] = mode(dat,2);

    for i = 1:length(C)
        if length(C{i})>1
            firstval = [];
            j = 1;
            while isempty(firstval) && j <= size(dat,2)
                if sum(dat(i,:) == dat(i,j)) == F(i)
                    firstval = dat(i,j);
                end
            end
            M(i) = firstval;
        end
    end
end
