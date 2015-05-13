function obj = montage(obj, varargin)
% obj = montage(obj, varargin)
%
% Creates montage of slices
% - Solid brain slices or contour outlines
% - Points or text labels or both
% - Flexible slice spacing, colors, marker sizes/styles, axis layout (one row/standard square)
% - axial or saggital orientation
%
% Takes all inputs of plot_points_on_slice.  Optional inputs:
%     {'noslice', 'nodraw'}, drawslice = 0;
%     'color', color = varargin{i+1}; varargin{i+1} = [];
%     'marker', marker = varargin{i+1}; varargin{i + 1} = [];
%     {'wh_slice'}, wh_slice = varargin{i+1}; 
%     {'close', 'closeenough', 'close_enough'}, close_enough =    varargin{i+1};
%     {'sagg','saggital','sagittal'}, orientation = 'sagittal';
%     {'MarkerSize', 'markersize'}, markersize = varargin{i+1};
%     {'MarkerFaceColor', 'markerfacecolor'}, facecolor = varargin{i+1};
%     'solid', disptype = 'solid';
%     'overlay', ovl = varargin{i + 1}; NOTE! DO NOT ENTER THIS HERE; ENTER
%     WHEN YOU INITIALIZE FMRIDISPLAY OBJECT
%     {'text', 'textcodes'}, textcodes = varargin{i + 1}; 
%     {'condf' 'colorcond'}, condf = varargin{i + 1};
%
% In addition:
% 'onerow' : arrange axes in one row
% 'slice_range' : [min max] values in mm for slices to plot
% 'wh_slice' or 'custom_coords' : followed my mm values for slices desired
%      e.g., for cluster/region centers, xyz = cat(1, cl.mm_center)
%      o2 = montage(o2, 'axial', 'wh_slice', xyz, 'onerow');
%      o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
% 'spacing' : followed by inter-slice spacing in mm
%
% Define new axes in existing figure, and use those for montage:
% axh = axes('Position', [0.05 0.4 .1 .5]);
% o2 = montage(o2, 'saggital', 'wh_slice', xyz(1,:), 'existing_axes', axh);
%
% Returns: 
% obj, an fmridisplay object
% obj = 
% 
%   fmridisplay
% 
%   Properties:
%             overlay: [1x105 char]
%               SPACE: [1x1 struct]
%     activation_maps: {[1x1 struct]}
%             montage: {[1x1 struct]}
%             surface: {}
%           orthviews: {}
%             history: {}
%     history_descrip: []
%     additional_info: ''
%
% Examples:
% See help fmridisplay

% Programmers' notes:
% 3/2012 : Fixed bug in slices displayed when choosing exactly 3 custom
% slices.  Transposed coordinates for voxel2mm.  (Tor)


% 
% initialize, if nothing passed in; but you would have to call overloaded
% method, fmridisplay.montage, to invoke this.
if nargin == 0 || isempty(obj) || ~isa(obj, 'fmridisplay')
    obj = fmridisplay;
end

if isempty(obj.SPACE) || ~isstruct(obj.SPACE.V)
    error('fmridisplay is not initialized correctly. run obj = fmridisplay; first and then pass in to this method.')
    
end

donewaxes = 1;  % default - new figure/axes
lightenstr = 'lighten';
myview = 'axial';
disptype = 'solid';  % or contour
doonerow = 0;
spacing = 6; % slice spacing, in mm
ovl = which('SPM8_colin27T1_seg.img');  % which('scalped_avg152T1.img');
textcodes = [];
texthandles = [];
slice_range = 'auto';
custom_coords = 0;
color = 'k';

% ------------------------------------------------------
% parse inputs
% ------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            
            case {'onerow'}, doonerow = 1;
                
            case 'slice_range', slice_range = varargin{i + 1};
                
            case {'spacing'}, spacing = varargin{i+1};
                
            case {'text', 'textcodes'} % do not pass on to slice plot...
                textcodes = varargin{i + 1};
                varargin{i+1} = [];
                varargin{i} = [];
                
            case {'sag', 'sagg','saggital','sagittal'}, myview = 'sagittal';
                
            case {'cor', 'coronal'}, myview = 'coronal';
                    
            case {'condf' 'colorcond'}, condf = varargin{i + 1};
                
            case 'overlay', ovl = varargin{i + 1};
                
            case {'wh_slice', 'custom_coords', 'custom_slices'}
                custom_coords = 1;
                custom_slice_mm = varargin{i + 1};
                
            case 'existing_axes'
                donewaxes = 0; 
                newax = varargin{i + 1};
                
                %otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% SETUP
% -----------------------------------------------
% Volume-level: Set up underlay

fprintf('Load underlay. ')
dat = spm_read_vols(obj.SPACE.V);

fprintf('Define axes. ')
setup_axes();

% fprintf('Lighten edges. ')
% dat = lighten_underlay_edges(dat, 10);

fprintf('Ready. \n')


% Do the work for each slice
% -----------------------------------------------

figure(slices_fig_h)

for i = 1:length(slice_vox_coords)
    
    axes(newax(i));
    
    wh_slice = slice_vox_coords(i);
    
    Z = display_slice(dat, wh_slice, obj.SPACE, myview, lightenstr);
    
    hold on
        
    axis off
end

% register info in fmridisplay object

obj.montage{end + 1} = struct('axis_handles', newax, 'orientation', myview, 'slice_mm_coords', slice_mm_coords);


% ------------------------------------------------------
% INLINE FUNCTIONS
% ------------------------------------------------------

    function setup_axes
        
        myviews = {'axial' 'coronal' 'sagittal'};   % for selecting SPM window
        whview = find(strcmp(myviews, myview));
        if isempty(whview), error('myview must be axial, coronal, or sagittal.'); end
        
        switch whview
            case 1
                if strcmp(slice_range, 'auto')
                    slice_range = [-45 70];
                end
                
                cen = [slice_range(1):spacing:slice_range(2)]';
                if custom_coords
                    cen = custom_slice_mm(:, 3);
                end
                slice_mm_coords = cen;
                cen = [zeros(length(cen), 2) cen];
                xyzvox = mm2voxel(cen', obj.SPACE.V.mat, 1);
                slice_vox_coords = xyzvox(:, 3);
                
            case 2
                if strcmp(slice_range, 'auto')
                    slice_range = [-100 65];
                end
                
                cen = [slice_range(1):spacing:slice_range(2)]';
                if custom_coords
                    cen = custom_slice_mm(:, 2);
                end
                slice_mm_coords = cen;
                cen = [zeros(length(cen),1) cen zeros(length(cen),1)];
                xyzvox = mm2voxel(cen', obj.SPACE.V.mat, 1);
                slice_vox_coords = xyzvox(:, 2);
                
            case 3
                if strcmp(slice_range, 'auto')
                    slice_range = [-70 70];
                end
                
                cen = [slice_range(1):spacing:slice_range(2)]';
                if custom_coords
                    cen = custom_slice_mm(:, 1);
                end
                slice_mm_coords = cen;
                cen = [ cen zeros(length(cen), 2)];
                xyzvox = mm2voxel(cen', obj.SPACE.V.mat, 1);
                slice_vox_coords = xyzvox(:, 1);
        end
        
        
        myviews2 = {'sagittal' 'coronal' 'axial' };  % for selecting coord
        whcoord = strmatch(myview, myviews2) ;
        
        % get text string base
        mystr = {'x = ' 'y = ' 'z = '};
        textbase = mystr{whcoord};
        
        % get optimal number of axes
        num_axes = size(cen, 1);
        rc = ceil(sqrt(num_axes)); % produces empty rows for num_axes=6, catch below. SG
        if mod(num_axes,rc)==0 && num_axes<rc^2  
            rr = ceil(num_axes/rc);
        else
            rr = rc;
        end
        
        if ~donewaxes
            % Break and return here if we are using existing axes.
            slices_fig_h = get(newax(1), 'Parent');
            return
        end
        
        slices_fig_h = figure; %create_figure(myview);
        set(slices_fig_h, 'Color', 'w');
        
        if doonerow
            ss = get(0, 'ScreenSize');
            set(gcf, 'Position', [round(ss(3)/12) round(ss(4)*.9) round(ss(3)*.9) round(ss(4)/7) ])
        end
        
        for i = 1:num_axes
            
            if doonerow
                newax(i) = subplot(1, num_axes, i);
            else
                newax(i) = subplot(rr, rc, i);
            end
            
            axis off;
        end
        
        enlarge_axes(gcf, 1.4);
        
    end % setup axes

end % main function

