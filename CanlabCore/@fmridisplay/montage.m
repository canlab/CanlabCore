function [obj, dat] = montage(obj, varargin)
% Creates montage of slices
%  - Solid brain slices or contour outlines
%  - Points or text labels or both
%  - Flexible slice spacing, colors, marker sizes/styles, axis layout (one row/standard square)
%  - axial or sagittal orientation
%
% :Usage:
% ::
%
%     obj = montage(obj, varargin)
%
% Takes all inputs of plot_points_on_slice.
%
% :Optional Inputs:
%
%  'existing_axes': use existing axes and figure, don't create new ones
%  
%  **'volume_data', dat**: 
%  Keyword 'volume_data' followed by loaded volume data for underlay, passed out in a previous call to fmridisplay.montage 
%
% {'nofigure', 'nofig', 'existing_figure'}
% 
%   **{'noslice', 'nodraw'}:**
%        drawslice = 0;
%
%   **'color':**
%        color = varargin{i+1}; varargin{i+1} = [];
%
%   **'marker':**
%        marker = varargin{i+1}; varargin{i + 1} = [];
%
%   **{'wh_slice'}:**
%        wh_slice = varargin{i+1};
%
%   **{'close', 'closeenough', 'close_enough'}:**
%        close_enough =    varargin{i+1};
%
%   **{'sagg','saggital','sagittal'}:**
%        orientation = 'sagittal';
%
%   **{'MarkerSize', 'markersize'}:**
%        markersize = varargin{i+1};
%
%   **{'MarkerFaceColor', 'markerfacecolor'}:**
%        facecolor = varargin{i+1};
%
%   **'solid':**
%        disptype = 'solid';
%
%   **'overlay':**
%        ovl = varargin{i + 1};
%        NOTE! DO NOT ENTER THIS HERE; ENTER WHEN YOU INITIALIZE FMRIDISPLAY OBJECT
%
%   **{'text', 'textcodes'}:**
%        textcodes = varargin{i + 1};
%
%   **{'condf' 'colorcond'}:**
%        condf = varargin{i + 1};
%
% In addition:
%
%   **'onerow':**
%        arrange axes in one row
%
%   **'slice_range':**
%        [min max] values in mm for slices to plot
%
%   **'wh_slice' or 'custom_coords':**
%        followed my mm values for slices desired
%      e.g., for cluster/region centers, xyz = cat(1, cl.mm_center)
%      o2 = montage(o2, 'axial', 'wh_slice', xyz, 'onerow');
%      o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow');
%
%   **'spacing':**
%        followed by inter-slice spacing in mm
%
%   **'brighten':**
%        followed by a brighten_factor. The maps become brighter 
%        if 0 < brighten_factor <= 1 and darker if -1 <= brighten_factor < 0.
%
% :Outputs:
%
%   **obj:**
%        an fmridisplay object
%
%   **dat:**
%        loaded volume data (3-D) for underlay image, so can be reused of
%        making many montages using the same underlay
%
% :Properties:
%
%  - overlay: [1x105 char]
%  - SPACE: [1x1 struct]
%  - activation_maps: {[1x1 struct]}
%  - montage: {[1x1 struct]}
%  - surface: {}
%  - orthviews: {}
%  - history: {}
%  - history_descrip: []
%  - additional_info: ''
%
% Examples:
% ::
%
%    o2 = fmridisplay; % create starting fmridisplay container object
%
% Define new axes in existing figure, and use those for montage:
% ::
%
%    axh = axes('Position', [0.05 0.4 .1 .5]);
%    o2 = montage(o2, 'saggital', 'wh_slice', xyz(1,:), 'existing_axes', axh);
%
%    o2 = montage(o2, 'saggital', 'slice_range', [-10 10], 'onerow');
%    o2 = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 4);
%    o2 = montage(o2, 'axial', 'slice_range', [-20 30], 'onerow', 'spacing', 8);
%    o2 = montage(o2, 'axial', 'wh_slice', xyz, 'onerow');
%
%    % Parasaggital only:
%    o2 = montage(o2, 'saggital', 'slice_range', [-4 4], 'onerow', 'spacing', 8);
%
% Add/remove blobs and points with fmridisplay.addblobs,
% fmridisplay.addpoints, fmridisplay.removeblobs, fmridisplay.removepoints
%
% :See also:
% fmridisplay, cluster_orthviews, montage_clusters and variants
%
% ..
%    Programmers' notes:
%    3/2012: Fixed bug in slices displayed when choosing exactly 3 custom
%    slices. Transposed coordinates for voxel2mm. (Tor)
% ..

if nargin == 0 || isempty(obj) || ~isa(obj, 'fmridisplay')
    obj = fmridisplay;
end
% initialize, if nothing passed in; but you would have to call overloaded
% method, fmridisplay.montage, to invoke this.


if isempty(obj.SPACE) || ~isstruct(obj.SPACE.V)
    error('fmridisplay is not initialized correctly. run obj = fmridisplay; first and then pass in to this method.')

end

% Multi-panel routing
% -------------------------------------------------------------------------
% The default montage(obj) (and montage(obj, <montagetype>) for any
% canlab_results_fmridisplay layout keyword) composes a MULTI-PANEL display via
% multiview() — the sagittal+axial 'compact' combo by default, matching
% fmri_data.montage. Explicit single-orientation calls ('axial'/'sagittal'/
% 'coronal', slice options, 'existing_axes' — including the per-orientation
% calls multiview/canlab_results_fmridisplay make internally) fall through to the
% single-orientation montage below. The orientation words are deliberately kept
% OUT of the routed set, so montage(obj,'sagittal') stays a single sagittal row.
mv_types = {'compact', 'compact2', 'compact3', 'full', 'full2', 'multirow', ...
    'allslices', 'blobcenters', 'regioncenters', 'full hcp', 'full hcp inflated', ...
    'full no surfaces', 'hcp grayordinates', 'hcp grayordinates subcortex', ...
    'hcp grayordinates compact', 'subcortex compact', 'subcortex full', ...
    'subcortex slices', 'subcortex 3d', 'leftright inout', 'leftright inout subcortex', ...
    'hcp inflated', 'hcp', 'hcp sphere', 'freesurfer inflated', 'freesurfer white', ...
    'freesurfer sphere', 'MNI152NLin2009cAsym white', 'MNI152NLin2009cAsym midthickness', ...
    'MNI152NLin2009cAsym pial', 'MNI152NLin6Asym white', 'MNI152NLin6Asym midthickness', ...
    'MNI152NLin6Asym pial', 'MNI152NLin6Asym sphere'};

% Keywords that request a SINGLE-orientation montage (existing behaviour) and so
% are NOT routed to multiview: the orientation words, the per-slice options, and
% the existing-axes/figure paths (which the multi-panel layouts use internally,
% preventing recursion). Anything else -> multi-panel multiview (default compact).
single_kw = {'axial', 'coronal', 'sag', 'sagg', 'saggital', 'sagittal', ...
    'wh_slice', 'slice_range', 'spacing', 'onerow', 'existing_axes', ...
    'volume_data', 'nofigure', 'nofig', 'existing_figure'};
is_single = false;
for i = 1:numel(varargin)
    if ischar(varargin{i}) && any(strcmpi(varargin{i}, single_kw)), is_single = true; break; end
end

if ~is_single
    mv_type = 'compact'; mv_extra = varargin;              % default multi-panel layout
    for i = 1:numel(varargin)                              % use an explicit montagetype if given
        if ischar(varargin{i}) && any(strcmpi(varargin{i}, mv_types))
            mv_type = varargin{i}; mv_extra(i) = []; break;
        end
    end
    if ~any(strcmpi(varargin, 'noverbose'))
        vname = inputname(1); if isempty(vname), vname = 'o2'; end
        fprintf('Composing multi-panel display. Equivalent call:  %s = multiview(%s, ''%s'');\n', ...
            vname, vname, mv_type);
    end
    obj = multiview(obj, mv_type, mv_extra{:});
    dat = [];   %#ok<NASGU>  % keep the 2nd output defined if requested
    return
end

donewfigure = true; % default - new figure and axes. even if donewaxes is on, this can be turned off to use existing figure
donewaxes = true;  % default - new axes
lightenstr = 'nolighten';
myview = 'axial';
disptype = 'solid';  % or contour
doonerow = 0;
spacing = 6; % slice spacing, in mm

% note: ovl is not functional here...set up in fmridisplay constructor 
%ovl = which('mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');  % symmetric ICBM 152-brain nonlinear iterative registration
%ovl= which('SPM8_colin27T1_seg.img');  % which('scalped_avg152T1.img');
ovl = which('fmriprep20_template.nii.gz');

textcodes = [];
texthandles = [];
slice_range = 'auto';
custom_coords = 0;
color = 'k';
doverbose = true;
brighten_factor = 0;
slice_mm_coords = [];
slices_fig_h = [];
slice_vox_coords = [];
load_volume_data = true;

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
                
            case {'dat', 'volume_data'}
                dat = varargin{i + 1};
                varargin{i+1} = [];
                varargin{i} = [];
                load_volume_data = false;
                
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
                varargin{i} = [];
                
            case {'nofigure', 'nofig', 'existing_figure'}
                donewfigure = false;
                varargin{i} = [];
                
            case 'noverbose'
                doverbose = false;
                
            case {'brighten'}, brighten_factor = varargin{i+1};
                
            % Remove these because they are not relevant here.    
            case 'nosymmetric'
                varargin{i} = [];
                
            case 'regioncenters'
                varargin{i} = [];
                
                %otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% SETUP
% -----------------------------------------------
% Volume-level: Set up underlay

if load_volume_data
    
    if doverbose, fprintf('Load underlay. '), end
    dat = spm_read_vols(obj.SPACE.V);
    
end

if doverbose, fprintf('Define axes. '), end
setup_axes();

% fprintf('Lighten edges. ')
% dat = lighten_underlay_edges(dat, 10);

if doverbose, fprintf('Ready. \n'), end


% Do the work for each slice
% -----------------------------------------------

figure(slices_fig_h) % note: this is slow for repeated calls...

% for i = 1:length(slice_vox_coords)
% 
%     axes(newax(i));
% 
%     wh_slice = slice_vox_coords(i);
% 
%     Z = display_slice(dat, wh_slice, obj.SPACE, myview, lightenstr);
% 
%     hold on
% 
%     axis off
% end



for i = 1:length(slice_vox_coords)
    axes(newax(i));

    wh_slice = slice_vox_coords(i);
    Z = display_slice(dat, wh_slice, obj.SPACE, myview, lightenstr);

    % ✅ Store the actual mm coordinate for this axis
    setappdata(newax(i), 'mm_coord', slice_mm_coords(i));

    hold on
    axis off

    % Hide the per-axes interaction toolbar (the "..." overflow button that
    % appears in the corner of every slice in recent MATLAB releases) and disable
    % the default pan/zoom/rotate interactions, which are meaningless for a slice.
    hide_axes_toolbar(newax(i));
end



% register info in fmridisplay object

obj.montage{end + 1} = struct('axis_handles', newax, 'orientation', myview, 'slice_mm_coords', slice_mm_coords);

% Pull existing blob layers onto this newly added montage, so blobs that were
% added BEFORE this montage (e.g. drawn onto surfaces first) also appear here.
% Mirrors the surface() pull-in. For the usual montage-then-addblobs order there
% are no layers yet, so this is a no-op. See VISUALIZATION_OVERHAUL_NOTES.md.
new_montage_idx = numel(obj.montage);
for k = 1:numel(obj.activation_maps)
    layer = obj.activation_maps{k};
    if ~isfield(layer, 'mapdata') || isempty(layer.mapdata), continue, end
    render_args = {};
    if isfield(layer, 'render_args') && ~isempty(layer.render_args), render_args = layer.render_args; end
    [blobhan, this_cmaprange] = render_blobs(layer, obj.montage{new_montage_idx}, obj.SPACE, render_args{:});
    if ~isempty(blobhan)
        obj.activation_maps{k}.blobhandles = [obj.activation_maps{k}.blobhandles; blobhan];
        if ~isfield(obj.activation_maps{k}, 'cmaprange') || isempty(obj.activation_maps{k}.cmaprange)
            obj.activation_maps{k}.cmaprange = this_cmaprange;
        end
    end
end

if brighten_factor  % if we are brightening - this is a bit slower...
    % set color map to enhance contrast
    datvec = dat(:);
    wh = datvec == 0;
    datvec = datvec(~wh);
    
    cmap = contrast(datvec);
    colormap(brighten(cmap, brighten_factor));
end

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
            % Using existing axes. Ensure we end up with exactly num_axes valid
            % axes: use the provided ones if enough were passed, otherwise EXPAND
            % from the first provided axis into a row of num_axes axes. The expand
            % path is what canlab_results_fmridisplay 'multirow' relies on — it
            % passes a single placeholder axis and lets montage lay out the slices
            % (previously this was gated on there being >1 axis in the figure, so
            % a lone axis with many slices errored: "Index exceeds ... must not
            % exceed 1"). SG 2016/10/26; fixed 2026.
            valid = false(1, numel(newax));
            for k = 1:numel(newax)
                valid(k) = isscalar(newax(k)) && ishandle(newax(k)) && strcmp(get(newax(k), 'type'), 'axes');
            end
            newax = newax(valid);
            if isempty(newax)
                error('fmridisplay:montage:badExistingAxes', ...
                    '''existing_axes'' did not contain any valid axes handles.');
            end

            if numel(newax) >= num_axes
                newax = newax(1:num_axes);
                axis(newax, 'off');
            else
                % Fewer axes than slices: expand from the first provided axis.
                init_pos = get(newax(1), 'Position');
                axis(newax(1), 'off');
                for i = 1:num_axes
                    x_pos = init_pos(1) + (i * (0.85 / num_axes));
                    newax(i) = axes('Position', [x_pos, init_pos(2), init_pos(3), init_pos(4)]);
                    axis(newax(i), 'image');
                end
            end

            slices_fig_h = ancestor(newax(1), 'figure');
            return
        end
        
        if donewfigure
            slices_fig_h = figure; %create_figure(myview);
        elseif canlab_is_uifigure(gcf)
            % 'nofigure' asked to reuse the current figure, but it is a uifigure
            % (the controller): subplot/axes cannot draw into it, so open a fresh
            % traditional figure instead (consistent with the other methods).
            slices_fig_h = figure;
        else
            slices_fig_h = gcf;
        end
        
        set(slices_fig_h, 'Color', 'w');
        
        if doonerow
            ss = get(0, 'ScreenSize');
            set(gcf, 'Position', [round(ss(3)/20) round(ss(4)*.5) round(ss(3)*.9) round(ss(4)/1.5) ])
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


function hide_axes_toolbar(ax)
% Hide the per-axes "..." interaction toolbar on the slice axes. Delegates to the
% shared canlab_hide_axes_toolbar so montages and surfaces behave identically.
canlab_hide_axes_toolbar(ax);
end
