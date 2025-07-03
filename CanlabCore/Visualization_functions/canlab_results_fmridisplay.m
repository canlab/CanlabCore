function o2 = canlab_results_fmridisplay(input_activation, varargin)
% :Usage:
% ::
%
%    canlab_results_fmridisplay(input_activation, [optional inputs])
%
% purpose:  This function display fmri results.
%
% :Input:
%
%   **input_activation:**
%        nii, img,
%
%        This image has the blobs you want to
%        display. You can also enter a cl "clusters" structure or
%        "region" object.
%
%        you can also get a thresholded image like the examples used here
%        from a number of places - by thresholding your results in SPM
%        and using "write filtered" to save the image, by creating masks
%        from meta-analysis or anatomical atlases, or by using
%        mediation_brain_results, robust_results_threshold,
%        robust_results_batch_script, threshold_imgs, or object
%        oriented tools including fmri_data and statistic_image objects.
%
% :Optional Inputs:
%
%   **'noblobs':**
%        do not display blobs
%
%   **'outline':**
%        display blob outlines
%
%   **'nooutline':**
%        do not display blob outlines (default)
%
%   **'addmontages':**
%        when entering existing fmridisplay obj, add new montages
%
%   **'noremove':**
%        do not remove current blobs when adding new ones
%
%   **'nofigure':**
%        do not create a new figure (for selected montage sets only)
%
%   **'outlinecolor:**
%        followed by new outline color
%
%   **'splitcolor':**
%        followed by 4-cell new split colormap colors (help fmridisplay or edit code for defaults as example)
%
%   **'montagetype':**
%        Note: for surface plotting MNI surface projection is available. See help render_on_surface for additional 
%        options to specify to enable the necessary transformations. Otherwise naive sampling based on naive surface 
%        vertex coordinates will be used, which in most cases will not correctly map to your data volume.
% 
%        'full'            Axial, coronal, and saggital slices, 4 cortical surfaces
%        'compact'         Midline saggital and two rows of axial slices [the default] 
%        'compact2'        A single row showing midline saggital and axial slices
%        'compact3'        One row of axial slices, midline sagg, and 4 HCP surfaces
%        'multirow'        A series of 'compact2' displays in one figure for comparing different images/maps side by side
%        'regioncenters'   A series of separate axes, each focused on one region
%        'full2'           for a slightly less full montage that avoids colorbar overlap issues
%        'full hcp'        for full montage, but with surfaces and volumes from HCP data
%        'full hcp inflated' for full montage using hcp inflated surfaces
%        'hcp grayordinates' for 4 surfaces and 18 zoomed in subcortical slices
%        'hcp grayordinates compact' for 4 surfaces and 4 zoomed in subcortical slices
%        'hcp grayordinates subcortex'
%                          for zoomed in subcortical slices
%        'hcp inflated'    for a connectome workbench style layout without
%                          volumetric slices
%        'freesurfer inflated' connectome workbench style layout (no volumetric slices) with fsaverage 164k surfaces.
%        'MNI152NLin[2009c|6]Asym [pial|midthickness|white]
%                          Connectome workbench style layout (no volumetric slices) using MNI152NLin2009cAsym or MNI152NLin6Asym 
%                          surfaces sampled at one of three depths. This will work with data that's already in the corresponding 
%                          MNI template space with naive mesh interpolation (no need for special surface projection 
%                          transformations), unlike all other surfaces (see help render_on_surface for details on projection 
%                          options otherwise). When projecting data to other surfaces you need a sampling depth though and these 
%                          MNI space surfaces can be helpful for deciding on a sampling depth to use in your projections (see 
%                          srcdepth argument to render_on_surface). Simply render on the surface corresponding to your data's 
%                          template space and the desired depth to see what data would be extracted from your volumes with that 
%                          srcdepth argument (pial, midthickness, white), and then include that with your srcdepth argument when 
%                          plotting your surface data.
%
%        'compact' [default] for single-figure parasagittal and axials slices.
%
%        'compact2': like 'compact', but fewer axial slices.
%
%        'multirow': followed by number of rows
%           e.g., o2 = canlab_results_fmridisplay([], 'multirow', 2);
%
%        {'blobcenters', 'regioncenters'}: Slices for the center of each blob/region
%        Note: this creates a new figure, tagged
%        'fmridisplay_regioncenters', and is not compatible with 'nofigure'
%        
%   **'noverbose':**
%        suppress verbose output, good for scripts/publish to html, etc.
%
%   **'overlay':**
%        specify anatomical image for montage (not surfaces), followed by
%        image name
%        e.g., o2 = canlab_results_fmridisplay([], 'overlay', 'icbm152_2009_symmetric_for_underlay.img')';
%
%         The default brain for overlays is the brain extracted MNI152NLin2009cAsym 
%         T1 1mm template from templatFlow.
%         For legacy brains based on Keuken et al. 2014 enter as arguments:
%         'overlay', which('keuken_2014_enhanced_for_underlay.img')
%         For legacy SPM8 single subject, enter as arguments:
%         'overlay', which('SPM8_colin27T1_seg.img')
% 
% Other inputs to addblobs and render_on_surface (fmridisplay methods) are allowed, e.g., 'cmaprange', [-2 2], 'trans'
%
% See help fmridisplay
% e.g., 'color', [1 0 0]
%
% You can also input an existing fmridisplay object, and it will use the
% one you have created rather than setting up the canonical slices.
%
% Try "brighten(.4) to make the images brighter.
%
% :Example Script:
% ::
%
%    input_activation = 'Pick_Atlas_PAL_large.nii';
%
%    % set up the anatomical overlay and display blobs
%    % (see the code of this function and help fmridisplay for more examples)
%
%    o2 = canlab_results_fmridisplay(input_activation);
%
%    %% ========== remove those blobs and change the color ==========
%
%    cl = mask2clusters(input_activation);
%    removeblobs(o2);
%    o2 = addblobs(o2, cl, 'color', [0 0 1]);
%
%    %% ========== OR
%
%    r = region(input_activation);
%    o2 = removeblobs(o2);
%    o2 = addblobs(o2, r, 'color', [1 0 0]);
%
%    %% ========== Create empty fmridisplay object on which to add blobs:
%    o2 = canlab_results_fmridisplay
%    o2 = canlab_results_fmridisplay([], 'compact2', 'noverbose');
%
%    %% ========== If you want to start over with a new fmridisplay object,
%    % make sure to clear o2, because it uses lots of memory
%
%    % This image should be on your path in the "canlab_canonical_brains" subfolder:
%
%    input_activation = 'pain-emotion_2s_z_val_FDR_05.img';
%    clear o2
%    close all
%    o2 = canlab_results_fmridisplay(input_activation);
%
%    %% ========== save PNGs of your images to insert into powerpoint, etc.
%    % for your paper/presentation
%
%    scn_export_papersetup(400);
%    saveas(gcf, 'results_images/pain_meta_fmridisplay_example_sagittal.png');
%
%    scn_export_papersetup(350);
%    saveas(gcf, 'results_images/pain_meta_fmridisplay_example_sagittal.png');
%
%    Change colors, removing old blobs and replacing with new ones:
%    o2 = canlab_results_fmridisplay(d, o2, 'cmaprange', [.3 .45], 'splitcolor', {[0 0 1] [.3 0 .8] [.9 0 .5] [1 1 0]}, 'outlinecolor', [.5 0 .5]);
%
%   %% ========== Legend control
%   There is a 'nolegend' option.
%   Colorbar legends are created in render_on_surface
%   You can access and control the handles like this:
%   set(obj.activation_maps{1}.legendhandle, 'Position', [[0.965 0.0994 0.01 0.4037]]);
%
%   %% ========== Colormap range control
%   Range is set automatically by default, and stored in
%   obj.activation_maps{wh_to_display}.cmaprange 
%   You can enter 'cmaprange', followed by inputs in the correct format, to
%   manually control this.
%   
% ..
%    Tor Wager
%    1/27/2012
% ..

% Toggle coordinates visibility
coordinates = 0;
if any(strcmp(varargin, 'coordinates'))
    coordinates = 1;
end

if ~which('fmridisplay.m')
    disp('fmridisplay is not on path.  it is in canlab tools, which must be on your path!')
    return
end

if nargin == 0
    o2 = canlab_results_fmridisplay(region(), 'noblobs', 'nooutline');
    return
end

if ischar(input_activation)
    
    if strcmp(input_activation, 'compact') || strcmp(input_activation, 'compact2') || strcmp(input_activation, 'compact3') || strcmp(input_activation, 'full') ...
            || strcmp(input_activation, 'multirow') || strcmp(input_activation, 'coronal') || strcmp(input_activation, 'sagittal') ...
            || strcmp(input_activation, 'full2') || strcmp(input_activation, 'full hcp')  || strcmp(input_activation, 'full hcp inflated') ...
            || strcmp(input_activation, 'hcp inflated') || strcmp(input_activation, 'freesurfer inflated') ... 
            || strcmp(input_activation, 'full no surfaces') ...
            || strcmp(input_activation, 'freesurfer sphere') || strcmp(input_activation, 'freesurfer white') ...
            || strcmp(input_activation, 'MNI152NLin6Asym white') || strcmp(input_activation, 'MNI152NLin6Asym midthickness') ...
            || strcmp(input_activation, 'MNI152NLin6Asym pial') || strcmp(input_activation, 'MNI152NLin2009cAsym white') ...
            || strcmp(input_activation, 'MNI152NLin2009cAsym midthickness') || strcmp(input_activation, 'MNI152NLin2009cAsym pial') || strcmp(input_activation,'hcp') ...
            || strcmp(input_activation,'hcp grayordinates') || strcmp(input_activation,'hcp grayordinates compact') || strcmp(input_activation,'hcp grayordinates subcortex') ...
            || strcmp(input_activation, 'allslices') || strcmp(input_activation, 'leftright inout') || strcmp(input_activation, 'leftright inout subcortex') ...
            || strcmp(input_activation, 'subcortex full') || strcmp(input_activation, 'subcortex compact') || strcmp(input_activation, 'subcortex 3d') || strcmp(input_activation, 'subcortex slices') 

        
        % Entered no data map; intention is not to plot blobs, just create underlay
        varargin{end + 1} = 'noblobs'; 
        varargin{end + 1} = 'nooutline';
        varargin{end + 1} = input_activation; % so code below finds it later
        % do nothing else for now - this is not an input image
        
    else
        % assume it is an input image
        
        cl = region(fmri_data(input_activation));  % mask2clusters(input_activation);
        
    end
    
elseif isstruct(input_activation) || isa(input_activation, 'region')
    cl = input_activation;
    if ~isa(input_activation, 'region'), cl = cluster2region(cl); end
    
elseif isa(input_activation, 'image_vector')
    cl = region(input_activation);
    
elseif isempty(input_activation)
    % do nothing for now
    
else
    error('I don''t recognize the format of input_activation.  It should be a thresholded mask, clusters, or region object');
end

% process input arguments
% --------------------------------------------
doblobs = true;
dooutline = false;
doaddmontages = false;
doremove = true;
outlinecolor = [0 0 0];
splitcolor = {[0 0 1] [0 .8 .8] [1 .4 .5] [1 1 0]}; % {[0 0 1] [.3 .6 .9] [.8 .3 0] [1 1 0]};  % more straight orange to yellow: {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}
montagetype = 'compact';
orientation = 'axial';
doverbose = true;
%overlay='SPM8_colin27T1_seg.img';
% overlay = 'fmriprep20_template.nii.gz';
overlay = canlab_get_underlay_image;
dofigure = true;

wh = strcmp(varargin, 'overlay');
if any(wh), wh = find(wh); overlay = varargin{wh(1) + 1};  varargin([wh wh+1]) = []; end

wh = strcmp(varargin, 'noblobs');
if any(wh), doblobs = false; varargin(wh) = []; end

wh = strcmp(varargin, 'nooutline');
if any(wh), dooutline = false; varargin(wh) = []; end

wh = strcmp(varargin, 'outline');
if any(wh), dooutline = true; varargin(wh) = []; end

wh = strcmp(varargin, 'addmontages');
if any(wh), doaddmontages = true; varargin(wh) = []; end

wh = strcmp(varargin, 'outlinecolor');
if any(wh), wh = find(wh); outlinecolor = varargin{wh(1) + 1}; end

wh = strcmp(varargin, 'splitcolor');
if any(wh), wh = find(wh); splitcolor = varargin{wh(1) + 1}; end

wh = strcmp(varargin, 'noremove');
if any(wh), doremove = false; varargin(wh) = []; end

wh = strcmp(varargin, 'full');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'multirow');
if any(wh), montagetype = varargin{find(wh)};
    nrows = varargin{find(wh) + 1};
    varargin{find(wh) + 1} = [];
        varargin(wh) = [];
end

wh = strcmp(varargin, 'full hcp');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'full hcp inflated');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'full no surfaces');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'hcp');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'hcp grayordinates');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'hcp grayordinates compact');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'hcp grayordinates subcortex');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'hcp inflated');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'hcp sphere');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'freesurfer sphere');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'freesurfer white');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'freesurfer inflated');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'MNI152NLin2009cAsym white');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'MNI152NLin2009cAsym midthickness');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'MNI152NLin2009cAsym pial');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'MNI152NLin6Asym white');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'MNI152NLin6Asym midthickness');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'MNI152NLin6Asym pial');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'MNI152NLin6Asym sphere');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'full2');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'compact');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'compact2');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'compact3');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'coronal');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'saggital');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'allslices');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'leftright inout');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'leftright inout subcortex');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'subcortex full');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'subcortex compact');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'subcortex 3d');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'subcortex slices');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end


% if these are run before coronal/saggital arguments are parsed they will
% overwrite the blob/regioncenters argument, so run these last.
wh = strcmp(varargin, 'blobcenters');
if any(wh)
    if ismember(montagetype,{'coronal','saggital'})
        orientation = montagetype;
    end
    montagetype = varargin{find(wh)}; 
    varargin(wh) = []; 
end

wh = strcmp(varargin, 'regioncenters');
if any(wh)
    if ismember(montagetype,{'coronal','saggital'})
        orientation = montagetype;
    end
    montagetype = varargin{find(wh)}; 
    varargin(wh) = []; 
end

wh = strcmp(varargin, 'noverbose');
if any(wh), doverbose = false; end

wh = strcmp(varargin, 'nofigure');
if any(wh), dofigure = false; varargin(wh) = []; end

wh = false(1, length(varargin));
for i = 1:length(varargin)
    wh(i) = isa(varargin{i}, 'fmridisplay');
    if wh(i), o2 = varargin{wh}; end
end
varargin(wh) = [];

grayord_xyz = [-24,-12,-6,6,12,24; -36,-17,-12,-5,0,14; -50,-34,-17,-10,6,10]';

%xyz = [-20 -10 -6 -2 0 2 6 10 20]';
xyz = [-40 -20 -10 -2 0 2 10 20 40]';
xyz(:, 2:3) = 0;

if isempty(input_activation)
    % we will skip the blobs, but process other optional input args
    doblobs = false;
    dooutline = false;
end

if ~exist('o2', 'var')
    
    % set up fmridisplay
    % --------------------------------------------
    % you only need to do this once
    % then you can add montages, add and remove blobs, add and remove points (for
    % meta-analysis), etc.
    
    if doverbose
        
        disp('Setting up fmridisplay objects');
        % disp('This takes a lot of memory, and can hang if you have too little.');
        
    end
    
    [opath, ofname, oext] = fileparts(overlay);
    if isempty(opath)
        % check for file in matlab path
        o2 = fmridisplay('overlay', which(overlay));
    else
        % complete path specified
        o2 = fmridisplay('overlay',overlay);
    end
    
    % You can customize these and run them from the command line
    
    switch montagetype
        
        case {'blobcenters', 'regioncenters'}
            % Make a series of montages at center of each region and add blobs to that

            if ~exist('cl', 'var')
                error('Must enter region object to use blobcenters option.');
            end
        
            xyz = cat(1, cl.mm_center);
        
            nr = floor(sqrt(length(cl)));
            nc = ceil(length(cl) ./ nr);
        
            [~, axh] = create_figure('fmridisplay_regioncenters', nr, nc, false, true); 
            set(axh, 'Visible', 'off');
        
            for i = 1:length(cl)
                if i == 1
                    [o2, dat] = montage(o2, orientation, ...
                        'wh_slice', xyz(i,:), ...
                        'onerow', 'existing_axes', axh(i), ...
                        'existing_figure', 'noverbose');
                else
                    o2 = montage(o2, 'volume_data', dat, orientation, ...
                        'wh_slice', xyz(i,:), ...
                        'onerow', 'existing_axes', axh(i), ...
                        'existing_figure', 'noverbose');
                end
        
                if coordinates
                    for i = 1:length(axh)
                        axpos = get(axh(i), 'Position');
                
                        % Bottom center of axis
                        box_x = axpos(1) + axpos(3)/2 - 0.05;
                        box_y = axpos(2) - 0.045;
                
                        % Keep within figure bounds
                        if box_y < 0
                            box_y = 0.005;
                        end
                
                        coordstr = sprintf('x=%.0f y=%.0f z=%.0f', xyz(i,1), xyz(i,2), xyz(i,3));
                        annotation(gcf, 'textbox', [box_x, box_y, 0.1, 0.04], ...
                            'String', coordstr, ...
                            'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', ...
                            'VerticalAlignment', 'top', ...
                            'FontSize', max(6, 10 - round(length(axh)/10)), ...
                            'Interpreter', 'none');
                    end
                end


            end
        
            wh_montages = 1:length(cl);
            
        case 'compact'
            % The default
            
            % saggital
            axh1 = axes('Position', [-0.02 0.4 .17 .17]);
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', [-4 0 0], 'onerow', 'noverbose', 'existing_axes', axh1);
            text(50, -50, 'left');
            drawnow
            
            axh2 = axes('Position', [-0.02 0.6 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'saggital', 'wh_slice', [4 0 0], 'onerow', 'noverbose', 'existing_axes', axh2);
            text(50, -50, 'right');
            drawnow
            
            % sagg center
            axh3 = axes('Position', [.08 0.5 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'saggital', 'wh_slice', [0 0 0], 'onerow', 'noverbose', 'existing_axes', axh3);
            drawnow;
            o2.montage{3}.slice_mm_coords;
         
            % axial bottom
            axh4 = axes('Position', [.10 0.4 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 12, 'noverbose', 'existing_axes', axh4);
            
            % axial top
            axh5 = axes('Position', [.12 0.54 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-46 50], 'onerow', 'spacing', 12, 'noverbose', 'existing_axes', axh5);
            

            % --- Apply coordinate titles if coordinates flag is on ---
            if coordinates
                label_map = containers.Map({'sagittal', 'coronal', 'axial'}, {'x', 'y', 'z'});
            
                for m = 1:numel(o2.montage)
                    mon = o2.montage{m};
                    if ~isfield(mon, 'axis_handles') || isempty(mon.axis_handles)
                        continue;
                    end
            
                    if ~isfield(mon, 'orientation') || ~isKey(label_map, mon.orientation)
                        continue;
                    end
                    label = label_map(mon.orientation);
            
                    for i = 1:numel(mon.axis_handles)
                        ax = mon.axis_handles(i);
                        if isgraphics(ax)
                            mm_coord = getappdata(ax, 'mm_coord');
                            if isempty(mm_coord), continue; end
            
                            % Get axis position in normalized figure units
                            axpos = get(ax, 'Position');  % [x y w h]
                            ann_x = axpos(1) + axpos(3) * 0.5;   % horizontal center
                            ann_y = axpos(2) + axpos(4) + 0.01;  % just above axis
            
                            % Ensure we're in bounds
                            if ann_x >= 0 && ann_x <= 1 && ann_y >= 0 && ann_y <= 1
                                annotation(gcf, 'textbox', [ann_x - 0.05, ann_y, 0.1, 0.03], ...
                                    'String', sprintf('%s=%.0f', label, mm_coord), ...
                                    'HorizontalAlignment', 'center', ...
                                    'VerticalAlignment', 'bottom', ...
                                    'EdgeColor', 'none', ...
                                    'FontSize', 14);
                            end
                        end
                    end
                end
            end

            wh_montages = [1 2 4 5];
            
            % Lines
            axes(o2.montage{3}.axis_handles)
            locs = [o2.montage{4}.slice_mm_coords; o2.montage{5}.slice_mm_coords];
            for i = 1:length(locs)
                %draw_vertical_line(locs(i));
                hh(i) = plot( [-105 65], [locs(i) locs(i)], 'b', 'LineWidth', 1);
            end
            
            % brighten(.5)
            sz = get(0, 'screensize');
            set(gcf, 'Color', 'w', 'Position', [sz(3).*.1 sz(4).*.9 sz(3).*.6 sz(4).*.6]);
         
        case 'multirow' 
            
            % Notes: for some reason, at least in Matlab 2017a, when you
            % use the existing figure the slices scale with the fig
            % position in size. When using montage to create a new figure,
            % they don't seem to do this. could be enlarge_axes, or ???
            
            if dofigure
                slices_fig_h = figure('Color', 'w');
            end
            
            ss = get(0, 'ScreenSize');
            myheightdivisor = 1.5; % 3/nrows;  % controls figure aspect ratio
            set(gcf, 'Position', [round(ss(3)/20) round(ss(4)*.5) round(ss(3)*.9) round(ss(4)/myheightdivisor) ])
            
            %shiftvals = [0:.17:nrows]; % more than we need, but works
            %shiftvals = [0:.24:nrows]; % more than we need, but works
            shiftvals = repmat([0:.24:.75], 1, ceil(nrows/4));  % repeat positions every 4, for adding new figures
            
            for i = 1:nrows
                
                % Can only put 4 on one figure, so create additional
                % figures as needed
                if i > 4 && rem(i, 4) == 1
                    
                    figure('Color', 'w');
                    set(gcf, 'Position', [round(i*ss(3)/20) round(ss(4)*.5) round(ss(3)*.9) round(ss(4)/myheightdivisor) ])
                    
                end
                    
                % saggital
                axh = axes('Position', [-0.02 .75-shiftvals(i) .17 .17]);  % [-0.02 0.15+shiftvals(i) .17 .17]);
                axh(2) = axes('Position', [.022 .854-shiftvals(i) .17 .17]);
                
                if i == 1
                    [o2, dat] = montage(o2, 'saggital', 'slice_range', [-2 2], 'spacing', 4, 'onerow', 'noverbose', 'existing_axes', axh);
                    
                else
                    o2 = montage(o2, 'volume_data', dat, 'saggital', 'slice_range', [-2 2], 'spacing', 4, 'onerow', 'noverbose', 'existing_axes', axh);
                    
                end
                drawnow
                
                % axial
                axh = axes('Position', [.022 0.8-shiftvals(i) .17 .17]);   % [.015 0.15+shiftvals(i) .17 .17]);
                o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-32 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
                drawnow
                
            end

            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
            
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end
            
                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;
            
                    % Determine coordinate label (x/y/z)
                    switch this_montage.orientation
                        case 'sagittal', label_prefix = 'x=';
                        case 'coronal',  label_prefix = 'y=';
                        case 'axial',    label_prefix = 'z=';
                        otherwise,       label_prefix = '';
                    end
            
                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end
            
                        % Axis position in normalized units (relative to figure)
                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
                        ann_y = pos(2) + pos(4) - 0.005;
            
                        % Skip if outside [0, 1] bounds (annotation cannot handle this)
                        if any([ann_x, ann_y, ann_x + 0.05, ann_y + 0.03] > 1) || ...
                           any([ann_x, ann_y] < 0)
                            continue;
                        end
            
                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', ...
                            'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', ...
                            'FontSize', 10);
                    end
                end
            end

        case 'full'
            % --- Saggital ---
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow', 'noverbose');
            shift_axes(-0.02, -0.04);
            sag_coords = xyz(:, 1);  % x-coordinates for saggital

            % --- Coronal ---
            sr_cor = [-40 50];
            cor_coords = (sr_cor(1):8:sr_cor(2))';
            axh = axes('Position', [-0.02 0.37 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'coronal', 'slice_range', sr_cor, 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);

            % --- Axial 1 ---
            sr_axi1 = [-40 50];
            axi_coords1 = (sr_axi1(1):8:sr_axi1(2))';
            axh = axes('Position', [-0.02 0.19 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', sr_axi1, 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);

            % --- Axial 2 ---
            sr_axi2 = [-44 50];
            axi_coords2 = (sr_axi2(1):8:sr_axi2(2))';
            axh = axes('Position', [-0.02 0.01 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', sr_axi2, 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);

            allaxh = findobj(gcf, 'Type', 'axes');
            for i = 1:(length(allaxh)-36)
                pos1 = get(allaxh(i), 'Position');
                pos1(1) = pos1(1) - 0.03;
                set(allaxh(i), 'Position', pos1);
            end

            % --- Coordinate Annotation Overlay ---
            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end
            
                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;
            
                    switch this_montage.orientation
                        case {'saggital', 'sagittal'}, label_prefix = 'x=';
                        case 'coronal',                label_prefix = 'y=';
                        case 'axial',                  label_prefix = 'z=';
                        otherwise,                     label_prefix = '';
                    end
            
                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end
            
                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
            
                        % 💡 Use slightly lower label height for sagittal to stay under surface views
                        if any(strcmp(this_montage.orientation, {'saggital', 'sagittal'}))
                            ann_y = pos(2) + pos(4) - 0.5;
                        else
                            ann_y = pos(2) + pos(4) - 0.01;
                        end
            
                        % Clamp to visible range
                        if ann_x < 0 || ann_x > 1 || ann_y < 0 || ann_y > 1
                            continue;
                        end
            
                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end


            % - Surfaces ---
            o2 = surface(o2, 'axes', [0.1 0.74 .25 .25], 'direction', 'hires left', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.3 0.74 .25 .25], 'direction', 'hires right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.5 0.74 .25 .25], 'direction', 'hires left', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.7 0.74 .25 .25], 'direction', 'hires right', 'orientation', 'lateral');

            wh_montages = [1 2 3 4];
            wh_surfaces = [1 2 3 4];

        case 'full2'

            % saggital
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow', 'noverbose');
            shift_axes(-0.02, -0.04);

            % coronal
            axh = axes('Position', [-0.02 0.37 .1 .17]);
            o2 = montage(o2, 'volume_data', dat, 'coronal', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);

            % axial
            axh = axes('Position', [-0.02 0.19 .1 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            allaxh = findobj(gcf, 'Type', 'axes');
            disp(length(allaxh));
            for i = 1:(length(allaxh)-36)
                pos1 = get(allaxh(i), 'Position');
                pos1(1) = pos1(1) - 0.03;
                set(allaxh(i), 'Position', pos1);
            end

            % --- Coordinate Annotation Overlay ---
            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end

                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;

                    switch this_montage.orientation
                        case {'saggital', 'sagittal'}, label_prefix = 'x=';
                        case 'coronal',  label_prefix = 'y=';
                        case 'axial',    label_prefix = 'z=';
                        otherwise,       label_prefix = '';
                    end

                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end

                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
                        ann_y = pos(2) + pos(4) - 0.005;

                        % Adjust position for sagittal slices to avoid surface overlap
                        if any(strcmp(this_montage.orientation, {'saggital', 'sagittal'}))
                            ann_y = pos(2) + pos(4) - 0.5;
                        end

                        % Ensure annotation is within visible bounds
                        if ann_x < 0 || ann_x > 1 || ann_y < 0 || ann_y > 1
                            continue;
                        end

                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', 'FontSize', 14, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

            % surface
            o2 = surface(o2, 'axes', [0.1 0.74 .25 .25], 'direction', 'hires left', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.3 0.74 .25 .25], 'direction', 'hires right', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.5 0.74 .25 .25], 'direction', 'hires left', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.7 0.74 .25 .25], 'direction', 'hires right', 'orientation', 'lateral');

            wh_montages = [1 2 3 4];
            wh_surfaces = [1 2 3 4];

        case 'compact2'  % creates a new figure

            [o2, dat] = montage(o2, 'axial', 'slice_range', [-32 50], 'onerow', 'spacing', 8, 'noverbose');
            
            % shift all axes down and right
            allaxh = o2.montage{1}.axis_handles;
            for i = 1:length(allaxh)
                pos1 = get(allaxh(i), 'Position');
                pos1(2) = pos1(2) - 0.08;
                pos1(1) = pos1(1) + 0.03;
                
                % enlarge a bit
                pos1(3:4) = pos1(3:4) + .02;
                
                set(allaxh(i), 'Position', pos1);
            end
            
            enlarge_axes(gcf, 1);
            axh = axes('Position', [-0.06 .75-.34 .29 .29]);  % [-0.02 0.15+shiftvals(i) .17 .17]);
            axh(2) = axes('Position', [-0.02 .854-.27 .29 .29]);

            o2 = montage(o2, 'volume_data', dat, 'saggital', 'slice_range', [-2 2], 'spacing', 4, 'onerow', 'noverbose', 'existing_axes', axh);

            % --- Coordinate Annotation Overlay ---
            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};

                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end

                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;

                    switch this_montage.orientation
                        case {'saggital', 'sagittal'}, label_prefix = 'x=';
                        case 'coronal',  label_prefix = 'y=';
                        case 'axial',    label_prefix = 'z=';
                        otherwise,       label_prefix = '';
                    end

                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end

                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
                        ann_y = pos(2) + pos(4) - 0.40;

                        % Adjust position for sagittal slices to avoid surface overlap
                        if any(strcmp(this_montage.orientation, {'saggital', 'sagittal'}))
                            ann_y = pos(2) + pos(4)  - 0.05;
                        end


                        % Skip if outside [0, 1] bounds (annotation cannot handle this)
                        if any([ann_x, ann_y, ann_x + 0.05, ann_y + 0.03] > 1) || ...
                           any([ann_x, ann_y] < 0)
                            continue;
                        end

                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

            wh_montages = [1 2];
            
            brighten(.4)
            

        case 'compact3'  % creates a new figure
            
            % Axial slices
            [o2, dat] = montage(o2, 'axial', 'slice_range', [-32 50], 'onerow', 'spacing', 8, 'noverbose');
            
            % shift all axes down and right
            allaxh = o2.montage{1}.axis_handles;
            for i = 1:length(allaxh)
                pos1 = get(allaxh(i), 'Position');
                pos1(2) = pos1(2) - 0.08;
                pos1(1) = pos1(1) + 0.03;
                
                % enlarge a bit
                pos1(3:4) = pos1(3:4) + .02;
                
                set(allaxh(i), 'Position', pos1);
            end
            
            enlarge_axes(gcf, 1);

            % Medial sagg slices
            axh = axes('Position', [-0.06 .75-.34 .29 .29]);  % [-0.02 0.15+shiftvals(i) .17 .17]);
            axh(2) = axes('Position', [-0.02 .854-.27 .29 .29]);

            o2 = montage(o2, 'volume_data', dat, 'saggital', 'slice_range', [-2 2], 'spacing', 4, 'onerow', 'noverbose', 'existing_axes', axh);

            % --- Coordinate Annotation Overlay ---
            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};

                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end

                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;

                    switch this_montage.orientation
                        case {'saggital', 'sagittal'}, label_prefix = 'x=';
                        case 'coronal',  label_prefix = 'y=';
                        case 'axial',    label_prefix = 'z=';
                        otherwise,       label_prefix = '';
                    end

                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end

                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
                        ann_y = pos(2) + pos(4) - 0.50;

                        % Adjust position for sagittal slices to avoid surface overlap
                        if any(strcmp(this_montage.orientation, {'saggital', 'sagittal'}))
                            ann_y = pos(2) + pos(4)  - 0.05;
                        end


                        % Skip if outside [0, 1] bounds (annotation cannot handle this)
                        if any([ann_x, ann_y, ann_x + 0.05, ann_y + 0.03] > 1) || ...
                           any([ann_x, ann_y] < 0)
                            continue;
                        end

                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

            wh_montages = [1 2];
            
            brighten(.4)

            % surfaces
            o2 = surface(o2, 'axes', [0.1+.1 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.27+.1 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.44+.1 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.61+.1 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'lateral');
            
            wh_surfaces = [1:8];

        case 'coronal'
            o2 = montage(o2, 'coronal', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose');
            wh_montages = 1;

            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end

                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;

                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end

                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
                        ann_y = pos(2) + pos(4) - 0.51;

                        if any([ann_x, ann_y, ann_x + 0.05, ann_y + 0.03] > 1) || any([ann_x, ann_y] < 0)
                            continue;
                        end

                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('y=%.0f', coords(a)), 'FitBoxToText', 'on', ...
                            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

        case {'saggital', 'sagittal'}
            o2 = montage(o2, montagetype, 'wh_slice', xyz, 'onerow', 'noverbose');
            wh_montages = 1;

            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end

                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;

                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end

                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
                        ann_y = pos(2) + pos(4) - 0.55;

                        if any([ann_x, ann_y, ann_x + 0.05, ann_y + 0.03] > 1) || any([ann_x, ann_y] < 0)
                            continue;
                        end

                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('x=%.0f', coords(a)), 'FitBoxToText', 'on', ...
                            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

        case 'allslices'

            % --- Sagittal montage ---
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow', 'noverbose');
            shift_axes(-0.02, -0.04);
            sag_coords = xyz(:, 1);  % x-coordinates
        
            % --- Coronal montage ---
            sr = [-40 50];
            cor_coords = (sr(1):8:sr(2))';
            axh = axes('Position', [-0.02 0.37 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'coronal', ...
                'slice_range', sr, 'onerow', 'spacing', 8, ...
                'noverbose', 'existing_axes', axh);
        
            % --- Axial montage ---
            sr = [-40 50];
            axi_coords = (sr(1):8:sr(2))';
            axh = axes('Position', [-0.02 0.19 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', ...
                'slice_range', sr, 'onerow', 'spacing', 8, ...
                'noverbose', 'existing_axes', axh);
        
            % --- Apply coordinate titles if flag is on ---
            if coordinates && isprop(o2, 'montage')
                coord_labels = {'x', 'y', 'z'};
                coord_sets = {sag_coords, cor_coords, axi_coords};
        
                for m = 1:min(3, numel(o2.montage))
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || isempty(this_montage.axis_handles)
                        continue;
                    end
        
                    axlist = this_montage.axis_handles;
                    coords = coord_sets{m};
                    label = coord_labels{m};
        
                    for i = 1:min(numel(axlist), numel(coords))
                        title(axlist(i), sprintf('%s=%.0f', label, coords(i)), 'FontSize', 10);
                    end
                end
            end
        
            wh_montages = [1 2 3];
            
        case 'full hcp'
            % saggital
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow', 'noverbose');
            shift_axes(-0.02, -0.04);surface right
            
            % coronal
            axh = axes('Position', [-0.02 0.37 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'coronal', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            % axial
            axh = axes('Position', [-0.02 0.19 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            axh = axes('Position', [-0.02 0.01 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-44 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            allaxh = findobj(gcf, 'Type', 'axes');
            disp(length(allaxh));
            for i = 1:(length(allaxh)-36)
                pos1 = get(allaxh(i), 'Position');
                pos1(1) = pos1(1) - 0.03;
                set(allaxh(i), 'Position', pos1);
            end


            % --- Coordinate Annotation Overlay ---
            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end
            
                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;
            
                    switch this_montage.orientation
                        case {'saggital', 'sagittal'}, label_prefix = 'x=';
                        case 'coronal',                label_prefix = 'y=';
                        case 'axial',                  label_prefix = 'z=';
                        otherwise,                     label_prefix = '';
                    end
            
                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end
            
                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
            
                        % 💡 Use slightly lower label height for sagittal to stay under surface views
                        if any(strcmp(this_montage.orientation, {'saggital', 'sagittal'}))
                            ann_y = pos(2) + pos(4) - 0.5;
                        else
                            ann_y = pos(2) + pos(4) - 0.01;
                        end
            
                        % Clamp to visible range
                        if ann_x < 0 || ann_x > 1 || ann_y < 0 || ann_y > 1
                            continue;
                        end
            
                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

            % surface
            o2 = surface(o2, 'axes', [0.1 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.3 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.5 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.7 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'lateral');
            
            wh_montages = [1 2 3 4];
            wh_surfaces = [1:8];


        case 'full hcp inflated'
            % saggital
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow', 'noverbose');
            shift_axes(-0.02, -0.04);
            
            % coronal
            axh = axes('Position', [-0.02 0.37 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'coronal', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            % axial
            axh = axes('Position', [-0.02 0.19 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            axh = axes('Position', [-0.02 0.01 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-44 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            allaxh = findobj(gcf, 'Type', 'axes');
            disp(length(allaxh));
            for i = 1:(length(allaxh)-36)
                pos1 = get(allaxh(i), 'Position');
                pos1(1) = pos1(1) - 0.03;
                set(allaxh(i), 'Position', pos1);
            end

            % --- Coordinate Annotation Overlay ---
            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end
            
                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;
            
                    switch this_montage.orientation
                        case {'saggital', 'sagittal'}, label_prefix = 'x=';
                        case 'coronal',                label_prefix = 'y=';
                        case 'axial',                  label_prefix = 'z=';
                        otherwise,                     label_prefix = '';
                    end
            
                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end
            
                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
            
                        % 💡 Use slightly lower label height for sagittal to stay under surface views
                        if any(strcmp(this_montage.orientation, {'saggital', 'sagittal'}))
                            ann_y = pos(2) + pos(4) - 0.5;
                        else
                            ann_y = pos(2) + pos(4) - 0.01;
                        end
            
                        % Clamp to visible range
                        if ann_x < 0 || ann_x > 1 || ann_y < 0 || ann_y > 1
                            continue;
                        end
            
                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

            % surface
            o2 = surface(o2, 'axes', [0.1 0.74 .25 .25], 'direction', 'hcp inflated right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.3 0.74 .25 .25], 'direction', 'hcp inflated left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.5 0.74 .25 .25], 'direction', 'hcp inflated right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.7 0.74 .25 .25], 'direction', 'hcp inflated left', 'orientation', 'lateral');
            
            wh_montages = [1 2 3 4];
            wh_surfaces = [1:8];


        case 'full no surfaces'
            % saggital
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow', 'noverbose');
            % shift_axes(-0.02, -0.04);
            
            % coronal
            axh = axes('Position', [-0.02 0.37 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'coronal', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            % axial
            axh = axes('Position', [-0.02 0.19 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            axh = axes('Position', [-0.02 0.01 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-44 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            allaxh = findobj(gcf, 'Type', 'axes');
            disp(length(allaxh));
            for i = 1:(length(allaxh)-36)
                pos1 = get(allaxh(i), 'Position');
                pos1(1) = pos1(1) - 0.03;
                set(allaxh(i), 'Position', pos1);
            end

            % --- Coordinate Annotation Overlay ---
            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end
            
                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;
            
                    switch this_montage.orientation
                        case {'saggital', 'sagittal'}, label_prefix = 'x=';
                        case 'coronal',                label_prefix = 'y=';
                        case 'axial',                  label_prefix = 'z=';
                        otherwise,                     label_prefix = '';
                    end
            
                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end
            
                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
            
                        % 💡 Use slightly lower label height for sagittal to stay under surface views
                        if any(strcmp(this_montage.orientation, {'saggital', 'sagittal'}))
                            ann_y = pos(2) + pos(4) - 0.5;
                        else
                            ann_y = pos(2) + pos(4) - 0.01;
                        end
            
                        % Clamp to visible range
                        if ann_x < 0 || ann_x > 1 || ann_y < 0 || ann_y > 1
                            continue;
                        end
            
                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

            wh_montages = [1 2 3 4];

        case 'hcp grayordinates'
            % saggital
            f1 = gcf;
            mainLayout = tiledlayout(f1,1,5);
            surfLayout = tiledlayout(mainLayout, 2, 2, 'Parent', mainLayout, 'TileSpacing', 'compact', 'Padding', 'none');
            surfLayout.Layout.Tile = 1; % Position the leftLayout in the first two tiles of the mainLayout
            surfLayout.Layout.TileSpan = [1, 2]; % Span two tiles
            
            % surface
            ax = {};
            for j = 1:4, ax{j} = nexttile(surfLayout); end
            o2 = surface(o2, 'axes', ax{1}, 'direction', 'hcp inflated right', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', ax{2}, 'direction', 'hcp inflated left', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', ax{3}, 'direction', 'hcp inflated left', 'orientation', 'medial');
            o2 = surface(o2, 'axes', ax{4}, 'direction', 'hcp inflated right', 'orientation', 'lateral');
            
            n_col = size(grayord_xyz,1);
            n_row = size(grayord_xyz,2);
            volLayout = tiledlayout(mainLayout, n_row, n_col, 'Parent', mainLayout, 'TileSpacing', 'tight', 'Padding', 'none');
            volLayout.Layout.Tile = 3; % Position the leftLayout in the third tile of the mainLayout
            volLayout.Layout.TileSpan = [1, 3]; % Span three tiles
            ax_vol=[];
            for j = 1:n_row, for k = 1:n_col, ax_vol(j,k) = nexttile(volLayout); end; end
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes',ax_vol(1,:));
            for j=1:n_col
                set(ax_vol(1,j),'XLim',[-100,30],'YLim',[-70,25]);
                title(ax_vol(1,j), sprintf('x=%.0f', grayord_xyz(j,1)), 'FontSize', 10);
            end

            o2 = montage(o2, 'volume_data', dat, 'coronal', 'wh_slice', grayord_xyz, 'onerow','noverbose', 'existing_axes', ax_vol(2,:));
            for j=1:n_col
                set(ax_vol(2,j),'XLim',[-40,40],'YLim',[-70,25]);
                title(ax_vol(2,j), sprintf('y=%.0f', grayord_xyz(j,2)), 'FontSize', 10);
            end

            o2 = montage(o2, 'volume_data', dat, 'axial', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes', ax_vol(3,:));
            for j=1:n_col
                set(ax_vol(3,j),'XLim',[-60,60],'YLim',[-100,30]);
                title(ax_vol(3,j), sprintf('z=%.0f', grayord_xyz(j,3)), 'FontSize', 10);
            end

            allaxh = findobj(gcf, 'Type', 'axes');

            wh_montages = [1 2 3];
            wh_surfaces = [1:4];

        case 'hcp grayordinates subcortex'
            % saggital
            f1 = gcf;
            n_col = size(grayord_xyz,1);
            n_row = size(grayord_xyz,2);
            volLayout = tiledlayout(f1, n_row, n_col, 'TileSpacing', 'tight', 'Padding', 'none');
            ax_vol=[];
            for j = 1:n_row, for k = 1:n_col, ax_vol(j,k) = nexttile(volLayout); end; end
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes',ax_vol(1,:));
            for j=1:n_col, set(ax_vol(1,j),'XLim',[-100,30],'YLim',[-70,25]); end

            % coronal
            o2 = montage(o2, 'volume_data', dat, 'coronal', 'wh_slice', grayord_xyz, 'onerow','noverbose', 'existing_axes', ax_vol(2,:));
            for j=1:n_col, set(ax_vol(2,j),'XLim',[-40,40],'YLim',[-70,25]); end

            % axial
            o2 = montage(o2, 'volume_data', dat, 'axial', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes', ax_vol(3,:));
            for j=1:n_col, set(ax_vol(3,j),'XLim',[-60,60],'YLim',[-100,30]); end

            allaxh = findobj(gcf, 'Type', 'axes');

            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end

                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;

                    switch this_montage.orientation
                        case {'saggital', 'sagittal'}, label_prefix = 'x=';
                        case 'coronal',  label_prefix = 'y=';
                        case 'axial',    label_prefix = 'z=';
                        otherwise,       label_prefix = '';
                    end

                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end

                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
                        ann_y = pos(2) + pos(4) + .12;

                        if any([ann_x, ann_y, ann_x + 0.05, ann_y + 0.03] > 1) || any([ann_x, ann_y] < 0)
                            continue;
                        end

                        fig_handle = ancestor(ax, 'figure');
                        if isempty(fig_handle) || ~ishandle(fig_handle)
                            continue;
                        end

                        annotation(fig_handle, 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

            wh_montages = [1 2 3];

        case 'hcp grayordinates compact'
            overlay = canlab_get_underlay_image;
            grayord_xyz = [-24,-12,-6,6,12,24; -36,-17,-12,-5,0,14; -50,-34,-17,-10,6,10]';
            % grayord_xyz = [-24,-12,-6,6,12,24; -36,-17,-12,-5,0,14; -50,-24,-17,-10,6,10]';
        
            o2 = fmridisplay('overlay', which(overlay));
        
            t0 = tiledlayout(5,5,'Padding','none','TileSpacing','tight');
            t1 = nexttile();
            t1.Layout.TileSpan = [2,2];
            axis off;
            t2 = nexttile();
            t2.Layout.TileSpan = [2,2];
            axis off;
            t3 = nexttile();
            t3.Layout.TileSpan = [2,2];
            t3.Layout.Tile = 11;
            axis off;
            t4 = nexttile();
            t4.Layout.TileSpan = [2,2];
            t4.Layout.Tile = 13;
            axis off;
            o2 = surface(o2, 'axes', t1, 'direction', 'hcp inflated right', 'orientation', 'medial','disableVis3d');     
            o2 = surface(o2, 'axes', t2, 'direction', 'hcp inflated left', 'orientation', 'lateral','disableVis3d');
            o2 = surface(o2, 'axes', t3, 'direction', 'hcp inflated left', 'orientation', 'medial','disableVis3d');     
            o2 = surface(o2, 'axes', t4, 'direction', 'hcp inflated right', 'orientation', 'lateral','disableVis3d');
            
            t5 = nexttile();
            t5.Layout.Tile = 21;
            axis off
            t6 = nexttile();
            t6.Layout.Tile = 22;
            axis off;
            t7 = nexttile();
            t7.Layout.Tile = 23;
            axis off;
            t8 = nexttile();
            t8.Layout.Tile = 24;
            axis off
        
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz(1,:), 'onerow', 'noverbose', 'existing_axes',t5);
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz(3,:), 'onerow', 'noverbose', 'existing_axes',t6);
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz(5,:), 'onerow', 'noverbose', 'existing_axes',t7);
            [o2, dat] = montage(o2, 'axial', 'wh_slice', grayord_xyz(2,:), 'onerow', 'noverbose', 'existing_axes',t8);
        
            title(t5,'X = -24','FontSize',10)
            title(t6,'X = -6','FontSize',10)
            title(t7,'X = 12','FontSize',10)
            title(t8, 'Z = -34','FontSize',10)
            % title(t8, 'Z = -24','FontSize',10)
            for t = [t5,t6,t7]
                set(t,'XLim',[-100,30],'YLim',[-70,25]);
            end
            set(t8,'XLim',[-60,60],'YLim',[-100,30]);
        
            wh_surfaces = [1:4];

        case 'subcortex compact'
            overlay = canlab_get_underlay_image;
            grayord_xyz = [-24,-12,-6,6,12,24; -36,-17,-12,-5,0,14; -50,-34,-17,-10,6,10]';
        
            o2 = fmridisplay('overlay', which(overlay));
        
            t0 = tiledlayout(5,5,'Padding','none','TileSpacing','tight');
            t1 = nexttile();
            t1.Layout.TileSpan = [2,2];
            axis off;
            t2 = nexttile();
            t2.Layout.TileSpan = [2,2];
            axis off;
            t3 = nexttile();
            t3.Layout.TileSpan = [2,2];
            t3.Layout.Tile = 11;
            axis off;
            t4 = nexttile();
            t4.Layout.TileSpan = [2,2];
            t4.Layout.Tile = 13;
            axis off;
            o2 = surface(o2, 'axes', t1, 'direction', 'caudate left');     
            o2 = surface(o2, 'axes', t2, 'direction', 'caudate right');
            o2 = surface(o2, 'axes', t3, 'direction', 'brainstem left');     
            o2 = surface(o2, 'axes', t4, 'direction', 'brainstem right');
            
            t5 = nexttile();
            t5.Layout.Tile = 21;
            axis off
            t6 = nexttile();
            t6.Layout.Tile = 22;
            axis off;
            t7 = nexttile();
            t7.Layout.Tile = 23;
            axis off;
            t8 = nexttile();
            t8.Layout.Tile = 24;
            axis off
        
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz(1,:), 'onerow', 'noverbose', 'existing_axes',t5);
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz(3,:), 'onerow', 'noverbose', 'existing_axes',t6);
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz(5,:), 'onerow', 'noverbose', 'existing_axes',t7);
            [o2, dat] = montage(o2, 'axial', 'wh_slice', grayord_xyz(2,:), 'onerow', 'noverbose', 'existing_axes',t8);
        
            title(t5,'X = -24','FontSize',10)
            title(t6,'X = -6','FontSize',10)
            title(t7,'X = 12','FontSize',10)
            title(t8, 'Z = -34','FontSize',10)
            for t = [t5,t6,t7]
                set(t,'XLim',[-100,30],'YLim',[-70,25]);
            end
            set(t8,'XLim',[-60,60],'YLim',[-100,30]);
        
            wh_surfaces = [1:4];

        case 'leftright inout subcortex'
            % saggital
            f1 = gcf;
            mainLayout = tiledlayout(f1,1,5);
            surfLayout = tiledlayout(mainLayout, 2, 2, 'Parent', mainLayout, 'TileSpacing', 'compact', 'Padding', 'none');
            surfLayout.Layout.Tile = 1;
            surfLayout.Layout.TileSpan = [1, 2];

            ax = {};
            for j = 1:4, ax{j} = nexttile(surfLayout); end
            o2 = surface(o2, 'axes', ax{1}, 'direction', 'bigbrain left');          
            o2 = surface(o2, 'axes', ax{2}, 'direction', 'bigbrain right');
            o2 = surface(o2, 'axes', ax{3}, 'direction', 'right_cutaway');
            o2 = surface(o2, 'axes', ax{4}, 'direction', 'left_cutaway');

            n_col = size(grayord_xyz,1);
            n_row = size(grayord_xyz,2);
            volLayout = tiledlayout(mainLayout, n_row, n_col, 'Parent', mainLayout, 'TileSpacing', 'tight', 'Padding', 'none');
            volLayout.Layout.Tile = 3;
            volLayout.Layout.TileSpan = [1, 3];
            ax_vol=[];
            for j = 1:n_row, for k = 1:n_col, ax_vol(j,k) = nexttile(volLayout); end; end
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes',ax_vol(1,:));
            for j=1:n_col
                set(ax_vol(1,j),'XLim',[-100,30],'YLim',[-70,25]);
                title(ax_vol(1,j), sprintf('x=%.0f', grayord_xyz(j,1)), 'FontSize', 10);
            end

            o2 = montage(o2, 'volume_data', dat, 'coronal', 'wh_slice', grayord_xyz, 'onerow','noverbose', 'existing_axes', ax_vol(2,:));
            for j=1:n_col
                set(ax_vol(2,j),'XLim',[-40,40],'YLim',[-70,25]);
                title(ax_vol(2,j), sprintf('y=%.0f', grayord_xyz(j,2)), 'FontSize', 10);
            end

            o2 = montage(o2, 'volume_data', dat, 'axial', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes', ax_vol(3,:));
            for j=1:n_col
                set(ax_vol(3,j),'XLim',[-60,60],'YLim',[-100,30]);
                title(ax_vol(3,j), sprintf('z=%.0f', grayord_xyz(j,3)), 'FontSize', 10);
            end

            wh_montages = [1 2 3];
            wh_surfaces = [1:4];


        case 'subcortex full'
            % saggital
            f1 = gcf;
            mainLayout = tiledlayout(f1,1,5);
            surfLayout = tiledlayout(mainLayout, 2, 2, 'Parent', mainLayout, 'TileSpacing', 'compact', 'Padding', 'none');
            surfLayout.Layout.Tile = 1;
            surfLayout.Layout.TileSpan = [1, 2];

            ax = {};
            for j = 1:4, ax{j} = nexttile(surfLayout); end
            o2 = surface(o2, 'axes', ax{1}, 'direction', 'caudate left');          
            o2 = surface(o2, 'axes', ax{2}, 'direction', 'caudate right');
            o2 = surface(o2, 'axes', ax{3}, 'direction', 'brainstem left');
            o2 = surface(o2, 'axes', ax{4}, 'direction', 'brainstem right');

            n_col = size(grayord_xyz,1);
            n_row = size(grayord_xyz,2);
            volLayout = tiledlayout(mainLayout, n_row, n_col, 'Parent', mainLayout, 'TileSpacing', 'tight', 'Padding', 'none');
            volLayout.Layout.Tile = 3;
            volLayout.Layout.TileSpan = [1, 3];
            ax_vol=[];
            for j = 1:n_row, for k = 1:n_col, ax_vol(j,k) = nexttile(volLayout); end; end
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes',ax_vol(1,:));
            for j=1:n_col
                set(ax_vol(1,j),'XLim',[-100,30],'YLim',[-70,25]);
                title(ax_vol(1,j), sprintf('x=%.0f', grayord_xyz(j,1)), 'FontSize', 10);
            end

            o2 = montage(o2, 'volume_data', dat, 'coronal', 'wh_slice', grayord_xyz, 'onerow','noverbose', 'existing_axes', ax_vol(2,:));
            for j=1:n_col
                set(ax_vol(2,j),'XLim',[-40,40],'YLim',[-70,25]);
                title(ax_vol(2,j), sprintf('y=%.0f', grayord_xyz(j,2)), 'FontSize', 10);
            end

            o2 = montage(o2, 'volume_data', dat, 'axial', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes', ax_vol(3,:));
            for j=1:n_col
                set(ax_vol(3,j),'XLim',[-60,60],'YLim',[-100,30]);
                title(ax_vol(3,j), sprintf('z=%.0f', grayord_xyz(j,3)), 'FontSize', 10);
            end

            wh_montages = [1 2 3];
            wh_surfaces = [1:4];


        case 'subcortex slices'
            f1 = gcf;
            mainLayout = tiledlayout(f1,1,5);

            n_col = size(grayord_xyz,1);
            n_row = size(grayord_xyz,2);
            volLayout = tiledlayout(mainLayout, n_row, n_col, 'Parent', mainLayout, 'TileSpacing', 'tight', 'Padding', 'none');
            volLayout.Layout.Tile = 3;
            volLayout.Layout.TileSpan = [1, 3];
            ax_vol=[];
            for j = 1:n_row, for k = 1:n_col, ax_vol(j,k) = nexttile(volLayout); end; end
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes',ax_vol(1,:));
            for j=1:n_col
                set(ax_vol(1,j),'XLim',[-100,30],'YLim',[-70,25]);
                title(ax_vol(1,j), sprintf('x=%.0f', grayord_xyz(j,1)), 'FontSize', 10);
            end

            o2 = montage(o2, 'volume_data', dat, 'coronal', 'wh_slice', grayord_xyz, 'onerow','noverbose', 'existing_axes', ax_vol(2,:));
            for j=1:n_col
                set(ax_vol(2,j),'XLim',[-40,40],'YLim',[-70,25]);
                title(ax_vol(2,j), sprintf('y=%.0f', grayord_xyz(j,2)), 'FontSize', 10);
            end

            o2 = montage(o2, 'volume_data', dat, 'axial', 'wh_slice', grayord_xyz, 'onerow', 'noverbose', 'existing_axes', ax_vol(3,:));
            for j=1:n_col
                set(ax_vol(3,j),'XLim',[-60,60],'YLim',[-100,30]);
                title(ax_vol(3,j), sprintf('z=%.0f', grayord_xyz(j,3)), 'FontSize', 10);
            end

            wh_montages = [1 2 3];
            wh_surfaces = [1:4];

        case 'hcp inflated'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'hcp inflated right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'hcp inflated left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'hcp inflated right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'hcp inflated left', 'orientation', 'lateral');
            
            wh_surfaces = [1:4];

        case 'hcp'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'surface right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'surface left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'surface right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'surface left', 'orientation', 'lateral');

            wh_surfaces = [1:4];

        case 'hcp sphere'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'hcp sphere right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'hcp sphere left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'hcp sphere right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'hcp sphere left', 'orientation', 'lateral');
            
            wh_surfaces = [1:4];

        case 'freesurfer inflated'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'freesurfer inflated right', 'orientation', 'medial', 'targetsurface', 'fsaverage_164k');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'freesurfer inflated left', 'orientation', 'medial', 'targetsurface', 'fsaverage_164k');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'freesurfer inflated right', 'orientation', 'lateral', 'targetsurface', 'fsaverage_164k');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'freesurfer inflated left', 'orientation', 'lateral', 'targetsurface', 'fsaverage_164k');
            
            wh_surfaces = [1:4];

        case 'freesurfer white'
            figure;
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'freesurfer white right', 'orientation', 'medial', 'targetsurface', 'fsaverage_164k');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'freesurfer white left', 'orientation', 'medial', 'targetsurface', 'fsaverage_164k');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'freesurfer white right', 'orientation', 'lateral', 'targetsurface', 'fsaverage_164k');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'freesurfer white left', 'orientation', 'lateral', 'targetsurface', 'fsaverage_164k');
            
            wh_surfaces = [1:4];

        case 'freesurfer sphere'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'freesurfer sphere right', 'orientation', 'medial', 'targetsurface', 'fsaverage_164k');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'freesurfer sphere left', 'orientation', 'medial', 'targetsurface', 'fsaverage_164k');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'freesurfer sphere right', 'orientation', 'lateral', 'targetsurface', 'fsaverage_164k');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'freesurfer sphere left', 'orientation', 'lateral', 'targetsurface', 'fsaverage_164k');
            
            wh_surfaces = [1:4];

        case 'MNI152NLin2009cAsym white'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'MNI152NLin2009cAsym white right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'MNI152NLin2009cAsym white left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'MNI152NLin2009cAsym white right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'MNI152NLin2009cAsym white left', 'orientation', 'lateral');
            
            wh_surfaces = [1:4];

        case 'MNI152NLin2009cAsym midthickness'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'MNI152NLin2009cAsym midthickness right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'MNI152NLin2009cAsym midthickness left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'MNI152NLin2009cAsym midthickness right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'MNI152NLin2009cAsym midthickness left', 'orientation', 'lateral');
            
            wh_surfaces = [1:4];

        case 'MNI152NLin2009cAsym pial'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'MNI152NLin2009cAsym pial right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'MNI152NLin2009cAsym pial left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'MNI152NLin2009cAsym pial right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'MNI152NLin2009cAsym pial left', 'orientation', 'lateral');
            
            wh_surfaces = [1:4];

        case 'MNI152NLin6Asym white'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'MNI152NLin6Asym white right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'MNI152NLin6Asym white left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'MNI152NLin6Asym white right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'MNI152NLin6Asym white left', 'orientation', 'lateral');
            
            wh_surfaces = [1:4];

        case 'MNI152NLin6Asym midthickness'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'MNI152NLin6Asym midthickness right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'MNI152NLin6Asym midthickness left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'MNI152NLin6Asym midthickness right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'MNI152NLin6Asym midthickness left', 'orientation', 'lateral');
            
            wh_surfaces = [1:4];

        case 'MNI152NLin6Asym pial'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'MNI152NLin6Asym pial right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'MNI152NLin6Asym pial left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'MNI152NLin6Asym pial right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'MNI152NLin6Asym pial left', 'orientation', 'lateral');
            
            wh_surfaces = [1:4];

        case 'MNI152NLin6Asym sphere'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'MNI152NLin6Asym sphere right', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'MNI152NLin6Asym sphere left', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'MNI152NLin6Asym sphere right', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'MNI152NLin6Asym sphere left', 'orientation', 'lateral');
            
            wh_surfaces = [1:4];

        case 'leftright inout'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'bigbrain left');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'bigbrain right');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'right_cutaway');
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'left_cutaway');
            
            wh_surfaces = [1:4];

        case 'subcortex 3d'
            axis off;
            o2 = surface(o2, 'axes', [0 0.5 .45 .45], 'direction', 'caudate left');
            o2 = surface(o2, 'axes', [0.4 0.5 .45 .45], 'direction', 'caudate right');
            o2 = surface(o2, 'axes', [0 0 .45 .45], 'direction', 'brainstem left');
            o2 = surface(o2, 'axes', [0.4 0 .45 .45], 'direction', 'brainstem right');

            wh_surfaces = [1:4];

        otherwise 
            error('illegal montage type. choose one of the following: blobcenters, regioncenters, compact, multirow, full, full2, compact2, compact3, coronal, saggital, allslices, full hcp, full hcp inflated, full no surfaces, hcp grayordinates, hcp grayordinates subcortex, hcp inflated, hcp sphere, freesurfer inflated, freesurfer white, freesurfer sphere, MNI152NLin2009cAsym white, MNI152NLin2009cAsym midthickness, MNI152NLin2009cAsym pial, MNI152NLin6Asym white, MNI152NLin6Asym midthickness, MNI152NLin6Asym pial, MNI152NLin6Asym sphere, inout leftright, inout leftright subcortex.');
    end
    
    % wh_montages = [1 2];

else % use existing o2 object to add montages

    if doverbose, disp('Using existing fmridisplay object'); end
    
    % Other inputs will be passed into addblobs
    existingmons = length(o2.montage);
    existingsurs = length(o2.surface);
    
    if doaddmontages
        % use same o2, but add montages
        switch montagetype
            case 'full'
            % saggital
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow', 'noverbose');
            shift_axes(-0.02, -0.04);        

            % coronal
            axh = axes('Position', [-0.02 0.37 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'coronal', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);

            % axial
            axh = axes('Position', [-0.02 0.19 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);

            axh = axes('Position', [-0.02 0.01 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-44 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            allaxh = findobj(gcf, 'Type', 'axes');
            disp(length(allaxh));
            for i = 1:(length(allaxh)-36)
                pos1 = get(allaxh(i), 'Position');
                pos1(1) = pos1(1) - 0.03;
                set(allaxh(i), 'Position', pos1);
            end

            % --- Coordinate Annotation Overlay ---
            if coordinates
                for m = 1:numel(o2.montage)
                    this_montage = o2.montage{m};
                    if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                        continue;
                    end
            
                    axlist = this_montage.axis_handles;
                    coords = this_montage.slice_mm_coords;
            
                    switch this_montage.orientation
                        case {'saggital', 'sagittal'}, label_prefix = 'x=';
                        case 'coronal',                label_prefix = 'y=';
                        case 'axial',                  label_prefix = 'z=';
                        otherwise,                     label_prefix = '';
                    end
            
                    for a = 1:min(numel(axlist), numel(coords))
                        ax = axlist(a);
                        if ~isgraphics(ax), continue; end
            
                        pos = get(ax, 'Position');
                        ann_x = pos(1) + pos(3)/2 - 0.025;
            
                        % 💡 Use slightly lower label height for sagittal to stay under surface views
                        if any(strcmp(this_montage.orientation, {'saggital', 'sagittal'}))
                            ann_y = pos(2) + pos(4) - 0.5;
                        else
                            ann_y = pos(2) + pos(4) - 0.01;
                        end
            
                        % Clamp to visible range
                        if ann_x < 0 || ann_x > 1 || ann_y < 0 || ann_y > 1
                            continue;
                        end
            
                        annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                            'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                            'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                            'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'VerticalAlignment', 'bottom');
                    end
                end
            end

            % surface
            o2 = surface(o2, 'axes', [0.1 0.74 .25 .25], 'direction', 'hires left', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.3 0.74 .25 .25], 'direction', 'hires right', 'orientation', 'medial');            
            o2 = surface(o2, 'axes', [0.5 0.74 .25 .25], 'direction', 'hires left', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.7 0.74 .25 .25], 'direction', 'hires right', 'orientation', 'lateral');

            wh_surfaces = existingsurs + [1 2 3 4];

            case 'compact'
                [o2, dat] = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 6, 'noverbose');
                axh = axes('Position', [0.05 0.4 .1 .5]);
                o2 = montage(o2, 'volume_data', dat, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh, 'noverbose');

                % --- Apply coordinate titles if coordinates flag is on ---
                if coordinates
                    label_map = containers.Map({'sagittal', 'coronal', 'axial'}, {'x', 'y', 'z'});
                
                    for m = 1:numel(o2.montage)
                        mon = o2.montage{m};
                        if ~isfield(mon, 'axis_handles') || isempty(mon.axis_handles)
                            continue;
                        end
                
                        if ~isfield(mon, 'orientation') || ~isKey(label_map, mon.orientation)
                            continue;
                        end
                        label = label_map(mon.orientation);
                
                        for i = 1:numel(mon.axis_handles)
                            ax = mon.axis_handles(i);
                            if isgraphics(ax)
                                mm_coord = getappdata(ax, 'mm_coord');
                                if isempty(mm_coord), continue; end
                
                                % Get axis position in normalized figure units
                                axpos = get(ax, 'Position');  % [x y w h]
                                ann_x = axpos(1) + axpos(3) * 0.5;   % horizontal center
                                ann_y = axpos(2) + axpos(4) + 0.01;  % just above axis
                
                                % Ensure we're in bounds
                                if ann_x >= 0 && ann_x <= 1 && ann_y >= 0 && ann_y <= 1
                                    annotation(gcf, 'textbox', [ann_x - 0.05, ann_y, 0.1, 0.03], ...
                                        'String', sprintf('%s=%.0f', label, mm_coord), ...
                                        'HorizontalAlignment', 'center', ...
                                        'VerticalAlignment', 'bottom', ...
                                        'EdgeColor', 'none', ...
                                        'FontSize', 14);
                                end
                            end
                        end
                    end
                end

            case 'compact2'
                [o2, dat] = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose');
                enlarge_axes(gcf, 1);
                axh = axes('Position', [-0.03 0.15 .2 1]);
                o2 = montage(o2, 'volume_data', dat, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh, 'noverbose');

                % shift all axes down and right
                shift_axes(+0.03, -0.10);

                %ss = get(0, 'ScreenSize');
                %set(gcf, 'Position', [round(ss(3)/12) round(ss(4)*.9) round(ss(3)*.8) round(ss(4)/5.5) ])

    
                % --- Coordinate Annotation Overlay ---
                if coordinates
                    for m = 1:numel(o2.montage)
                        this_montage = o2.montage{m};
    
                        if ~isfield(this_montage, 'axis_handles') || ~isfield(this_montage, 'slice_mm_coords')
                            continue;
                        end
    
                        axlist = this_montage.axis_handles;
                        coords = this_montage.slice_mm_coords;
    
                        switch this_montage.orientation
                            case {'saggital', 'sagittal'}, label_prefix = 'x=';
                            case 'coronal',  label_prefix = 'y=';
                            case 'axial',    label_prefix = 'z=';
                            otherwise,       label_prefix = '';
                        end
    
                        for a = 1:min(numel(axlist), numel(coords))
                            ax = axlist(a);
                            if ~isgraphics(ax), continue; end
    
                            pos = get(ax, 'Position');
                            ann_x = pos(1) + pos(3)/2 - 0.025;
                            ann_y = pos(2) + pos(4) - 0.40;
    
                            % Adjust position for sagittal slices to avoid surface overlap
                            if any(strcmp(this_montage.orientation, {'saggital', 'sagittal'}))
                                ann_y = pos(2) + pos(4)  - 0.05;
                            end
    
    
                            % Skip if outside [0, 1] bounds (annotation cannot handle this)
                            if any([ann_x, ann_y, ann_x + 0.05, ann_y + 0.03] > 1) || ...
                               any([ann_x, ann_y] < 0)
                                continue;
                            end
    
                            annotation(get(ax, 'Parent'), 'textbox', [ann_x, ann_y, 0.05, 0.03], ...
                                'String', sprintf('%s%.0f', label_prefix, coords(a)), ...
                                'FitBoxToText', 'on', 'EdgeColor', 'none', ...
                                'HorizontalAlignment', 'center', 'FontSize', 10, ...
                                'VerticalAlignment', 'bottom');
                        end
                    end
                end                     
                    
            otherwise error('illegal montage type. choose full, compact, or compact2 when adding to existing montage set.')
        end
        
        wh_montages = existingmons + [1 2];

    else
        if doremove
            o2 = removeblobs(o2);
        end
        
        wh_montages = 1:existingmons;
        wh_surfaces = 1:existingsurs;

    end
    
end

% Safely extract axis handles
% sag_axes = getfield(o2.montage{1}, 'axis_handles', {});
% cor_axes = getfield(o2.montage{2}, 'axis_handles', {});
% axi_axes = getfield(o2.montage{3}, 'axis_handles', {});
% 
% % Determine how many titles we can safely assign
% n_sag = min(numel(sag_axes), size(xyz, 1));
% n_cor = min(numel(cor_axes), size(xyz, 1));
% n_axi = min(numel(axi_axes), size(xyz, 1));
% 
% % Apply titles to each montage view
% for i = 1:n_sag
%     title(sag_axes(i), sprintf('x=%.0f', xyz(i,1)), 'FontSize', 10);
% end
% for i = 1:n_cor
%     title(cor_axes(i), sprintf('y=%.0f', xyz(i,2)), 'FontSize', 10);
% end
% for i = 1:n_axi
%     title(axi_axes(i), sprintf('z=%.0f', xyz(i,3)), 'FontSize', 10);
% end


% if coordinates && exist('xyz', 'var') && ~isempty(xyz)
%     n_sag = numel(sag_axes);
%     n_cor = numel(cor_axes);
%     n_axi = numel(axi_axes);
% 
%     for i = 1:n_sag
%         title(sag_axes(i), sprintf('x=%.0f', xyz(i,1)), 'FontSize', 10);
%     end
%     for i = 1:n_cor
%         title(cor_axes(i), sprintf('y=%.0f', xyz(i,2)), 'FontSize', 10);
%     end
%     for i = 1:n_axi
%         title(axi_axes(i), sprintf('z=%.0f', xyz(i,3)), 'FontSize', 10);
%     end
% end

% if coordinates && isfield(o2, 'montage')
%     coord_labels = {'x', 'y', 'z', 'z'};  % label for each montage
% 
%     for m = 1:numel(o2.montage)
%         if ~isfield(o2.montage{m}, 'axis_handles'), continue; end
% 
%         axlist = o2.montage{m}.axis_handles;
%         coords = slice_coords{m};
%         lab = coord_labels{m};
% 
%         for i = 1:min(numel(axlist), numel(coords))
%             title(axlist(i), sprintf('%s=%.0f', lab, coords(i)), 'FontSize', 10);
%         end
%     end
% end






% % Only apply coordinate titles if 'coordinates' is specified
% allax = findobj(gcf, 'Type', 'axes');
% % Sort by vertical position (bottom to top) and then flip for montage order
% [~, sortIdx] = sort(arrayfun(@(h) h.Position(2), allax));
% ax_sorted = allax(sortIdx);
% ax_sorted = flipud(ax_sorted);  % match montage order
% 
% if coordinates
%     num_coords = size(xyz, 1);
%     for i = 1:min(num_coords, numel(ax_sorted))
%         ax = ax_sorted(i);
%         try
%             if num_coords <= numel(axh)
%                 subplot_handle = axh(i); % assumes axh already sorted
%                 axes(subplot_handle);
%                 title(['x=' num2str(xyz(i,1)) ', y=' num2str(xyz(i,2)) ', z=' num2str(xyz(i,3))]);
%             end
%         catch
%             % Fail silently
%         end
%     end
% end


% Now we can add blobs
% --------------------------------------------

% they are added to all montages by default, but you can specify selected
% montages if you want to as well.

% it's easy to remove them as well:
% o2 = removeblobs(o2);

if doblobs
    if exist('wh_surfaces', 'var')
        o2 = addblobs(o2, cl, 'splitcolor', splitcolor, 'wh_montages', wh_montages, 'wh_surfaces', wh_surfaces, 'nolegend', varargin{:});
    else
        o2 = addblobs(o2, cl, 'splitcolor', splitcolor, 'wh_montages', wh_montages, varargin{:});
    end
end

if dooutline
    o2 = addblobs(o2, cl, 'color', outlinecolor, 'outline', 'wh_montages', wh_montages, 'no_surface');
end


% ------------------------------------------------------
% INLINE FUNCTIONS
% ------------------------------------------------------

    function shift_axes(x_offset, y_offset)

        % usage:   function shift_axes(-0.03, -0.18);
        % purpose: To shift axes according of the habdles.
        %
        % input:   Function parameters are (x, y) offset relative to current position

        % shift all axes according to offset values
        allaxh = findobj(gcf, 'Type', 'axes');
        for i = 1:length(allaxh)
            pos1 = get(allaxh(i), 'Position');
            pos1(2) = pos1(2) + y_offset;
            pos1(1) = pos1(1) + x_offset;
            set(allaxh(i), 'Position', pos1);
        end
    end % shift_axes

    function [xn, yn] = ds2nfu(ax, xd, yd)
        % Helper function for 'coordinates' argument - Michael Sun, PhD
        % 06/09/2025

        % Convert data space (xd, yd) into normalized figure units (xn, yn)
        axpos = get(ax, 'Position');
        xlim = get(ax, 'XLim');
        ylim = get(ax, 'YLim');
    
        xn = axpos(1) + (xd - xlim(1)) / diff(xlim) * axpos(3);
        yn = axpos(2) + (yd - ylim(1)) / diff(ylim) * axpos(4);
    end



end  % function
