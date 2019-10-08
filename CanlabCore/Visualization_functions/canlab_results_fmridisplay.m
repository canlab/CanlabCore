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
%   **'nooutline':**
%        do not display blob outlines
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
%        'full' for full montages of axial and sagg slices.
%
%        'full hcp' for full montage, but with surfaces and volumes from
%        HCP data
%
%        'compact' [default] for single-figure parasagittal and axials slices.
%
%        'compact2': like 'compact', but fewer axial slices.
%
%        'multirow': followed by number of rows
%           e.g., o2 = canlab_results_fmridisplay([], 'multirow', 2);
%
%        {'blobcenters', 'regioncenters'}: Slices for the center of each
%        blob/region
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
%         The default brain for overlays is based on Keuken et al. 2014
%         For legacy SPM8 single subject, enter as arguments:
%         'overlay', which('SPM8_colin27T1_seg.img')
% 
% Other inputs to addblobs (fmridisplay method) are allowed, e.g., 'cmaprange', [-2 2], 'trans'
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
% ..
%    Tor Wager
%    1/27/2012
% ..

if ~which('fmridisplay.m')
    disp('fmridisplay is not on path.  it is in canlab tools, which must be on your path!')
    return
end

if nargin == 0
    o2 = canlab_results_fmridisplay(region(), 'noblobs', 'nooutline');
    return
end

if ischar(input_activation)
    
    if strcmp(input_activation, 'compact') || strcmp(input_activation, 'compact2') || strcmp(input_activation, 'full') ...
            || strcmp(input_activation, 'multirow') || strcmp(input_activation, 'coronal') || strcmp(input_activation, 'sagittal')
        
        % Entered no data map; intention is not to plot blobs, just create underlay
        varargin{end + 1} = 'noblobs'; 
        varargin{end + 1} = 'nooutline';
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
dooutline = true;
doaddmontages = false;
doremove = true;
outlinecolor = [0 0 0];
splitcolor = {[0 0 1] [0 .8 .8] [1 .4 .5] [1 1 0]}; % {[0 0 1] [.3 .6 .9] [.8 .3 0] [1 1 0]};  % more straight orange to yellow: {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}
montagetype = 'compact';
doverbose = true;
%overlay='SPM8_colin27T1_seg.img';
overlay = 'keuken_2014_enhanced_for_underlay.img';
dofigure = true;

wh = strcmp(varargin, 'overlay');
if any(wh), wh = find(wh); overlay = varargin{wh(1) + 1};  varargin([wh wh+1]) = []; end

wh = strcmp(varargin, 'noblobs');
if any(wh), doblobs = false; varargin(wh) = []; end

wh = strcmp(varargin, 'nooutline');
if any(wh), dooutline = false; varargin(wh) = []; end

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

wh = strcmp(varargin, 'blobcenters');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'regioncenters');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'full hcp');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'compact');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'compact2');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'coronal');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'saggital');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'allslices');
if any(wh), montagetype = varargin{find(wh)}; varargin(wh) = []; end

wh = strcmp(varargin, 'noverbose');
if any(wh), doverbose = false; varargin(wh) = []; end

wh = strcmp(varargin, 'nofigure');
if any(wh), dofigure = false; varargin(wh) = []; end

wh = false(1, length(varargin));
for i = 1:length(varargin)
    wh(i) = isa(varargin{i}, 'fmridisplay');
    if wh(i), o2 = varargin{wh}; end
end
varargin(wh) = [];

xyz = [-20 -10 -6 -2 0 2 6 10 20]';
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
        disp('This takes a lot of memory, and can hang if you have too little.');
        
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
            %Make a series of montages at center of each region and add blobs to that:
            
            if ~exist('cl', 'var'), error('Must enter region object to use blobcenters option.'); end
            
            xyz = cat(1, cl.mm_center);
            
            %           onerowstr = [];
            %           if length(cl) < 20, onerowstr = 'onerow'; end
            
            orientation = 'axial';
            
            % Make a grid - determine subplots
            nr = floor(sqrt(length(cl)));
            nc = ceil(length(cl) ./ nr);
            
            [~, axh] = create_figure('fmridisplay_regioncenters', nr, nc, false, true); 
            
            set(axh,'Visible','off'); % turn off axis grid for all axes
             
            for i = 1:length(cl)
                
                %axh(i) = subplot(nr, nc, i);
                %axes(axh(i))
                
                if i == 1
                    [o2, dat] = montage(o2, orientation, 'wh_slice', xyz(i, :), 'onerow', 'existing_axes', axh(i), 'existing_figure', 'noverbose');
                else
                    o2 = montage(o2, 'volume_data', dat, orientation, 'wh_slice', xyz(i, :), 'onerow', 'existing_axes', axh(i), 'existing_figure', 'noverbose');
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

            % surface
            o2 = surface(o2, 'axes', [0.1 0.74 .25 .25], 'direction', 'hires left', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.3 0.74 .25 .25], 'direction', 'hires right', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.5 0.74 .25 .25], 'direction', 'hires left', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.7 0.74 .25 .25], 'direction', 'hires right', 'orientation', 'lateral');

            wh_montages = [1 2 3 4];
            wh_surfaces = [1 2 3 4];


        case 'compact2'  % creates a new figure
            %subplot(2, 1, 1);
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
            %axh = axes('Position', [0.0 0.08 .15 1]);
            axh = axes('Position', [-0.02 .75-.3 .17 .17]);  % [-0.02 0.15+shiftvals(i) .17 .17]);
            axh(2) = axes('Position', [.022 .854-.3 .17 .17]);
            
            %o2 = montage(o2, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh, 'noverbose');
            o2 = montage(o2, 'volume_data', dat, 'saggital', 'slice_range', [-2 2], 'spacing', 4, 'onerow', 'noverbose', 'existing_axes', axh);

            %ss = get(0, 'ScreenSize');
            %set(gcf, 'Position', [round(ss(3)/12) round(ss(4)*.9) round(ss(3)*.8) round(ss(4)/5.5) ]) % this line messes p the
            %images, makes it too big an overlapping
            wh_montages = [1 2];
            
            brighten(.4)
            
        case 'coronal'
            % coronal
            o2 = montage(o2, 'coronal', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose');
             wh_montages = 1;

        case  'saggital'
            o2 = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow', 'noverbose');
            
             wh_montages = 1;

        case 'allslices'
            
            [o2, dat] = montage(o2, 'saggital', 'wh_slice', xyz, 'onerow', 'noverbose');
            shift_axes(-0.02, -0.04);
            
            axh = axes('Position', [-0.02 0.37 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'coronal', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            % axial
            axh = axes('Position', [-0.02 0.19 .17 .17]);
            o2 = montage(o2, 'volume_data', dat, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose', 'existing_axes', axh);
            
            wh_montages = [1 2 3];
            
        case 'full hcp'
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

            % surface
            o2 = surface(o2, 'axes', [0.1 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'medial');
            o2 = surface(o2, 'axes', [0.3 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'medial');          
            o2 = surface(o2, 'axes', [0.5 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'lateral');
            o2 = surface(o2, 'axes', [0.7 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'lateral');
            
            wh_montages = [1 2 3 4];
            wh_surfaces = [1:8];

        otherwise error('illegal montage type. choose full or compact.');
    end
    
    % wh_montages = [1 2];

else
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
                
            case 'compact2'
                [o2, dat] = montage(o2, 'axial', 'slice_range', [-40 50], 'onerow', 'spacing', 8, 'noverbose');
                enlarge_axes(gcf, 1);
                axh = axes('Position', [-0.03 0.15 .2 1]);
                o2 = montage(o2, 'volume_data', dat, 'saggital', 'wh_slice', [0 0 0], 'existing_axes', axh, 'noverbose');

                % shift all axes down and right
                shift_axes(+0.03, -0.10);

                %ss = get(0, 'ScreenSize');
                %set(gcf, 'Position', [round(ss(3)/12) round(ss(4)*.9) round(ss(3)*.8) round(ss(4)/5.5) ])
                
                
            otherwise error('illegal montage type. choose full or compact when adding to existing montage set.')
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

% Now we can add blobs
% --------------------------------------------

% they are added to all montages by default, but you can specify selected
% montages if you want to as well.

% it's easy to remove them as well:
% o2 = removeblobs(o2);

if doblobs
    if exist('wh_surfaces', 'var')
        o2 = addblobs(o2, cl, 'splitcolor', splitcolor, 'wh_montages', wh_montages, 'wh_surfaces', wh_surfaces, varargin{:});
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


end  % function
