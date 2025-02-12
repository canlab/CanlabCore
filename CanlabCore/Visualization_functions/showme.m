function [o1, o2, o3, o4 roi] = showme(lbl, varargin)
    % showme- Looks up and visualizes regions of interest on
    % surfaces, volumes, flatmaps, and coronal slabs. Searches for terms
    % through each level of an atlas's labels. If atlas not specified, will
    % default to canlab2024.
    %
    % Syntax:  [o1, o2, o3, o4, roi] = showme(lbl, varargin)
    %
    % Inputs:
    %    lbl - Cell array of region labels to be plotted
    %    varargin - Additional parameters for customization (e.g., showSurfaces, sourcespace, inflation, showVolumes, volRadius, showFlatmap, showSlabs, atlas)
    %        'showSurfaces' - Boolean flag to display surfaces (default: true)
    %        'labels' - Cell string array to label surfaces (default: {})
    %        'inflation' - Type of surface inflation ('none', 'inflated', 'sphere', default: 'inflated')
    %        'showVolumes' - Boolean flag to display volumetric slices (default: true)
    %        'volRadius' - Radius for volumetric display in mm (default: 8)
    %        'showFlatmap' - Boolean flag to display flatmap (default: true)
    %        'showSlabs' - Boolean flag to display coronal slabs (default: true)
    %        'sourcespace' - Source space for the surfaces (default: 'MNI152NLin2009cAsym')
    %        'atlas' - Atlas data structure (default: load_atlas('canlab2024'))
    %
    % Outputs:
    %    o1 - Handle to the surface display object
    %    o2 - Handle to the volumetric slices display object
    %    o3 - Handle to the flatmap 
    %    o4 - Handle to the flatmap 
    %    roi - Region of interest data object
    %
    % Example: 
    %    [o1, o2, o3, 04, roi] = showme('insula');
    %    [o1, o2, o3, 04, roi] = showme('insula', 'showSurfaces', true, 'inflation', 'inflated', 'showVolumes', true, 'volRadius', 10, 'showFlatmap', true, 'showSlabs', true, 'sourcespace', 'MNI152NLin2009cAsym', 'atlas', load_atlas('canlab2024'));
    %   
    %   % Then, mask with your own data
    %   tthr=threshold(output, .05, 'fdr');
    %   o1=addblobs(o1, apply_mask(tthr, roi)); % Checkout your activations on the surface 
    %   o2=addblobs(o2, apply_mask(tthr, roi)); % Checkout your activations on volumes
    %
    % See also: fmridisplay, create_figure, atlas2region, load_atlas
    %
    % Author: Michael Sun, Ph.D. 7/24/2024


    % Parse optional inputs
    parser = inputParser;
    addParameter(parser, 'showSurfaces', true);
    addParameter(parser, 'inflation', 'inflated');
    addParameter(parser, 'showVolumes', true);
    addParameter(parser, 'volRadius', 8);
    addParameter(parser, 'showFlatmap', false);
    addParameter(parser, 'showSlabs', false);
    addParameter(parser, 'sourcespace', 'MNI152NLin2009cAsym');
    addParameter(parser, 'atlas', load_atlas('canlab2024'));
    addParameter(parser, 'labels', {});
    addParameter(parser, 'colormap', []);
    addParameter(parser, 'noplot', false);
    parse(parser, varargin{:});
    showSurfaces = parser.Results.showSurfaces;
    inflation = parser.Results.inflation;
    sourcespace = parser.Results.sourcespace;
    
    showVolumes = parser.Results.showVolumes;
    volRadius = parser.Results.volRadius;

    showFlatmap = parser.Results.showFlatmap;
    showSlabs = parser.Results.showSlabs;
    
    
    atl = parser.Results.atlas;

    % atl.probability_maps=[]; % Show deterministic maps only.

    mylabels = format_strings_for_legend(parser.Results.labels);

    cmap = parser.Results.colormap;

    noplot = parser.Results.noplot;

    if noplot==true
        display('noplot is toggled. Plotting suppressed.')
        showSurfaces=false;
        showVolumes=false;
        showFlatmap=false;
        showSlabs=false;
    end

    [o1 o2 o3 o4] = deal([]); 

    img=[];
    if isa(lbl, 'region')
        roi=lbl;
    elseif isa(lbl, 'atlas')
        roi=lbl;

    else     % If string/text:
        disp(lbl)
        lbl=cellstr(lbl);
    
        labels = {'labels', 'labels_2', 'labels_3', 'labels_4', 'labels_5', 'label_descriptions'};
        roi = {};
        
        % for i = 1:length(labels)
        for r = 1:numel(lbl)
            for i = 1:length(labels)
                try
                    roi{end+1} = atlas2region(atl.select_atlas_subset(lbl(r), labels{i}, 'flatten'));
                    roi{end}.Z=repmat(r, 1, roi{end}.numVox);
                    break; % Exit the loop if successful
                catch
                    if i == length(labels)
                        error('No regions identified in atlas based on search term');
                    end
                    warning('No regions identified in %s, trying next', labels{i});
                end
            end
        end

        roi=[roi{:}];
    end

    if isempty(img)
        img=fmri_data;
        img.dat=zeros(numel(img.dat),1);
    end

    
    if isa(roi, 'atlas')
        num_regions=numel(roi.labels);
        mm_center = atlas2region(roi).mm_center;
    else
        mm_center = roi.mm_center;
        num_regions=numel(roi);
    end

    mm_center = [mm_center-volRadius; mm_center-(volRadius/2); mm_center; mm_center+(volRadius/2); mm_center+volRadius];

    if isempty(cmap)
        cmap = hsv(num_regions);
    end
    
    
    % Show surfaces
    if showSurfaces
        o1 = fmridisplay();
        [~, axh] = create_figure('surface_plots', 1, 5, false, true);
        
        % Retrieve all axes handles in the figure
        allAxesHandles = findall(gcf, 'Type', 'axes');
        
        % Sort axes handles based on their creation order (position order in the figure)
        [~, sortIdx] = sort(arrayfun(@(h) h.Position(2), allAxesHandles));
        sortedAxesHandles = allAxesHandles(sortIdx);
        
        % Match the sorted handles to axh
        actualAxHandles = sortedAxesHandles;
        actualAxHandles = flipud(actualAxHandles);
    
        % Turn off the axis lines, ticks, and labels for all axes
        set([actualAxHandles.XAxis], 'Visible', 'off');
        set([actualAxHandles.YAxis], 'Visible', 'off');

        switch inflation
            case 'none'

                o1 = surface(o1, 'axes', [0.1 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'medial', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.3 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'medial', 'sourcespace', sourcespace);          
                o1 = surface(o1, 'axes', [0.5 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'lateral', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.7 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'lateral', 'sourcespace', sourcespace);

            case 'inflated'
                o1 = surface(o1, 'axes', [0.1 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'medial', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.3 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'medial', 'sourcespace', sourcespace);          
                o1 = surface(o1, 'axes', [0.5 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'lateral', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.7 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'lateral', 'sourcespace', sourcespace);

                o1 = surface(o1, 'axes', [0.1 0.44 .25 .25], 'direction', 'hcp inflated right', 'orientation', 'medial', 'sourcespace', sourcespace);          
                o1 = surface(o1, 'axes', [0.3 0.44 .25 .25], 'direction', 'hcp inflated left', 'orientation', 'lateral', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.5 0.44 .25 .25], 'direction', 'hcp inflated left', 'orientation', 'medial', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.7 0.44 .25 .25], 'direction', 'hcp inflated right', 'orientation', 'lateral', 'sourcespace', sourcespace);
            case 'sphere'
                o1 = surface(o1, 'axes', [0.1 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'medial', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.3 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'medial', 'sourcespace', sourcespace);          
                o1 = surface(o1, 'axes', [0.5 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'lateral', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.7 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'lateral', 'sourcespace', sourcespace);
            
                % Show sphere instead..THIS CURRENTLY DOESN'T SEEM TO WORK
                o1 = surface(o1, 'axes', [0.1 0.44 .25 .25], 'direction', 'hcp sphere right', 'orientation', 'medial', 'sourcespace', sourcespace);          
                o1 = surface(o1, 'axes', [0.3 0.44 .25 .25], 'direction', 'hcp sphere left', 'orientation', 'lateral', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.5 0.44 .25 .25], 'direction', 'hcp sphere left', 'orientation', 'medial', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.7 0.44 .25 .25], 'direction', 'hcp sphere right', 'orientation', 'lateral', 'sourcespace', sourcespace);
            
            case 'leftright inout'
                o1 = canlab_results_fmridisplay([], 'montagetype', 'leftright inout');
                % render_on_surface(roi, o1);
            % case 'left_cutaway'
            %     o1 = canlab_results_fmridisplay([], 'montagetype', 'inout leftright');
            %     render_on_surface(roi, o1);
            otherwise

                o1 = surface(o1, 'axes', [0.1 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'medial', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.3 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'medial', 'sourcespace', sourcespace);          
                o1 = surface(o1, 'axes', [0.5 0.74 .25 .25], 'direction', 'surface left', 'orientation', 'lateral', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.7 0.74 .25 .25], 'direction', 'surface right', 'orientation', 'lateral', 'sourcespace', sourcespace);

                o1 = surface(o1, 'axes', [0.1 0.44 .25 .25], 'direction', 'hcp inflated right', 'orientation', 'medial', 'sourcespace', sourcespace);          
                o1 = surface(o1, 'axes', [0.3 0.44 .25 .25], 'direction', 'hcp inflated left', 'orientation', 'lateral', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.5 0.44 .25 .25], 'direction', 'hcp inflated left', 'orientation', 'medial', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.7 0.44 .25 .25], 'direction', 'hcp inflated right', 'orientation', 'lateral', 'sourcespace', sourcespace);
                
                % Show sphere instead..THIS CURRENTLY DOESN'T SEEM TO WORK
                o1 = surface(o1, 'axes', [0.1 0.14 .25 .25], 'direction', 'hcp sphere right', 'orientation', 'medial', 'sourcespace', sourcespace);          
                o1 = surface(o1, 'axes', [0.3 0.14 .25 .25], 'direction', 'hcp sphere left', 'orientation', 'lateral', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.5 0.14 .25 .25], 'direction', 'hcp sphere left', 'orientation', 'medial', 'sourcespace', sourcespace);
                o1 = surface(o1, 'axes', [0.7 0.14 .25 .25], 'direction', 'hcp sphere right', 'orientation', 'lateral', 'sourcespace', sourcespace);
        end
        
        o1=addblobs(o1, roi, 'indexmap', cmap, 'trans', 'transvalue', .5, 'interp', 'nearest');

        % colors = scn_standard_colors(2);

        % Calculate the YTick positions to be centered within each color segment
        x_positions = linspace(0, 1, numel(mylabels) + 1); % +1 for the edges
        x_positions = (x_positions(1:end-1) + x_positions(2:end)) / 2; % Midpoints
        
        bar1axis = axes('Position', [.35 -0.30 .38 .42]);
        % bar1axis = axes('Position', [0 0 0 0]);
        colormap(bar1axis, cmap)
        colorbar_han = colorbar(bar1axis, 'northoutside');
        set(bar1axis, 'Visible', 'off');
        colorbar_han.Ticks=x_positions;
        colorbar_han.TickLabels=mylabels;

        drawnow, snapnow;
    end

    if showVolumes
        % Show volumetric slices

        o2 = fmridisplay();
        [~, axh] = create_figure('volumetric_slices', 3, 5, false, true);
        
        % Retrieve all axes handles in the figure
        allAxesHandles = findall(gcf, 'Type', 'axes');
        
        % Sort axes handles based on their creation order (position order in the figure)
        [~, sortIdx] = sort(arrayfun(@(h) h.Position(2), allAxesHandles));
        sortedAxesHandles = allAxesHandles(sortIdx);
        
        % Match the sorted handles to axh
        actualAxHandles = sortedAxesHandles;
        actualAxHandles = flipud(actualAxHandles);
    
        % Turn off the axis lines, ticks, and labels for all axes
        set([actualAxHandles.XAxis], 'Visible', 'off');
        set([actualAxHandles.YAxis], 'Visible', 'off');

        o2 = montage(o2,'saggital','wh_slice', mm_center, 'existing_axes', axh(1:5), 'existing_figure');
        o2 = montage(o2,'coronal','wh_slice', mm_center,'existing_axes', axh(6:10), 'existing_figure');
        o2 = montage(o2,'axial','wh_slice', mm_center, 'existing_axes', axh(11:15), 'existing_figure');

        % colors = scn_standard_colors(2);
        % num_regions=numel(roi.labels);
        % o2=addblobs(o2, roi, 'indexmap', cmap, 'trans', 'transvalue', .5, 'interp', 'nearest');
        o2=addblobs(o2, roi, 'indexmap', cmap, 'interp', 'nearest');

        drawnow, snapnow;


        % Get the axis position and set the title at the bottom
        
        for h = 1:5
            sagtitleText = ['x=' num2str(mm_center(h,1))];
            cortitleText = ['y=' num2str(mm_center(h,2))];
            axititleText = ['z=' num2str(mm_center(h,3))];
    
            % Easiest, stupidest way to do this without screwing up the
            % positioning:
            subplot(3,5,h)
            title(sagtitleText);
            subplot(3,5,h+5)
            title(cortitleText);
            subplot(3,5,h+10)
            title(axititleText);
    
        end
    end

    % At the moment, flat maps and coronal slabs don't seem to work at all.
    if showFlatmap
        figure;
        axis off;
        addbrain('flat surfaces', 'noverbose');
        allsurfaceHandles = findall(gcf, 'Type', 'surface');
    
        o3=render_on_surface(roi, allsurfaceHandles);
        sgtitle('Rendering on Flat Maps not yet working -- DEBUG ME')
        drawnow, snapnow
    end

    if showSlabs

        figure;
        axis off;
        addbrain('coronal_slabs', 'noverbose');
        allsurfaceHandles = findall(gcf, 'Type', 'surface');
    
        o4=render_on_surface(roi, allsurfaceHandles);
        sgtitle('Rendering on Coronal Slabs not yet working -- DEBUG ME')
        drawnow, snapnow
    end


end