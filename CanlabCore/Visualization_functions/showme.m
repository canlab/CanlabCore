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
    %    [o1, o2, o3, 04, roi] = lookat('insula');
    %    [o1, o2, o3, 04, roi] = lookat('insula', 'showSurfaces', true, 'inflation', 'inflated', 'showVolumes', true, 'volRadius', 10, 'showFlatmap', true, 'showSlabs', true, 'sourcespace', 'MNI152NLin2009cAsym', 'atlas', load_atlas('canlab2024'));
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
    addParameter(parser, 'showFlatmap', true);
    addParameter(parser, 'showSlabs', true);
    addParameter(parser, 'sourcespace', 'MNI152NLin2009cAsym');
    addParameter(parser, 'atlas', load_atlas('canlab2024'));
    parse(parser, varargin{:});
    showSurfaces = parser.Results.showSurfaces;
    inflation = parser.Results.inflation;
    sourcespace = parser.Results.sourcespace;
    
    showVolumes = parser.Results.showVolumes;
    volRadius = parser.Results.volRadius;

    showFlatmap = parser.Results.showFlatmap;
    showSlabs = parser.Results.showSlabs;
    
    
    atl = parser.Results.atlas;

    img=[];
    % pain_pathways=load_atlas('painpathways')
    % atl=pain_pathways;
    % atl=[];

    % If Region


    % If string/text:

    if isempty(img)
        img=fmri_data;
        img.dat=zeros(numel(img.dat),1);
    end

    % if isempty(atl)
    %     atl=load_atlas('canlab2024');
    % end

    % tbl(contains(tbl.Network, net_word), :);
    % reg = atl.select_atlas_subset(lbl).threshold(0.2);

    disp(lbl)
    lbl=cellstr(lbl);
    
    labels = {'labels', 'labels_2', 'labels_3', 'labels_4', 'labels_5', 'label_descriptions'};
    roi = [];
    
    for i = 1:length(labels)
        try
            roi = atl.select_atlas_subset(lbl, labels{i});
            break; % Exit the loop if successful
        catch
            if i == length(labels)
                error('No regions identified in atlas based on search term');
            end
            warning('No regions identified in %s, trying next', labels{i});
        end
    end


    
    mm_center = atlas2region(roi).mm_center;
    mm_center = [mm_center-volRadius; mm_center-(volRadius/2); mm_center; mm_center+(volRadius/2); mm_center+volRadius];
    

    
    figure;
    % Show surfaces
    if showSurfaces
        o1 = fmridisplay();
        [~, axh] = create_figure('volumetric_slices', 1, 5, false, true);
        
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
        
        % colors = scn_standard_colors(2);
        num_regions=numel(roi.labels);
        o1=addblobs(o1,roi, 'indexmap', lines(num_regions), 'trans', 'transvalue', .5, 'interp', 'nearest');

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
        num_regions=numel(roi.labels);
        o2=addblobs(o2, roi, 'indexmap', lines(num_regions), 'trans', 'transvalue', .5, 'interp', 'nearest');

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