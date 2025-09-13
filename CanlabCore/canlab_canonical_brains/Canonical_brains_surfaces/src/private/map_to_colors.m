function [c_lh_colored, c_rh_colored] = map_to_colors(c_lh, c_rh, o2, varargin)
    map_function = @(c, x1, x2, y1, y2)  y1 + (c - x1) * (y2 - y1) ./ (x2 - x1);

    if ismember(varargin,'doindexmap')
        nvals = min([256, length(unique([c_lh;c_rh]))]);
        colors = scn_standard_colors(nvals);
        cm = cell2mat(colors(:));
        uniq_regions = unique([c_lh;c_rh]);
        uniq_regions(uniq_regions == 0) = [];
        n_regions = length(uniq_regions);
        cm = repmat(cm, ceil(n_regions/size(cm,1)),1);
        cm = [0.5,0.5,0.5; cm]; % add gray
        colormap(o2.surface{1}.axis_handles,cm)
        clim = [0,n_regions];

        wh_lh = c_lh == 0 | isnan(c_lh);  
        wh_rh = c_rh == 0 | isnan(c_rh);  

        surf = get(o2.surface{4}.axis_handles,'Children');
        c_gray = zeros(size(get(surf(2), 'Vertices'),1),1);

        c_lh(mod(c_lh,1) > 0) = 0;
        c_lh_colored = c_lh;

        c_rh(mod(c_rh,1) > 0) = 0;
        c_rh_colored = c_rh;
    else
    
        % map to colors
        clim = [prctile([c_rh; c_lh],10),prctile([c_rh; c_lh],90)];
        nvals = 256;
        minposcolor = [1 0.4 0.5];
        maxposcolor = [1 1 0];
        minnegcolor = [0 0 1];
        maxnegcolor = [0 0.8 0.8];
        pos_colormap = colormap_tor(minposcolor, maxposcolor);
        neg_colormap = colormap_tor(minnegcolor, maxnegcolor);
        
        [cm, kpos, kneg] = hotcool_split_colormap(nvals, clim, o2.surface{1}.axis_handles, pos_colormap(1, :), pos_colormap(end, :), neg_colormap(1, :), neg_colormap(end, :));



        c_lh_colored = zeros(size(c_lh));
        whpos = c_lh > 0;
        %     kpos = 61;   % which block of 256 colors; depends on colormap
        cpos = map_function(c_lh(whpos), 0, clim(2), (kpos-1)*nvals+1, kpos*nvals); % map into indices in hot cm range of colormap
        c_lh_colored(whpos) = cpos;
        
        whneg = c_lh < 0;
        %     kneg = 55;   % which block of 256 colors
        cneg = map_function(c_lh(whneg), clim(1), 0, (kneg-1)*nvals+1, kneg*nvals); % map into indices in cool cm range of colormap
        c_lh_colored(whneg) = cneg;
    
        wh_lh = c_lh == 0 | isnan(c_lh);                    % save these to replace with gray-scale later
    
        surf = get(o2.surface{4}.axis_handles,'Children');
        c_gray = get(surf(2), 'FaceVertexCData');
        c_gray(wh_lh) = (c_gray(wh_lh) - min(c_gray(wh_lh))) ./ range(c_gray(wh_lh)) .* nvals;




        c_rh_colored = zeros(size(c_rh));
        whpos = c_rh > 0;
        %     kpos = 61;   % which block of 256 colors; depends on colormap
        cpos = map_function(c_rh(whpos), 0, clim(2), (kpos-1)*nvals+1, kpos*nvals); % map into indices in hot cm range of colormap
        c_rh_colored(whpos) = cpos;
        
        whneg = c_rh < 0;
        %     kneg = 55;   % which block of 256 colors
        cneg = map_function(c_rh(whneg), clim(1), 0, (kneg-1)*nvals+1, kneg*nvals); % map into indices in cool cm range of colormap
        c_rh_colored(whneg) = cneg;
                
        wh_rh = c_rh == 0 | isnan(c_rh);                    % save these to replace with gray-scale later
        % Get the original grayscale values
        % -----------------------------------------------------------------------
        % Exclude colored values
        surf = get(o2.surface{1}.axis_handles,'Children');
        c_gray = get(surf(2), 'FaceVertexCData');
        c_gray(wh_rh) = (c_gray(wh_rh) - min(c_gray(wh_rh))) ./ range(c_gray(wh_rh)) .* nvals;
    end
        
    
    c_lh_colored(wh_lh) = c_gray(wh_lh);
    


    c_rh_colored(wh_rh) = c_gray(wh_rh);
    
    %c_colored = map_function(c);    % Map to colormap indices (nvals = starting range, nvals elements)
    
    % FIX: rescale range so we don't map into gray range at edges due to
    % interpolation
    %     xx = c_colored(~wh); % in-blob vertices
    %     xx = 1 + nvals + max(xx) .* (xx - min(xx)) ./ (range(xx) + nvals);
    %     c_colored(~wh) = xx;
    %cm(256, :) = [1 0 1];
    
    
end