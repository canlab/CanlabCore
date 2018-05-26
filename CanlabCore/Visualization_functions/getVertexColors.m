function [c, alld] = getVertexColors(xyz, v, actcolor, varargin)
% :Usage:
% ::
%
%    [c, alld] = getVertexColors(xyz, v, actcolor, [basecolor], [mind], 'vert', [xyz2], [actcolor2], 'vert', [xyz3], [actcolor3])
%
% given a point list of XYZ mm coordinates (3 columns)
% and a list of vertices in an isosurface, 
% returns FaceVertexCData color values for brain near points and brain not near points.
% c is vertex color specification, 3 columns indicating RGB values
%
% :Inputs:
%
%   **xyz:**
%        a 3-vol list of vertices to color
%
%   **v:**
%        can be a matrix of vertices
%        or a handle to a patch object containing vertices
%        if it's a handle, this function sets the color to interp
%        and the FaceVertexCData to the color matrix c
%
%   **actcolor:**
%        [r g b] activation color
%
%   **basecolor:**
%        [r g b] baseline color - optional.
%
%   **mind:**
%        optional - min distance to color vertex
%        Vertices within mind of an xyz coordinate will be colored
%
%   **colorscale
%       optional.  followed by vector of values by which to multiply input
%       color
%       these are scaled to be between .3 and one.
%
%       if entered, this will make the colors vary by, for example, Z score
%       so Z-scores are an acceptable input.
%
%       cscale should be in the same coordinate order as xyz
%
%       for ADDITIONAL clusters, repeat the 'colorscale', Z argument pair in the function call
%
%       YOU CAN ALSO pass true RGB values for each xyz coordinate in: 'colorscale', rgblist, 
%       IF cscale is a 3-vector, it specifies the ACTUAL colors, and is not scaled to .3 - 1
%
%   Following basecolor and mind:
%
%   additional xyz coordinate lists, with syntax:
%   vert', xyz2 [your xyz input], [r g b] color for xyz plot
%
%   also, you can enter 'ovlcolor' followed by [r g b] for overlaps between xyz sets
%   colors will ONLY appear in the overlap color if they share actual coordinates in common, 
%   not necessarily if surface vertices are within the specified distance from both sets of coords.
%
% :Examples:
% ::
%
%    % to get a good brain surface, try this:
%    figure
%    p = patch('Faces', faces, 'Vertices', vertices, 'FaceColor', [.5 .5 .5], ...
%            å  'EdgeColor', 'none', 'SpecularStrength', .2, 'FaceAlpha', 1, 'SpecularExponent', 200);
%    lighting gouraud;
%    camlight right
%    axis image;
%    myLight = camlight(0, 0);
%    set(myLight, 'Tag', 'myLight');
%    set(gcf, 'WindowButtonUpFcn', 'lightFollowView');
%    lightfollowview
%    drawnow
%
% ..
% by Tor Wager  August 25, 2002
% ..


    global do_fsavg_right do_fsavg_left ras
    
    mind = 3;
    basecolor = [.5 .5 .5];
    alld = [];
    xyza = xyz;
    cscale = [];
    allda = [];
    vv = [];
    doalph = 0; 
    cscale = [];
    alphascale = {};
    do_fsavg_left = false;
    do_fsavg_right = false;
    
    % -----------------------------------------------------------------------
    % * set up input arguments
    % -----------------------------------------------------------------------
    
    if ~isempty(varargin), basecolor = varargin{1}; end
    
    if length(varargin) > 1, mind = varargin{2}; end
    
    ind = 1;
    for i = 3:length(varargin)
        if strcmp(varargin{i}, 'vert')
            vv{ind} = varargin{i+1};

            % intersections
            xyzb{ind, 1} = intersect(xyz, vv{ind}, 'rows');
            for j = ind-1:-1:1
                xyzb{ind, j} = intersect(vv{j}, vv{ind}, 'rows');
                xyza = intersect(xyza, xyzb{ind, j}, 'rows');
            end

            cc{ind} = varargin{i+2};
            ind = ind+1;
        elseif strcmp(varargin{i}, 'ovlcolor')
            ocol = varargin{i+1};
            
        elseif strcmp(varargin{i}, 'alphascale')
            doalph = 1;
            alphascale{1} = varargin{i + 1};
            
        elseif strcmp(varargin{i}, 'allcolor')
            acol = varargin{i+1};
            
        elseif strcmp(varargin{i}, 'colorscale')
            
            cscale{end+1} = varargin{i+1};
            
            % now do this ahead of time
            %             if min(size(cscale{end})) == 1
            %                 % scale colors (may be necessary) - only if single vector, not RGB values
            %                 cscale{end} = cscale{end} ./ max(cscale{end});
            %             end
            if any(cscale{end} < 0), error('Some color scale values are less than zero.'), end
        elseif strcmp(varargin{i}, 'fsavg_left')
            ras = load(which('lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat'));
            do_fsavg_left = true;
        elseif strcmp(varargin{i}, 'fsavg_right')
            % uses freesurfer inflated brain with Thomas Yeo group's RF_ANTs mapping
            % from MNI to Freesurfer. (https://doi.org/10.1002/hbm.24213)
            ras = load(which('rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat'));
            do_fsavg_right = true;
        end
    end

    if isempty(cscale)
        cscale{1} = ones(size(xyz, 1), 1);
        for i = 1:length(vv)    % additional vertices
            cscale{end+1} = ones(size(vv{i}, 1), 1);
        end
    end

    %if ~isempty(cscale) don't need this.
    %    if length(cscale) < length(vv)
    %        cscale{length(vv)} = [];
    %    end
    %end

    
    if ishandle(v)
        p = v;
        v = get(p, 'Vertices');
        c = get(p, 'FaceVertexCData');   % get existing colors from surface
    else
        % uh-oh, not a handle!
        warning('Figure handle missing: Figure was closed?')
        p = findobj('Type', 'patch');
        v = get(p(1), 'Vertices');
        c = get(p, 'FaceVertexCData');   % get existing colors from surface
    end
    
    
    % if FaceVertexCData is v x 1 matrix, it is mapped to the figure colormap
    % if it is v x 3, it is true color.
    % We can replace the colormapped FaceVertexCData with its true-color
    % equivalent, making it colormap independent.
    % This will allow us to use the existing values and average in new colors
    % for blobs
    if  ishandle(p) && size(v,1) == size(c, 1) && size(c, 2) == 1
        cm = get(gcf, 'Colormap');
        c2 = round((c./max(c)) .* length(cm));
        c2(c2==0)=1;
        c2 = cm(c2, :);
        c = c2;  % this is the true-color equivalent
        clear cm c2
        set(p, 'FaceVertexCData', c);
    end
    
    bad = any(size(c) - size(v)); % this could occur for various reasons?
    
    if bad
        c = repmat(basecolor, size(v, 1), 1);
    end
    
% c = true-color vertices x 3
% v = xyz coords of vertices x 3
% cscale: mapped true colors or index for  xyz, cell with cscale{1} = xyz coords to color x 3 or x 1

    % -----------------------------------------------------------------------
    % * main xyz color change
    % -----------------------------------------------------------------------

    t1 = clock;
    fprintf('Main color vertices: ')

    c = change_colors(c, xyz, v, mind, cscale, actcolor, p, alphascale);
    drawnow();



    % -----------------------------------------------------------------------
    % * additional optional vertices
    % -----------------------------------------------------------------------
    if exist('vv', 'var')
        for j = 1:length(vv)
            fprintf('\nAdditional vertices: ')

            % cscale:
            % pass in vector of ones length coords for solid color
            % or scalar vals for color mapping
            % pass in cell array

            % cc{j} is color, cscale{j+1} is scaling vals for coords
            c = change_colors(c, vv{j}, v, mind, cscale(j+1), cc{j}, p, alphascale);
        end
        drawnow();
    end


    % -----------------------------------------------------------------------
    % * figure out which vertices should be colored with ocol (overlap color)
    % -----------------------------------------------------------------------

    if exist('ocol', 'var') && exist('xyzb', 'var') && ~isempty(cat(1, xyzb{:}))

        alld = zeros(size(v, 1), 1);  % keeps track of overlap vertices
        xyzb = cat(1, xyzb{:});

        t1 = clock;
        fprintf('\nOverlap vertices: ')

        cscaletmp = {ones(size(xyzb, 1), 1)};
        c = change_colors(c, xyzb, v, mind, cscaletmp, ocol, p, alphascale);
        drawnow();

        fprintf('%3.0f.done in %3.0f s\n', i, etime(clock, t1))
    end


    % -----------------------------------------------------------------------
    % * figure out which vertices should be colored with acol (all color)
    % -----------------------------------------------------------------------
    if exist('acol', 'var') && ~isempty(xyza) && length(vv)>1

        xyzall = xyza;
        for i = 2:length(vv)
            xyzall = intersect(xyzall, vv{i}, 'rows');
        end

        t1 = clock;
        fprintf('\nAll overlap vertices: ')

        cscaletmp = {ones(size(xyzall, 1), 1)};
        c = change_colors(c, xyzall, v, mind, cscaletmp, acol, p, alphascale);
        drawnow();

        fprintf('%3.0f.done in %3.0f s\n', i, etime(clock, t1))

    end


    % -----------------------------------------------------------------------
    % * final color change
    % -----------------------------------------------------------------------


    if exist('p', 'var')
        set(p, 'FaceColor', 'interp')
        set(p, 'FaceVertexCData', c)
        drawnow()
    end

    %lightFollowView
end






% -----------------------------------------------------------------------
% * SUB-FUNCTIONS
% -----------------------------------------------------------------------



function c = change_colors(c, coords, v, mind, cscale, actcolor, p, alphascale)

% c is list of colors. by default, usually [.5 .5 .5] for all colors
% (basecolor)
    global do_fsavg_left do_fsavg_right ras

    if isempty(coords), disp('Coords is empty. Nothing to plot.'), return, end

    % select surface vertices that are close to coords
    % vertices that are far will be omitted
    % smallv is a reduced set of vertices
    % --------------------------------------------------------------------
    cmax = max(coords, [], 1);
    cmin = min(coords, [], 1);
    fprintf('%3.0f vertices.  selecting: ', size(v, 1));
    
    wh = any(v - repmat(cmax, size(v, 1), 1) > mind, 2);
    wh2 = any(repmat(cmin, size(v, 1), 1) - v > mind, 2);

    % list vertices to test and possibly change color
    whverts = (1:size(v, 1))';
    whverts(wh | wh2) = [];     % indices of smallv in big (original) list
    if do_fsavg_left || do_fsavg_right
        smallv = ras.ras';
    else
        smallv = v;
    end
    smallv(wh | wh2,:) = [];     % vertices--restricted list

    fprintf('%3.0f\n', size(whverts, 1));

    if isempty(smallv), return, end

    % select coords that are close to surface vertices5
    % coords, cscale{1}, and alphascale{1} all have same indices
    % ---------------------------------------------------------------------
    cmax = max(smallv, [], 1);
    cmin = min(smallv, [], 1);
    fprintf('%3.0f coords.  selecting: ', size(coords, 1));
    
    % omit wh and wh2 - outside scope of this coord set
    wh = any(coords - repmat(cmax, size(coords, 1), 1) > mind, 2);
    wh2 = any(repmat(cmin, size(coords, 1), 1) - coords > mind, 2);
    
    % if cscale is matrix, must select these values of cscale as well!
    if length(cscale) > 0 && size(cscale{1}, 1) == size(coords, 1)
        cscale{1}(wh | wh2,:) = [];
    end
    
    if length(alphascale) > 0 && size(alphascale{1}, 1) == size(coords, 1)
        alphascale{1}(wh | wh2,:) = [];
    end
    
    coords(wh | wh2,:) = [];
    
    if isempty(coords), return, end

    nc = size(coords, 1);
    fprintf('%3.0f\n', nc);

    % break up coords into list and run
    % ----------------------------------
    
    % break up coords into list
    xyz2 = {}; indx = 1;
    for kk = 1:1000:nc
        setwh{indx} = (kk:min(nc, kk + 1000 - 1))';
        xyz2{indx} = coords(setwh{indx},:);

        indx = indx + 1;
    end

    fprintf('Running %3.0f sets of coordinates: 000', length(xyz2));

    indxval = 1;
    wh_coords_near_surface = false(size(coords, 1), 1);
    
    for setno = 1:length(xyz2)
        fprintf('\b\b\b%03d', setno);

        for i = 1:size(xyz2{setno}, 1)
            % i indexes coordinate; e.g., 1000 iterations for a typical
            % coord set
            % find vertices that are within range of point i in set  setno
            % vertex_indices are indices into the ORIGINAL list 
            vertex_indices = find_in_radius(xyz2, setno, i, smallv, mind, whverts);

            % two modes: 
            % - if cscale{1} is a matrix, treats as rgb values, and put in
            % color stored in cscale.  actcolor is ignored.
            % - if cscale{1} is a vector, treat it as a scaling value for actcolor
            % In either case, indxval should index location of coordinate in FULL
            % list (corresponding to full list in cscale)
            
            % cscale{1}(i) and alphascale{1}(i) contain unique values for each coord i
            % these are used to map colors for single point i to multiple
            % nearby vertices vertex_indices
            
            if length(cscale) > 0 && ~isempty(cscale{1})
                mycscale = cscale{1}(setwh{setno}(i), :);  % color for this point, index into full list of coords/scalevals
            else
                mycscale = [];
            end
            
            if length(alphascale) > 0 && ~isempty(alphascale{1})
                myalphascale = alphascale{1}(setwh{setno}(i), :);  % alpha for this point, index into full list of coords/scalevals
            
                % adjust by cube of mind, to adjust for overlap in vertices
                % affected by adjacent coords - done now in cluster_surf
                %myalphascale = myalphascale ./ mind^3;
                
            else
                myalphascale = [];
            end
            
            c = color_change_vertices(c, mycscale, actcolor, vertex_indices, myalphascale);

            if ~isempty(vertex_indices), wh_coords_near_surface(indxval) = 1; end
            
            indxval = indxval + 1;
        end

        if exist('p', 'var')
            set(p, 'FaceColor', 'interp')
            set(p, 'FaceVertexCData', c)
            drawnow
        end
    end
end






function z = dist_tmp(w, p)

    %
    % if isstr(w)
    %   switch (w)
    %     case 'deriv', 
    %       z = '';
    %     otherwise
    %       error('Unrecognized code.')
    %   end
    %   return
    % end

    % CALCULATION
    if nargin == 1
        p = w;
        w = w';
    end

    [S, R] = size(w);
    [R2, Q] = size(p);
    if (R ~= R2), error('Inner matrix dimensions do not match.'), end

    z = zeros(S, Q);
    if (Q<S)
        p = p';
        copies = zeros(1, S);
        for q=1:Q
            z(:,q) = sum((w-p(q+copies,:)).^2, 2);
        end
    else
        w = w';
        copies = zeros(1, Q);
        for i=1:S
            z(i,:) = sum((w(:,i+copies)-p).^2, 1);
        end
    end
    z = sqrt(z);

end






function [vertex_indices, d] = find_in_radius(xyz2, setno, i, smallv, mind, whverts)
    % output: indices of vertices in BIG list, and distances
    
    % get vertices v within box -- fast method
    wh = find(all(abs(bsxfun(@minus, xyz2{setno}(i,:), smallv)) <= mind, 2));
    d = dist_tmp(smallv(wh,:), xyz2{setno}(i,:)');

    % convert back to big list
    vertex_indices = whverts(wh(d < mind));
end




function c = color_change_vertices(c, mycscale, actcolor, vertex_indices, myalphascale)

% pass in/out larger list c, but update only a few coordinates
% indexed by vertex_indices with a particular color value.  

% cscale is ordered from least color
% change to most for graphical/aesthetic purposes, so this function is
% called repetitively to index lower i (less change) to higher i (more
% change) successively.

% two modes: if cscale{1} is a matrix, treats as rgb values, and put in
% color stored in cscale
% if cscale{1} is a vector, treat it as a scaling value for actcolor
% In either case, the ith index point should correspond to the
% coordinate being worked on.

n = length(vertex_indices);

if n
    if length(mycscale) > 1
        % if we have rgb values rather than scaling values
        % this occurs if heatmap = yes and colorscale = no, we pass in rgb values
        % actcolor is ignored
        
        if isempty(myalphascale)
            % solid colors
            c(vertex_indices,:) = repmat(mycscale, n, 1);
        else
            w = myalphascale;  % the weight for new color vs. old.
            c(vertex_indices, :) = (1-w) .* c(vertex_indices, :) + w .* repmat(mycscale, n, 1);
        end
        
    else
        % this will scale the activation color in proportion to cscale
        % for each voxel
        
        w = mycscale;  % the weight for new color vs. old.  1 = all new, 0 = all old
        c(vertex_indices, :) = (1-w) .* c(vertex_indices, :) + w .* repmat(actcolor, n, 1);
    end
end

end

