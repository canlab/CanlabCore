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
% Inputs:
% xyz   a 3-vol list of vertices to color
% v     can be a matrix of vertices
%       or a handle to a patch object containing vertices
%       if it's a handle, this function sets the color to interp
%       and the FaceVertexCData to the color matrix c
% actcolor
%       [r g b] activation color
% basecolor
%       [r g b] baseline color - optional.
% mind  optional - min distance to color vertex
%       Vertices within mind of an xyz coordinate will be colored
% colorscale
%       optional.  followed by vector of values by which to multiply input
%       color
%       these are scaled to be between .3 and one.
%       if entered, this will make the colors vary by, for example, Z score
%       so Z-scores are an acceptable input.
%       cscale should be in the same coordinate order as xyz
%       for ADDITIONAL clusters, repeat the 'colorscale', Z argument pair in the function call
%
%       YOU CAN ALSO pass true RGB values for each xyz coordinate in: 'colorscale', rgblist, 
%       IF cscale is a 3-vector, it specifies the ACTUAL colors, and is not scaled to .3 - 1
%
%
% following basecolor and mind:
% additional xyz coordinate lists, with syntax:
% 'vert', xyz2 [your xyz input], [r g b] color for xyz plot
%
% also, you can enter 'ovlcolor' followed by [r g b] for overlaps between xyz sets
%   colors will ONLY appear in the overlap color if they share actual coordinates in common, 
%   not necessarily if surface vertices are within the specified distance from both sets of coords.
%
% to get a good brain surface, try this:
%figure
%p = patch('Faces', faces, 'Vertices', vertices, 'FaceColor', [.5 .5 .5], ...
% 'EdgeColor', 'none', 'SpecularStrength', .2, 'FaceAlpha', 1, 'SpecularExponent', 200);
%lighting gouraud;camlight right
%axis image; myLight = camlight(0, 0);set(myLight, 'Tag', 'myLight');
%set(gcf, 'WindowButtonUpFcn', 'lightFollowView');lightfollowview
%drawnow
%
% ..
% by Tor Wager  August 25, 2002
% ..


function [c, alld] = getVertexColors_old_backup(xyz, v, actcolor, varargin)

    mind = 3;
    basecolor = [.5 .5 .5];
    alld = [];
    xyza = xyz;
    cscale = [];
    allda = [];
    vv = [];

    % -----------------------------------------------------------------------
    % * set up input arguments
    % -----------------------------------------------------------------------
    doalph = 0; cscale = [];
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
        elseif strcmp(varargin{i}, 'alphaone')
            doalph = 1;
        elseif strcmp(varargin{i}, 'allcolor')
            acol = varargin{i+1};
        elseif strcmp(varargin{i}, 'colorscale')
            cscale{end+1} = varargin{i+1};
            if min(size(cscale{end})) == 1
                % scale colors (may be necessary) - only if single vector, not RGB values
                cscale{end} = cscale{end} ./ max(cscale{end});
            end
            if any(cscale{end} < 0), error('Some color scale values are less than zero.'), end

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

    bad = any(size(c) - size(v));
    if bad
        c = repmat(basecolor, size(v, 1), 1);
    end




    % -----------------------------------------------------------------------
    % * main xyz color change
    % -----------------------------------------------------------------------

    t1 = clock;
    fprintf('Main color vertices: ')

    c = change_colors(c, xyz, v, mind, cscale, actcolor, p);
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
            c = change_colors(c, vv{j}, v, mind, cscale(j+1), cc{j}, p);
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
        c = change_colors(c, xyzb, v, mind, cscaletmp, ocol, p);
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
        c = change_colors(c, xyzall, v, mind, cscaletmp, acol, p);
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



function c = change_colors(c, coords, v, mind, cscale, actcolor, p)

    if isempty(coords), disp('Coords is empty. Nothing to plot.'), return, end

    % select vertices that are even close
    cmax = max(coords, [], 1);
    cmin = min(coords, [], 1);
    fprintf('%3.0f vertices.  selecting: ', size(v, 1));
    wh = any(v - repmat(cmax, size(v, 1), 1) > mind, 2);
    wh2 = any(repmat(cmin, size(v, 1), 1) - v > mind, 2);

    % list vertices to test and possibly change color
    whverts = (1:size(v, 1))';
    whverts(wh | wh2) = [];     % vertex indices in big list
    smallv = v;
    smallv(wh | wh2,:) = [];     % vertices--restricted list

    fprintf('%3.0f\n', size(whverts, 1));

    if isempty(smallv), return, end

    % select coords that are even close
    % ----------------------------------
    cmax = max(smallv, [], 1);
    cmin = min(smallv, [], 1);
    fprintf('%3.0f coords.  selecting: ', size(coords, 1));
    wh = any(coords - repmat(cmax, size(coords, 1), 1) > mind, 2);
    wh2 = any(repmat(cmin, size(coords, 1), 1) - coords > mind, 2);
    
    
    % if cscale is matrix, must select these values of cscale as well!
    if size(cscale{1}, 1) == size(coords, 1)
        cscale{1}(wh | wh2,:) = [];
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
            % find vertices that are within range of point i in set  setno
            vertex_indices = find_in_radius(xyz2, setno, i, smallv, mind, whverts);

            % two modes: if cscale{1} is a matrix, treats as rgb values, and put in
            % color stored in cscale.  if cscale{1} is a vector, treat it as a scaling value for actcolor
            % In either case, indxval should index location of coordinate in FULL
            % list (corresponding to full list in cscale)
            
            c = color_change_vertices(c, mind, indxval, cscale, actcolor, vertex_indices);

            if ~isempty(vertex_indices), wh_coords_near_surface(indxval) = 1; end
            
            indxval = indxval + 1;
        end

        if exist('p', 'var')
            set(p, 'FaceColor', 'interp')
            set(p, 'FaceVertexCData', c)
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




function c = color_change_vertices(c, mind, i, cscale, actcolor, vertex_indices)

    % two modes: if cscale{1} is a matrix, treats as rgb values, and put in
    % color stored in cscale
    % if cscale{1} is a vector, treat it as a scaling value for actcolor
    % In either case, the ith index point should correspond to the
    % coordinate being worked on.

    n = length(vertex_indices);

    if n
        if min(size(cscale{1})) > 1  % if we have rgb values rather than scaling values
            % this occurs if heatmap = yes and colorscale = no, we pass in rgb values
            c(vertex_indices,:) = repmat(cscale{1}(i,:), n, 1);
        else
            c(vertex_indices,:) = repmat(actcolor.*cscale{1}(i,:), n, 1);
        end
    end

end

