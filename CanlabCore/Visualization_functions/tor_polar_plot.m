function [hh, hhfill, ang] = tor_polar_plot(vals, colors, names, varargin)
% Make polar line plot(s)
%
% :Usage:
% ::
%
%    hh = tor_polar_plot(vals, colors, names, ['nofigure'])
%
% :Inputs:
%
%   **vals:**
%        cell array, one cell per plot
%
%        in each cell, matrix of observations x variables
%        plots one line for each variable.
%
%   **names:**
%        is cell array, one cell per plot
%
%        contains cell array of names for each condition
%
% :Optional Inputs:
%
%   **'nofigure':**
%        suppress figure
%
%   **'nonneg':**
%        make all values non-negative by subtracting min value from all values in series (plot)
%
%   **'nofill':**
%        Do not fill in polygons
%
%   **'nonumbers':**
%        Suppress numbers (values on polar axis)
%
%   **'linesonly':**
%        Suppress all plots but lines; useful for adding to existing plots
%
%   **'fixedrange':**
%        Set min and max of circles numbers (values on polar axis)
%        Follow by range vector: [min_val max_val]
% :Output:
%
%   **hh:**
%        Handles to line objects
%
%   **hhfill:**
%        Handles to fill objects
%
%   **ang:**
%        Angles for polar plot.  
%        If v is a data vector, plot new lines with polar(ang, [v v(1)]);
%
% :Examples:
% ::
%
%    tor_polar_plot({w+1}, {'r' 'b'}, setnames(1))
%
% Note: Dark grey inner line is zero point if nonneg option is used.
% Otherwise, zero is the origin.


dofigure = 0;
dofill = 1;
dononneg = 0;
donumbers = 1;
dofixrange = 0;
dolinesonly = false;

p = length(vals); % plots

nr = 1; nc = p;  % rows/cols in subplots
if p > 4, nr = 2; nc = ceil(p./2); end

hh = {};
hhfill = {};

if any(strcmp(varargin, 'nofigure'))
    % suppress figure
    dofigure = 0;
end

if any(strcmp(varargin, 'dofigure'))
    % include figure [now default]
    dofigure = 1;
end

if dofigure
    create_figure('tor_polar', nr, nc)
end

if any(strcmp(varargin, 'nonneg'))
    dononneg = 1;
end

if any(strcmp(varargin, 'nofill'))
    dofill = 0;
end

if any(strcmp(varargin, 'fixedrange'))
    dofixrange = 1;
    fixedrange=varargin{find(strcmp(varargin, 'fixedrange'))+1};
end

if any(strcmp(varargin, 'nonumbers'))
    donumbers = 0;
end

if any(strcmp(varargin, 'linesonly'))
    dolinesonly = true;
end

for s = 1:p % s indexes plot series
    
    if dofigure, subplot(nr, nc, s); end
    
    set(gca, 'FontSize', 18)
    
    [k, numvars] = size(vals{s});  % obs points and variables to plot (each gets a line)
    
    % make non-negative
    % this changes the interpretation of the origin of the plot
    if dononneg && ~dofixrange
        if any(vals{s}(:) < 0)
            origmin(s) = abs(min(vals{s}(:)));  % this is the new zero point
            vals{s} = vals{s} - min(vals{s}(:));
        else
            origmin(s) = 0;
        end
    elseif dofixrange
        vals{s} = vals{s} - fixedrange(1);
        origmin(s) = abs(fixedrange(1));
    end
    
    % Polar plots of values
    % ------------------------------------------------------------
    
    for i = 1:numvars
        
        if i == 1, hold on, end % must turn off to get polar text automatically. we now do our own.
        
        v = vals{s}(:, i)';
        ang = linspace(0, 2*pi, k + 1);
        
        hh{s}(i) = polar(ang, [v v(1)]);
        
        % extend colors by replicating if needed, i.e., if we want them all the same color and only entered one 
        if i > length(colors), colors{end + 1} = colors{end}; end
        
        set(hh{s}(i), 'LineWidth', 2, 'Color', colors{i})
        hold on
        
        if dofill
            [Xf, Yf] = pol2cart(ang, [v v(1)]);
            hhfill{s}(i) = fill(Xf,Yf, colors{i});
            set(hhfill{s}(i), 'FaceAlpha', .1, 'EdgeColor', colors{i});
        end
        
    end
    
    if dolinesonly, return, end
    
    % lines and text
    % ------------------------------------------------------------
    
    % Circles
    if dofixrange
        maxval=fixedrange(2)-fixedrange(1);
    else
        maxval = max(max(vals{s})) + .1 * (max(max(vals{s})));
    end
    
    h = circle([0 0], maxval);
    set(h, 'Color', [.3 .3 .3], 'LineWidth', 2);
    
    if dononneg || dofixrange
        minval = origmin(s);
    else
        minval = 0;
    end
    
    circlevals = linspace(minval, maxval, 4); % 4 is no. or circles
    
    for i = 1:length(circlevals)
        h = circle([0 0], circlevals(i));
        if i == 1
            set(h, 'Color', [.3 .3 .3], 'LineWidth', 2);
        else
            
            set(h, 'Color', [.5 .5 .5], 'LineStyle', ':');
        end
    end
    
    % Spokes, and set angles
    % ------------------------------------------------------------
    
    %for i = 1:k
        
        ang = linspace(0, 2*pi, k + 1);
        
        [X, Y] = pol2cart(ang, maxval);
        z = zeros(size(X));
        lineh = line([z; X], [z; Y]);
        set(lineh, 'Color', 'k', 'LineStyle', ':');
        
    %end
    
    % Axis label names
    % ------------------------------------------------------------
    
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    
    myang = ang(1:end-1) * 360 / (2*pi); % Radians to degrees
    
    % adjust angles so text on left will not be upside down
    wh = myang > 90 & myang < 270;
    myang(wh) = myang(wh) - 180;
    
    clear myextent

    % set font size now - harmonize with image_similarity_plot - so extent
    % calc will be correct
    mytextsize = 30 ./ (k.^.3);
                
    if ~isempty(names) && ~isempty(names{s})
        
        if k ~= length(names{s})
            error('Length of names cell does not match data. Check inputs!'); 
        end
            
        for j = 1:length(names{s})
            x = X(j);
            %if x < 0, x = x - range(xlim)*.01*length(names{s}{j}); end
            
            y = Y(j);
            %y = y + sign(y) * range(ylim)* .04;
            
            texth(j) = text(x, y, names{s}{j}, 'FontSize', mytextsize);
            
        end
        
        % adjust axis to make room for text labels
        % Must be done before text is shifted to account for extent
        axis equal
        
        for j = 1:length(names{s})
            myextent_j = get(texth(j), 'Extent');
            myextent(1, j) = myextent_j(3) + .02 * maxval; %.1 * myextent_j(3); %maxval;
        end
        mylim = get(gca, 'XLim');
        mylim = mylim + [-max(myextent) max(myextent)];
        set(gca, 'XLim', mylim);
        
    % Final calculation of extent and rotation setting
        for j = 1:length(names{s})
            % For text shifting, get extent before rotating.
            % axis equal/scaling and size matters when getting extent, too, so make sure correct
            % size ahead of doing this!!
            myextent_j = get(texth(j), 'Extent');
            myextent(1, j) = myextent_j(3) + .02 * maxval; %.1 * myextent_j(3); %maxval;
            
            set(texth(j), 'Rotation', myang(j))
            
        end
    end
    
    % adjust max position so text will end at outer circle
    % for inverted text (needs to be shifted)
    textstart = repmat(maxval, 1, k);
    textstart(wh) = textstart(wh) + myextent(wh);
    textstart(~wh) = textstart(~wh) + .02 * maxval;
    
    % Re-calculate X and Y position accounting for shifted text
    for j = 1:k
        
        [X(j), Y(j)] = pol2cart(ang(j), textstart(j));
        
        mypos = get(texth(j), 'Position');
        
        mypos(1:2) = [X(j) Y(j)];
        
        set(texth(j), 'Position', mypos);
        
    end
    
    
    
    % Value labels (text) for polar guide lines
    % ------------------------------------------------------------
    if donumbers
        
        len = length(circlevals);
        [Xc, Yc] = pol2cart(linspace(0, -.7, len), circlevals);
        for i = 1:len
            text(Xc(i), Yc(i), sprintf('%3.2f', circlevals(i) - minval), 'FontSize', 14, 'Color', [.4 .4 .4]);
        end
        
    end
    
    %axis equal
    axis off
    
    
end % subplot

end % function

