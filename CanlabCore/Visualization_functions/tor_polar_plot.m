function [hh, hhfill] = tor_polar_plot(vals, colors, names, varargin)
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
% :Output:
%
%   **hh:**
%        Handles to line objects
%
%   **hhfill:**
%        Handles to fill objecta
%
% :Examples:
% ::
%
%    tor_polar_plot({w+1}, {'r' 'b'}, setnames(1))
%
% Note: Dark grey inner line is zero point if nonneg option is used.
% Otherwise, zero is the origin.


dofigure = 1;
dofill = 1;
dononneg = 0;
donumbers = 1;

p = length(vals); % plots

nr = 1; nc = p;  % rows/cols in subplots
if p > 4, nr = 2; nc = ceil(p./2); end

hh = {};
hhfill = {};

if any(strcmp(varargin, 'nofigure'))
    % suppress figure
    dofigure = 0;
else
    create_figure('tor_polar', nr, nc)
end

if any(strcmp(varargin, 'nonneg'))
    dononneg = 1;
end

if any(strcmp(varargin, 'nofill'))
    dofill = 0;
end

if any(strcmp(varargin, 'nonumbers'))
    donumbers = 0;
end

for s = 1:p
    
    if dofigure, subplot(nr, nc, s); end
    
    set(gca, 'FontSize', 18)
    
    [k, numvars] = size(vals{s});  % obs points and variables to plot (each gets a line)
    
    % make non-negative
    % this changes the interpretation of the origin of the plot
    if dononneg
        if any(vals{s}(:) < 0)
            origmin(s) = abs(min(vals{s}(:)));  % this is the new zero point
            vals{s} = vals{s} - min(vals{s}(:));
        else
            origmin(s) = 0;
        end
    end
    
    % Polar plots of values
    % ------------------------------------------------------------
    
    for i = 1:numvars
        if i == 1, hold on, end % must turn off to get polar text automatically. we now do our own.
        
        v = vals{s}(:, i)';
        ang = linspace(0, 2*pi, k + 1);
        
        hh{s}(i) = polar(ang, [v v(1)]);
        
        set(hh{s}(i), 'LineWidth', 2, 'Color', colors{i})
        hold on
        
        if dofill
            [Xf, Yf] = pol2cart(ang, [v v(1)]);
            hhfill{s}(i) = fill(Xf,Yf, colors{i});
            set(hhfill{s}(i), 'FaceAlpha', .1, 'EdgeColor', colors{i});
        end
        
    end
    
    % lines and text
    % ------------------------------------------------------------
    
    % Circles
    maxval = max(max(vals{s})) + .1 * (max(max(vals{s})));
    h = circle([0 0], maxval);
    set(h, 'Color', [.3 .3 .3], 'LineWidth', 2);
    
    if dononneg
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
    
    for i = 1:k
        
        ang = linspace(0, 2*pi, k + 1);
        
        [X, Y] = pol2cart(ang, maxval);
        z = zeros(size(X));
        lineh = line([z; X], [z; Y]);
        set(lineh, 'Color', 'k', 'LineStyle', ':');
        
    end
    
    % Axis label names
    % ------------------------------------------------------------
    
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    
    if ~isempty(names) && ~isempty(names{s})
        
        for j = 1:length(names{s})
            x = X(j);
            %if x < 0, x = x - range(xlim)*.01*length(names{s}{j}); end
            
            y = Y(j);
            %y = y + sign(y) * range(ylim)* .04;
            
            texth(j) = text(x, y, names{s}{j}, 'FontSize', 14);
            
            myang = ang(j) * 360 / (2*pi);
            if myang > 90 && myang < 270
                myang = myang - 180;
                
                [newx, newy] = pol2cart(ang(j), maxval + .023*length(names{s}{j}));
                set(texth(j), 'Position', [newx newy]);
                
            end
            set(texth(j), 'Rotation', myang)
            
        end
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
    
    axis equal
    axis off
    
    
end % subplot

end % function

