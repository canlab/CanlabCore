function tor_polar_plot(vals, colors, names, varargin)
% Make polar line plot(s)
%
% Usage:
% tor_polar_plot(vals, colors, names, ['nofigure'])
%
% vals = cell array, one cell per plot
% in each cell, matrix of observations x variables
% plots one line for each variable.
%
% names is cell array, one cell per plot
% contains cell array of names for each condition
%
% e.g.,
% tor_polar_plot({w+1}, {'r' 'b'}, setnames(1))

dofigure = 1;

p = length(vals); % plots

nr = 1; nc = p;  % rows/cols in subplots
if p > 4, nr = 2; nc = ceil(p./2); end

if any(strcmp(varargin, 'nofigure'))
    % suppress figure
    dofigure = 0;
else
    create_figure('tor_polar', nr, nc)
end

for s = 1:p
    
    if dofigure, subplot(nr, nc, s);, end
    
    set(gca, 'FontSize', 18)
    
    k = size(vals{s}, 2);  % variables to plot (each gets a line)
    
    % make non-negative
    % this changes the interpretation of the origin of the plot
    if any(vals{s}(:) < 0)
        vals{s} = vals{s} - min(vals{s}(:));
    end
    
    for i = 1:k
        if i == 1, hold on, end % must turn off to get polar text automatically. we now do our own.
        
        v = vals{s}(:, i)';
        
        ang = linspace(0, 2*pi, length(v) + 1);
        hh = polar(ang, [v v(1)]);
        
        set(hh, 'LineWidth', 2, 'Color', colors{i})
        hold on
        
    end
    
    % lines and text
    
    % Circles
    maxval = max(max(vals{s}));
    h = circle([0 0], maxval);
    set(h, 'Color', 'k');
    
    h = circle([0 0], maxval ./ 2);
    set(h, 'Color', 'k', 'LineStyle', ':');
    
    % spokes
    [X, Y] = pol2cart(ang, maxval);
    z = zeros(size(X));
    lineh = line([z; X], [z; Y]);
    set(lineh, 'Color', 'k', 'LineStyle', ':');
    
    % names
    
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    
    if ~isempty(names) && ~isempty(names{s})
        
        for j = 1:length(names{s})
            x = X(j);
            if x < 0, x = x - range(xlim)*.02*length(names{s}{j}); end
            
            y = Y(j);
            y = y + sign(y) * range(ylim)* .04;
            
            texth(j) = text(x, y, names{s}{j}, 'FontSize', 14);
        end
    end
    
    axis equal
    axis off
    
    
end % subplot

end % function

