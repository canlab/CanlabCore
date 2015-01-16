function han = plot_horizontal_line(y,color)
    % function han = plot_horizontal_line(y,color)
    
    if nargin < 2, color = 'k'; end
    xlim = get(gca,'XLim');
    han = plot(xlim, [y y], color);

    return
    