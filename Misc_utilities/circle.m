function [h, fillh] = circle(center, radius, varargin)
    % draws a circle
    %[h, fillh] = circle(center, radius, ['fill'], [color string or numbers])
    %
    % 'fill' requires the geom2d toolbox, by David Legland


    x = pi*[0:.5:2];
    y = [0  1  0 -1  0  1  0;
        1  0  1  0 -1  0  1] * radius;
    pp = spline(x, y);
    circ_data = ppval(pp, linspace(0, 2*pi, 101));
    h = plot(circ_data(1,:) + center(1), circ_data(2,:) + center(2));
    %axis equal

    if length(varargin) > 0 && strcmp(varargin{1}, 'fill')

        xy = circleAsPolygon([center radius], 20);
        fillh = fillPolygon(xy);

        delete(h)
        h = [];
        
        if length(varargin) > 1
            set(fillh, 'FaceColor', varargin{2})
        end

    end

end