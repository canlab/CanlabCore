function [h, fillh] = circle(center, radius, varargin)
% draws a circle at x,y center = [center(1) center(2)]
%
% :Usage:
% ::
%
%     [h, fillh] = circle(center, radius, ['fill'], [color string or numbers])
%
% Examples:
%
% figure; hold on;
% circle([0 0], 1);
% circle([1 1], .5, 'fill');
% circle([-.5 .3], .7, 'fill', [1 1 0]);
%
% Tor Wager. Latest update: 12/2018

dofill = false;
fillalpha = .5;
linecolor = [.2 .2 .2];
fillcolor = [.3 .3 1];

if length(varargin) > 0 && strcmp(varargin{1}, 'fill')
    dofill = true;
end

if length(varargin) > 1
    fillcolor = varargin{2};
end
    
[xunit, yunit] = circle_coords(center, radius);

h = plot(xunit, yunit, 'Color', linecolor);


    if dofill

        fillh = fill(xunit, yunit, fillcolor, 'FaceAlpha', fillalpha);

        delete(h)
        h = [];
        
%         if length(varargin) > 1
%             set(fillh, 'FaceColor', varargin{2})
%         end

    end

end

% Old code:
%     x = pi*[0:.5:2];
%     y = [0  1  0 -1  0  1  0;
%         1  0  1  0 -1  0  1] * radius;
%     pp = spline(x, y);
%     circ_data = ppval(pp, linspace(0, 2*pi, 101));
%     h = plot(circ_data(1,:) + center(1), circ_data(2,:) + center(2));
%     %axis equal
%         xy = circleAsPolygon([center radius], 20);
%         fillh = fillPolygon(xy);


function [xunit, yunit] = circle_coords(center, r)

th = 0:pi/50:2*pi;

xunit = r * cos(th) + center(1);
yunit = r * sin(th) + center(2);

end

