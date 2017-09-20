function out = riverplot_line(x, y, color, thickness, steepness, varargin)
% out = riverplot_line(x, y, color, thickness, steepness, varargin)
%
% Draw sigmoidal line from 2-coord point x to 2-coord point y
% x = [xcoord ycoord] for left-hand point
% y = [xcoord ycoord] for right-hand point
%
% x = rect1.topright
% y = rect2.topleft

% Programmers' notes:
% 8/21/2017 Stephan Geuter
% changed steepnes treatment computation. Steepness coefficient is now used
% in sigmoid function to control steepnes and curvature at tangets (box
% edges). shift in xstartpoints is removed. Makes ribbons smoother.
%
%


% x, y inputs can be [x,y,z] coordinate triplets or clusters.
if isstruct(x), x = x.mm_center(1:2); end
if isstruct(y), y = y.mm_center(1:2); end

% make middle points with sigmoid function

xdiff = y(1) - x(1);  % difference in x position
ydiff = y(2) - x(2);

% xstartpoints = [x(1) + xdiff * steepness(1) y(1) - xdiff * steepness(1)]';
xstartpoints = [x(1) y(1)]';


% sigmoid reference curve
% sigmoid = inline('p(1) .* ( 1 ./ (1 + p(2)*exp(-p(3)*x)) )','p','x');
sigmoid = @(p,x) (p(1) .* ( 1 ./ (1 + p(2)*exp(-p(3)*x)) ));
xx = linspace(-5, 5, 50); % generate 50 points with standard sigmoid  
%                         % Range determines steepness bounds, larger range
%                         = more steep
yy = sigmoid([1 1 1+steepness(1)], xx)'; % add steepness coefficient to have ribbons meet the box edges


xmidpoints = linspace(xstartpoints(1), xstartpoints(end), 50)';
ymidpoints = x(2) + yy * ydiff;

xcoords = [x(1); xstartpoints(1); xmidpoints; xstartpoints(end); y(1)];
ycoords = [x(2); x(2); ymidpoints; y(2); y(2)];

%     % make x, y bend percents
%     if length(steepness) == 1
%         steepness = repmat(steepness, 1, 2);
%     end
% 
%     % make 4 coords, adding 2 middle points, so we can bend
%     xdiff = y(1) - x(1);  % difference in x position
%     xmidpoints = [x(1) + xdiff * steepness(1) x(1) + xdiff * .5  y(1) - xdiff * steepness(1)]';
%     
%     xcoords = [x(1); xmidpoints; y(1)];
% 
%     ydiff = (y(2) - x(2)) ./ 2;
%     ycoords = [x(2); x(2); x(2) + ydiff; y(2); y(2)];

%     plot(xcoords, ycoords, 'go');
%     
%      if any(steepness)
% 
%     %[p,sse,fit] = nonlin_fit(ycoords, xcoords - mean(xcoords), 'start',[1 1 1]);
%     
% sigmoid = inline('p(1) .* ( 1 ./ (1 + p(2)*exp(-p(3)*x)) )','p','x');
% xcoords2 = linspace(xcoords(1), xcoords(end), 100);  
% 
% ycoords2 = sigmoid([1 1 1], xcoords2 - mean(xcoords2));
% ycoords2 = ycoords(1) + ydiff * ycoords2;
% 
% plot(xcoords2, ycoords2, 'rx');

%         % bow out: curved line
%         n = length(xcoords);
% 
%         nsamples = [];
%         if length(varargin) > 0
%             nsamples = varargin{1};
%         end
%         
%         if isempty(nsamples)
%             nsamples = 10 * n;
%         end
% 
%         t = 1:n;
%         ts = 1:((n-1)/(nsamples-1)):n;          % spline grid
% 
%         xcoords = spline(t, xcoords, ts);
%         ycoords = spline(t, ycoords, ts);
% 
%     end

    h = plot(xcoords, ycoords,'Color',color,'LineWidth',thickness);

    out.h = h;
    out.xcoords = xcoords;
    out.ycoords = ycoords;

end % function