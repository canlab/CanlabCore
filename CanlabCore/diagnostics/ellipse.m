function [x,y] = ellipse(x,v1,v2,c,varargin)
% Gives x and y coordinates for an ellipse, given x coordinates,
% at a distance of c
%
% :Usage:
% ::
%
%     [x,y] = ellipse(x,v1,v2,c,[sorting method])
%
% :Inputs:
%
%   Based on the formula for an ellipse, x^2/v1^2 + y^2/v2^2 = c
%
%   **c:**
%        is the distance from the origin
%
%   **v1:**
%        is the x half-length
%
%   **v2:**
%        is the y half-length
%
%   **x:**
%        is a vector of points covering the x coordinates in the ellipse
%
% :Sorting methods: 
%   - sort by x, produces elliptical line in plot
%   - sort by y, produces horizontal lines in plot
%
% :Examples:
% ::
%
%    [x,y]=ellipse((randn(1000,1)),1.5,2.5,1); figure; hh = plot(x(2:end-1),y(2:end-1),'r-')
%    rotate(hh,[0 90],45)  % rotate around z-axis by 45 degrees
%    x2 = get(hh,'XData'); y2 = get(hh,'YData'); hold on; plot(x2,y2,'bx');
%    rotate(hh,[0 90],-45)  % rotate original ellipse back
%
%    % fill
%    fill(x,y,'r','FaceAlpha',.2)
%
% ..
%    tor wager
% ..

if length(varargin) > 0, sortby = varargin{1};, else, sortby = 1;, end
if size(x,2) > size(x,1), x = x';, end
if sortby == 1, x = sort(x);, end


y = ( v2.^2 .* (c - x.^2 ./ v1.^2) ) .^ .5;

wh = imag(y) ~= 0;
x(wh) = [];
y(wh) = [];

if sortby == 1,     % all this is so it makes a nice line plot
    xy = [x y (length(x):-1:1)'];
    xy = sortrows(xy,3);
    x2 = xy(:,1);
    y2 = xy(:,2);
end

y = [y2;-y];
x = [x2;x];

xy = [x y];

switch sortby
case 1
    %xy = sortrows(xy,1);
    % sort first!
case 2
    xy = sortrows(xy,2);
otherwise
    %xy = sort(xy,1);    % this produces a sigmoid function!
end

x = xy(:,1); y = xy(:,2);

% complete the ellipse
y = [y; y(1)];
x = [x; x(1)];

return
