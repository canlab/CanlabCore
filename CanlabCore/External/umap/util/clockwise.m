%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [x, y]=clockwise(xy)
x=xy(:,1);
y=xy(:,2);
cx = mean(x);
cy = mean(y);
a = atan2(y - cy, x - cx);
[~, order] = sort(a);
x = x(order);
y = y(order);
%assert(all(ismember([x y],xy,'rows')));
end
