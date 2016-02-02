function hh = nmdsfig_fill(varargin)
% Purpose: to take information about objects in multidimensional space
% and draw colored contours around them.
% Used within nmdsfig.m
%
% Usage:
%
% 1) If c is a structure from cluster_nmdsfig, which is compatible with
% nmdsfig_tools, then:
% ::
%
%    hh = nmdsfig_fill(c)
%
% Fields used are classes, groupSpace, and colors (see code for
% details)
%
% 2) Pass in arguments directly:
% ::
%
%    hh = nmdsfig_fill(classes, positions, colors)
%
% classes is n x 1 vector of group assignments (or all ones for one
% group)
%
% positions is n x 2 matrix of x, y coordinates
% colors is a cell containing as many colors as classes, {'r' 'g' 'b' ...}
% or {[1 0 0] [0 1 0] [0 0 1] ...}
%
% :Examples:
% ::
%
%    load nmdsfig_output
%    hh = nmdsfig_fill(c)
%    set(findobj('Type','Line'), 'Color', 'k')
%
% ..
%    by tor wager, feb 07
%
%    Other notes:
%    fills area around coordinates of regions in a class (group) in a color
%    using spline interpolation and other stuff.
%    designed for cluster imaging in nmdsfig figures.
%    Uses fill_area_around_points, within which
%    Many choices can be set that control how fills are done.
%    tor recommends .2 for borderscale...
% ..

hh = [];
if nargin == 0, help nmdsfig_fill, return, end

if isstruct(varargin{1})
    c = varargin{1};
    clas = c.ClusterSolution.classes;
    positions = c.GroupSpace;
    
    if isfield(c, 'colors')
        colors = c.colors;
    else
        colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
    end
else
    clas = varargin{1};
    positions = varargin{2};
    colors = varargin{3};
end


   for i = 1:max(clas)
       
       coords = positions(clas == i, 1:2);
       
       if ischar(colors{i})
        h = fill_area_around_points(coords(:,1), coords(:,2), .2, colors{i}(1));
       else
           h = fill_area_around_points(coords(:,1), coords(:,2), .2, colors{i});
       end
       
       if ~isempty(h) && all(ishandle(h))
           delete(h(2)) % line
       
            hh(i) = h(1);
       end
       
   end
   
   drawnow

end
