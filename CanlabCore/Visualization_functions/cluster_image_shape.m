function [hpatch, cl] = cluster_image_shape(cl, varargin)
% :Usage:
% ::
%
%    [hpatch, cl] = cluster_image_sphere(cl or [k x 3 list of mm coords], varargin)
%
% Images spheres at cluster centers
% Combine with addbrain and cluster_tools('connect3d') 
% or cluster_nmdsfig_glassbrain
% to get 3-D plots of connected blobs
%
% :Optional Inputs:
%
%   {'color', 'colors'}, mycolor = varargin{i+1};
%
%   'radius', myradius = varargin{i+1};
%
% :Outputs:
%
%   **hpatch:**
%        patch handles
%
%   **cl:**
%        new cl with sphere coordiates in XYZmm and XYZ
%
% :Examples:
% ::
%
%    function [hpatch, cl] = cluster_image_sphere(cl)
%
%    % With optional arguments:
%    [hpatch, cl] = cluster_image_sphere(cl, 'color', 'b', 'radius', 10)
%    [hpatch, cl] = cluster_image_sphere(cl, 'color', {'r' 'g' 'b' etc}, 'radius', 10)
%
% :Example: Given an MNI coordinate, plot a sphere on a brain surface
% ::
%
%    my_mm_coord = [40, 46, 22]';
%    create_figure('surface')
%    cl = [];
%    cl.XYZmm = my_mm_coord;
%    cl.mm_center = my_mm_coord';
%    V = spm_vol(which('brainmask.nii'));
%    cl.M = V.mat;
%    [hpatch, cl] = cluster_image_sphere(cl, 'color', 'g', 'radius', 10)
%    p = addbrain;
%    set(p, 'FaceAlpha', 1);
%    axis image
%    view(135, 30);
%    lighting gouraud;
%    lightRestoreSingle;
%    material dull;
%
%    % Turn xyz mm coordinates into clusters and image them
%
%    my_mm_coord = [40 46 22; 50 26 40; 45 36 50; 60 12 0]
%    [hpatch, cl] = cluster_image_sphere(my_mm_coord, 'color', 'b', 'radius', 4);
%
% ..
%    Tor Wager, July 2007 (original version)
%
%    Programmers' notes:
%    Tor Wager, July 2007
%    updated April 2011 for flexible radius
%    updated 12/2012 for xyz to spheres
% ..

hpatch = [];
mycolor = 'r';
myradius = 8;
shape = 'sphere';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case {'color', 'colors'}, mycolor = varargin{i+1}; varargin{i+1} = [];
            case 'cube', shape = 'cube';
            case 'sphere', shape = 'sphere';
            case 'radius', myradius = varargin{i+1};
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if isa(cl, 'region')
    cl = region2struct(cl);
end

% if cl is actually a series of xyz coordinates...
if ~isstruct(cl)
    cl = create_cl_struct(cl);
end

if length(myradius) == 1
    myradius = repmat(myradius, 1, length(cl));
end

% get colors for each sphere
if length(mycolor) == 1 || ~iscell(mycolor)
    mycolor = repmat({mycolor}, 1, length(cl));
end

%newcl = convert_cl_xyz_to_sphere(cl(1), myradius);
for i = 1:length(cl)
    if strcmp(shape,'sphere')
        newcl(i) = convert_cl_xyz_to_sphere(cl(i), myradius(i));
    elseif strcmp(shape,'cube')
        newcl(i) = convert_cl_xyz_to_cube(cl(i), myradius(i));
    else
        error('programmers error: shape = ''%s''',shape)
    end
    
    hpatch(i) = image_isosurface(newcl(i), mycolor{i});
end

cl = newcl;

end


function cl = create_cl_struct(xyz)

V = spm_vol(which('brainmask.nii'));
if size(xyz, 2) ~= 3, xyz = xyz'; end
if size(xyz, 2) ~= 3, error('Bad input xyz coordinates?'); end

cl = [];
for i = 1:size(xyz, 1)
    my_mm_coord = xyz(i, :);
    cl(i).XYZmm = my_mm_coord';
    cl(i).mm_center = my_mm_coord;
    cl(i).M = V.mat;
end

end


function cl = convert_cl_xyz_to_sphere(cl, myradius)

    newXYZmm = round(points_in_sphere(cl(1).mm_center, myradius)');

    cl(1).Z = ones(1, size(newXYZmm,2));
    cl(1).XYZmm = newXYZmm;
    cl(1).XYZ = mm2voxel(cl(1).XYZmm, cl(1).M, 1)';
    cl(1).numVox = size(cl(1).XYZ, 2);

end


function cl = convert_cl_xyz_to_cube(cl, myradius)

c = round(cl(1).mm_center);
r = round(myradius);
twor = 2 * r;
newXYZmm = [];
for x = c(1) - r : c(1) + twor
    for y = c(2) - r : c(2) + twor
        for z = c(3) - r : c(3) + twor
            newXYZmm = [newXYZmm [x y z]'];
        end
    end
end

cl(1).Z = ones(1, size(newXYZmm,2));
cl(1).XYZmm = newXYZmm;
cl(1).XYZ = mm2voxel(cl(1).XYZmm, cl(1).M, 1)';
cl(1).numVox = size(cl(1).XYZ, 2);

end


function hpatch = image_isosurface(cl, mycolor)

    % controls padding to make sure we cover whole area
    padval = 1;
        
    % controls smoothness, etc.
    mythresh = .5;
    mysmoothbox = 3;
    mygaussstd = 1;

    mm_coords = cl(1).XYZmm;

    xyzmin = min(mm_coords') - padval;     % minus/plus for padding
    xyzmax = max(mm_coords') + padval;

    [X, Y, Z] = meshgrid(xyzmin(1):xyzmax(1), xyzmin(2):xyzmax(2), xyzmin(3):xyzmax(3));


    % construct volume data for area to image

    xvox = mm_coords(1,:) - xyzmin(1) + 1;
    yvox = mm_coords(2,:) - xyzmin(2) + 1;
    zvox = mm_coords(3,:) - xyzmin(3) + 1;

    V = zeros(size(X));

    for i = 1:size(xvox, 2)
        V(xvox(i), yvox(i), zvox(i)) = cl(1).Z(i);
    end

    % not needed if we have all mm points in sphere
    %VI = interp3(cl(1).XYZmm(1,:), cl(1).XYZmm(2,:), cl(1).XYZmm(3,:), cl(1).Z, X, Y, Z);

    V = smooth3(V, 'gaussian', mysmoothbox, mygaussstd);
    FV = isosurface(X,Y,Z,V, mythresh);

    hpatch = patch(FV);
    set(hpatch, 'EdgeColor', 'none', 'FaceColor', mycolor);

    drawnow

    %lighting gouraud

end
