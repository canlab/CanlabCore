function [p, mesh_struct, my_isosurface, my_isocap] = isosurface(obj, varargin)
% Create and visualize an isosurface created from the boundaries in an image object.
%
% [p, mesh_struct, my_isosurface, my_isocap] = isosurface(obj, varargin)
%
% - p is a patch handle to graphics patch object showing surface and isocaps
% - mesh_struct has meshgrid with x, y, z coordinates in mm, and volume data
% - my_isosurface and my_isocap are structures with .faces and .vertices,
%   which you can visualize using patch()
%
% Options:
% case 'thresh',  mythresh = varargin{i + 1}; varargin{i + 1} = [];
% 
% case 'nosmooth', dosmooth = false;
% case 'smoothbox', mysmoothbox = varargin{i + 1}; varargin{i + 1} = [];
% case 'sd', mygaussstd = varargin{i + 1}; varargin{i + 1} = [];
% 
% case 'color', mycolor = varargin{i + 1}; varargin{i + 1} = [];
% 
% 
% case 'xlim', xlim = varargin{i + 1}; varargin{i + 1} = [];
% case 'ylim', ylim = varargin{i + 1}; varargin{i + 1} = [];
% case 'zlim', zlim = varargin{i + 1}; varargin{i + 1} = [];
% 
% Examples:
% % ------------------------------------------------------
%
% % An example cortical brain surface
% ------------------------------------------------------------------------
%
% anat = fmri_data(which('keuken_2014_enhanced_for_underlay.img'), 'noverbose');
% figure;
% p = isosurface(anat, 'thresh', 140, 'nosmooth');
% set(p, 'FaceAlpha', 1);
% view(132, 6);
% lightRestoreSingle;
%
% Make a brain-bottom cutaway:
% ------------------------------------------------------------------------
%
% create_figure('cutaway');
% p = isosurface(anat, 'thresh', 140, 'nosmooth', 'zlim', [-Inf 20]);
% view(132, 30);
%
% A coronal cutaway around the nucleus accumbens:
% ------------------------------------------------------------------------
%
% create_figure('cutaway');
% p = isosurface(anat, 'thresh', 140, 'nosmooth', 'ylim', [-Inf 10]);
% view(132, 30);
%
% An often-used CANlab 3-d cutaway
% ------------------------------------------------------------------------
% create_figure('cutaway');
% anat = fmri_data(which('keuken_2014_enhanced_for_underlay.img'), 'noverbose');
% p = isosurface(anat, 'thresh', 140, 'nosmooth', 'ylim', [-Inf -30]);
% p2 = isosurface(anat, 'thresh', 140, 'nosmooth', 'xlim', [-Inf 0], 'YLim', [-30 Inf]);
% alpha 1 ; lightRestoreSingle; view(135, 30); colormap gray;
% p3 = addbrain('limbic hires');
% set(p3, 'FaceAlpha', .6, 'FaceColor', [.5 .5 .5]);
% delete(p3(3)); p3(3) = [];
% lightRestoreSingle;
% surface_handles = [p p2 p3];
%
% Another example using a different group average image:
% ------------------------------------------------------------------------
% dat = fmri_data(which('icbm152_2009_symm_enhanced_for_cutaways.nii'), 'noverbose');
% figure;
% [p, mesh_struct] = isosurface(dat, 'percentagethresh', 80);
% alpha 1 ; camlight, lightRestoreSingle; view(135, 30); colormap gray
%
% See also: tor_3d, cluster_cutaway, cluster_image_shape

% ------------------------------------------------------
% defaults
% ------------------------------------------------------

dosmooth = true;
mycolor = [.5 .5 .5];
mysmoothbox = 3;
mygaussstd = 1;
mythresh = 0;
percentagethresh = [];

[xlim, ylim, zlim] = deal([-Inf Inf]);

% ------------------------------------------------------
% parse inputs
% ------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            % reserved keywords
            
            case 'thresh',  mythresh = varargin{i + 1}; varargin{i + 1} = [];
                
            case 'percentagethresh',  percentagethresh = varargin{i + 1}; varargin{i + 1} = [];
                    
            case 'nosmooth', dosmooth = false;
            case 'smoothbox', mysmoothbox = varargin{i + 1}; varargin{i + 1} = [];
            case 'sd', mygaussstd = varargin{i + 1}; varargin{i + 1} = [];
                
            case 'color', mycolor = varargin{i + 1}; varargin{i + 1} = [];
                
                
            case 'xlim', xlim = varargin{i + 1}; varargin{i + 1} = [];
            case 'ylim', ylim = varargin{i + 1}; varargin{i + 1} = [];
            case 'zlim', zlim = varargin{i + 1}; varargin{i + 1} = [];
                
                
            case 'noverbose'
                doverbose = false;
                
                %otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% ------------------------------------------------------
% Get cutaways
% ------------------------------------------------------

x = obj.volInfo.xyzlist(:, 1);
y = obj.volInfo.xyzlist(:, 2);
z = obj.volInfo.xyzlist(:, 3);

% Convert limits to voxels
zeros = [0 0];
xlimvox = mm2voxel([xlim; zeros; zeros], obj.volInfo.mat);
xlimvox = xlimvox(:, 1)';

ylimvox = mm2voxel([zeros; ylim; zeros], obj.volInfo.mat);
ylimvox = ylimvox(:, 2)';

zlimvox = mm2voxel([zeros; zeros; zlim], obj.volInfo.mat);
zlimvox = zlimvox(:, 3)';

% kludge in case x is flipped
if diff(xlimvox) < 0, xlimvox = xlimvox(end:-1:1); end

% Exclude relevant voxels; see also subvolume.m
whx = x < xlimvox(1) | x > xlimvox(2);
why = y < ylimvox(1) | y > ylimvox(2);
whz = z < zlimvox(1) | z > zlimvox(2);

obj.dat(whx | why | whz, :) = 0;


% ------------------------------------------------------
% Isosurface
% ------------------------------------------------------

[~, ~, mesh_struct] = reconstruct_image(obj);                % get volume data for slices
voldata = mesh_struct.voldata;

% If we want isocaps, we must cut volume off - otherwise it will be solid surface
v0 = voldata == 0;

wh_z = ~squeeze(all(all(v0, 1), 2));
zl = [min(find(wh_z)) max(find(wh_z))];

wh_y = ~squeeze(all(all(v0, 1), 3));
yl = [min(find(wh_y)) max(find(wh_y))];

wh_x = ~squeeze(all(all(v0, 2), 3));
xl = [min(find(wh_x)) max(find(wh_x))];

mesh_struct.X = mesh_struct.X(wh_x, wh_y, wh_z);
mesh_struct.Y = mesh_struct.Y(wh_x, wh_y, wh_z);
mesh_struct.Z = mesh_struct.Z(wh_x, wh_y, wh_z);

voldata = voldata(wh_x, wh_y, wh_z);

% mylimits = [xl yl zl];
% [mesh_struct.X, mesh_struct.Y, mesh_struct.Z, voldata] = subvolume(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, voldata, mylimits);
%     
    
if dosmooth
    V = smooth3(voldata, 'gaussian', mysmoothbox, mygaussstd);
else
    V = voldata;
end

if ~isempty(percentagethresh)
    mythresh = prctile(V(:), percentagethresh);
end

my_isosurface = isosurface(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, V, mythresh);

p = patch('Faces',my_isosurface.faces,'Vertices',my_isosurface.vertices,'FaceColor', mycolor, ...
    'EdgeColor','none','SpecularStrength', .2,'FaceAlpha', .3,'SpecularExponent', 200);

drawnow

% Isocaps, if needed

% pp = [];
my_isocap = isocaps(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, V, mythresh);
p(end + 1) = patch(my_isocap, 'FaceColor', 'interp','EdgeColor', 'none', 'FaceAlpha',1);


% ------------------------------------------------------
% Colors and lighting
% ------------------------------------------------------

%set(p, 'FaceColor', mycolor);

lighting gouraud
camlight right
camlight left
material dull

axis image
axis vis3d

drawnow

end