function [p, mesh_struct] = isosurface(obj, varargin)
% Create and visualize an isosurface created from the boundaries in an image object.
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
% % Example cortical brain surface
%
% anat = fmri_data(which('keuken_2014_enhanced_for_underlay.img'));
% figure;
% p = isosurface(anat, 'thresh', 140, 'nosmooth');
% set(p, 'FaceAlpha', 1);
% view(132, 6);
% lightRestoreSingle;
%
% Make a brain-bottom cutaway:
% p = isosurface(anat, 'thresh', 140, 'nosmooth', 'zlim', [-Inf 20]);
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

[xlim, ylim, zlim] = deal([-Inf Inf]);

% ------------------------------------------------------
% parse inputs
% ------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            
            case 'thresh',  mythresh = varargin{i + 1}; varargin{i + 1} = [];
                
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

%vsub = subvolume(

% ------------------------------------------------------
% Isosurface
% ------------------------------------------------------

[~, ~, mesh_struct] = reconstruct_image(obj);                % get volume data for slices
voldata = mesh_struct.voldata;

if dosmooth
    V = smooth3(voldata, 'gaussian', mysmoothbox, mygaussstd);
else
    V = voldata;
end

Fvoldata = isosurface(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, V, mythresh);

p = patch('Faces',Fvoldata.faces,'Vertices',Fvoldata.vertices,'FaceColor', mycolor, ...
    'EdgeColor','none','SpecularStrength', .2,'FaceAlpha', .3,'SpecularExponent', 200);

drawnow

% Isocaps, if needed

% pp = [];
isocap = isocaps(mesh_struct.X, mesh_struct.Y, mesh_struct.Z, V, mythresh);
p(end + 1) = patch(isocap, 'FaceColor', 'interp','EdgeColor', 'none', 'FaceAlpha',1);


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