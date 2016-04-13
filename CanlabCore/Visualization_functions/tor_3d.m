function [D,Ds,hdr,p,coords,X,Y,Z] = tor_3d(varargin)
% :Usage:
% ::
%
%    [D,Ds,hdr,p,coords,X,Y,Z] = tor_3d(varargin)
%
% made to use single_subj_T1.img from SPM99
%
% :Options:
%
%   **'data':**
%        followed by image data, must also use 'hdr'; data is full image volume
%
%   **'hdr':**
%        followed by hdr structure (use read_hdr)
%
%   **'coords':**
%        3 element row vector of x,y,z coordinates at which to cut away
%
%   **'figure':**
%        create a new figure to plot on
%
%   **'whichcuts':**
%        followed by incisions; choices are x y z w (whole)
%
%        example:'whichcuts','xyz' (xyz is default)
%        order of output handles (p(i)) is 'wyzx'
%
%        New: Special methods:
%          - 'coronal slice right'
%          - 'coronal slice left'
%          - 'coronal slice'
%
%   **'filename':**
%        followed by filename of image data to load, in single quotes
%
%        should be analyze img format, without the .img extension
%
%        cluster imaging assumes neurological orientation, but should work anyway.
%
%   **'revx':**
%        reverse x cut direction, so cut in from left instead of right
%
%   **'topmm':**
%        topmm = varargin{i+1};
%
%   **'intensity_threshold':**
%        percentile of data above which is considered in-object
%        higher = more sparse object
%
% :Outputs:
%
%   **D:**
%        img data
%
%   **Ds:**
%        smoothed data
%
%   **hdr:**
%        header
%
%   **p:**
%        image handles
%
%   **coords:**
%        coordinates
%
%   **X, Y, Z:**
%        are the millimeter, origin centered reference frame for the head isosurface
%
% :Examples:
% ::
%
%    [D,Ds,hdr,p,coords] = tor_3d('figure','data',D,'hdr',hdr,'whichcuts','yzx');
%    [D,Ds,hdr,p,coords] = tor_3d('figure');
%    [D,Ds,hdr,headhandle,coords] = tor_3d('whichcuts','z', 'coords', [-Inf
%    -Inf Inf], 'filename', 'T1_face_exemplar', 'intensity_threshold', 80); set(gcf, 'Color', 'w'); axis on
%
%    % Special slice example:
%    [D,Ds,hdr,handle,coords]  = tor_3d('whichcuts', 'coronal slice right', 'coords', [0 12 0], 'topmm', 100);
%    lightRestoreSingle
%    set(handle(1), 'FaceColor', [.5 .5 .5])
%
% ..
%    Made by Tor Wager, 10/3/2001 last modified 10/19/01
% ..


% ..
%    Set up arguments and default values
% ..
D = [];
Ds =[];
p = zeros(1,2);
coords = [0 0 0];
whichcuts = 'xyz';
filename = 'spm2_single_subj_T1_scalped'; 'single_subj_T1';    % default structural image to use
topmm = 60;                     % image only this high in mm to avoid image distortions above head.
bottommm = nan;                 % bottom cutoff in mm
dolight = 1;
revx = 0;
intensity_threshold = [];

for i = 1:nargin
    if isstr(varargin{i})
        switch varargin{i}
        case 'data', D = varargin{i+1};
        case 'hdr', hdr = varargin{i+1};
        case 'coords', coords = varargin{i+1};
        case 'figure', h3dfig = figure; set(h3dfig,'Tag','myFig');
            colormap(gray(100)),set(gcf,'Color','k'),axis off,set(gcf,'Position',[184   115   672   543])
        case 'whichcuts',whichcuts = varargin{i+1};
        case 'filename',filename = varargin{i+1};
        case 'nolight', dolight = 0;
        case 'revx', revx = 1;
        case 'topmm', topmm = varargin{i+1};
        case 'bottommm', bottommm = varargin{i+1}; 
            
        case 'intensity_threshold', intensity_threshold = varargin{i + 1};
                
        end
    end
end

%------------------------------------------------------------------------------
% Load the image file, if necessary
%------------------------------------------------------------------------------
if isempty(D)
    fullpath = which([filename '.img']);
    if isempty(fullpath), error(['Cannot find file: ' filename '.img']);
    else disp(['Loading structural image: ' fullpath]);drawnow
    end
    [array,hdr] = readim2(filename);
    % rotate 270 degrees so that front of head is positive y
    % (works with canonical SPM images, at least, so this orientation is ok.)
    for i = 1:size(array,3)
        D(:,:,i) = rot90(rot90(rot90((array(:,:,i)))));
    end

end

if isempty(Ds)
    Ds = smooth3(D);
end


%------------------------------------------------------------------------------
% Define X Y Z coordinates of head data
%------------------------------------------------------------------------------
% Define X, Y, Z relative to the origin, so head will be centered on origin
% Multiply by voxel size in header so that coordinates are in mm from the origin
% Make sure origin is set properly in header for this to work accurately.

% NOTE: spm5 no longer sets things in origin, apparently, so use spm_vol
% to get origin
VV = spm_vol(fullpath); 
origin = VV.mat(1:3, 4);
voxsize = (diag(VV.mat(1:3, 1:3)));

[M N P] = size(D);
%[X Y Z] = meshgrid(((1:N)-hdr.origin(1))*hdr.xsize, ((1:M)-hdr.origin(2))*hdr.ysize, ((1:P)-hdr.origin(3))*hdr.zsize);
[X Y Z] = meshgrid(((1:N)*voxsize(1)+origin(1)), ((1:M)*voxsize(2)+origin(2)), ((1:P)*voxsize(3)+origin(3)));

% define threshold
if isempty(intensity_threshold)
    surface_threshold = mean(D(:));  % original way, pre-2009
else
    surface_threshold = prctile(D(:), intensity_threshold);  %mean(mean(mean(D)));
end

%------------------------------------------------------------------------------
% Make patches for cutaway
%------------------------------------------------------------------------------
index = 1;

if strcmp(whichcuts, 'coronal slice right')
    % Special methods
    [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[0 nan coords(2)-2 coords(2) bottommm topmm], surface_threshold);

elseif strcmp(whichcuts, 'coronal slice left')
    [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan 0 coords(2)-10 coords(2) bottommm topmm], surface_threshold);

elseif strcmp(whichcuts, 'coronal slice')
    [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan nan coords(2)-10 coords(2) bottommm topmm], surface_threshold);

else

    if any(whichcuts == 'w')
        [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan nan nan nan bottommm topmm], surface_threshold);          % whole head
        index = index + 2;
    end

    % xyz inset cutaway - seems to work only if x does not come first
    if any(whichcuts == 'y')
        [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan nan nan coords(2) bottommm topmm], surface_threshold);
        index = index + 2;
    end
    if any(whichcuts == 'z')
        [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan nan nan nan bottommm coords(3)], surface_threshold);
        index = index + 2;
    end
    if any(whichcuts == 'x')
        if revx, [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[coords(1) nan nan nan bottommm topmm], surface_threshold);
        else [p(index),p(index+1)] = imagePatch(X,Y,Z,D,Ds,[nan coords(1) nan nan bottommm topmm], surface_threshold);
        end
        index = index + 2;
    end

end

%------------------------------------------------------------------------------
% Set lighting conditions
%------------------------------------------------------------------------------
if dolight && length(p) > 1
    myLight = standardMRIlighting('full',p(1:2));
elseif dolight
    myLight = standardMRIlighting('full',p);
end

if length(p) > 2, standardMRIlighting('reflectance',p(3:4)); end
if length(p) > 4, standardMRIlighting('reflectance',p(5:6)); end

% set(myLight,'Tag','myLight')  % done in stMRIlighting.

rotate3d off
	%------------------------------------------------------------------------------
    % set callback to light follow the camera
	%------------------------------------------------------------------------------
	if exist('lightFollowView') == 2
        set(gcf, 'WindowButtonUpFcn', 'lightFollowView');
    else
        warning('Cannot find lightFollowView.m to set light position.')
    end


%------------------------------------------------------------------------------
% Sub-function for imaging patch
%------------------------------------------------------------------------------
function [p1,p2] = imagePatch(X,Y,Z,D,Ds,inValues, surface_threshold)

    

    if isempty(inValues)
        inValues = [nan nan nan nan nan nan];
    end

    if isempty(X) | isempty(Y) | isempty(Z)
        % no xyz coordinates
        FV2 = isosurface(Ds,surface_threshold);
        IS2 = isocaps(D,surface_threshold);

        % Draw figure patches
        try
            p1 = patch(FV2,'FaceColor',[.8,.5,.4],'EdgeColor','none','FaceAlpha',1,'SpecularExponent',200,'SpecularStrength',.2);
            p2 = patch(IS2,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
        catch
            p1 = patch(FV2,'FaceColor',[1,.75,.65],'EdgeColor','none','SpecularExponent',200,'SpecularStrength',.2);
            p2 = patch(IS2,'FaceColor','interp','EdgeColor','none');
        end
        
    else

        % Define subvolume for cutaway view
        [x y z A] = subvolume(X,Y,Z,D,inValues);
        [x y z As] = subvolume(X,Y,Z,Ds,inValues);
        FV2 = isosurface(x,y,z,As,surface_threshold);
        IS2 = isocaps(x,y,z,A,surface_threshold);

        % Draw figure patches
        try
            p1 = patch(FV2,'FaceColor',[.8,.5,.4],'EdgeColor','none','FaceAlpha',1,'SpecularExponent',200,'SpecularStrength',.2);
            p2 = patch(IS2,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
        catch
            p1 = patch(FV2,'FaceColor',[.8,.5,.4],'EdgeColor','none','FaceAlpha',1,'SpecularExponent',200,'SpecularStrength',.2);
            p2 = patch(IS2,'FaceColor','interp','EdgeColor','none');
        end
        % isonormals(A,p1)
    end

    drawnow
    return



