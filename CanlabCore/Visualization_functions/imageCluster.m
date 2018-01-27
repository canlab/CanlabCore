function [out,cl] = imageCluster(varargin)
% :Usage:
% ::
%
%    [out,cl] = imageCluster(arguments as specified below)
%
% Images a cluster isosurface on an existing 3D head or brain plot
%
% works with tsu (Talairach Space Utility) main window as current figure
% old: clusters = getappdata(gcf, 'clusters');
% 
% Inputs (in any order): keyword followed by input argument
%
%   **'cluster':**
%        followed by cluster to image, from SPM or TSU.
%
%   **'getclusters':**
%        no other args necessary - starts gui for cluster selection
%        function returns all clusters.  select with clusters(i)
%
%   **'getfigclusters':**
%        get clusters from TSU main figure.  must be current figure.
%
%   **'figure':**
%        create a new figure to image color on
%
%   **'color':**
%        followed by color value - either text or vector
%
%   **'alpha':**
%        followed by transparency value for cluster surface, 0-1
%        1 is opaque, 0 is completely transparent
% 
% :Control of smoothing:
%
%   **'heightthresh':**
%        followed by cutoff threshold post-smooth, in percentage of min Z value in cl
%          - enter a number between 0 and 1
%
%   **'fwhm':**
%        followed by smoothing kernel FWHM (Gaussian)
%
%   **'kernelsize':**
%        followed by box size for kernel support (5 5 5 is  default)
%
% :Outputs:
%
%   **out:**
%        Patch handle
%
%   **cl:**
%        cluster struct
%
% :Uses: cl.XYZmm and cl.voxSize
%
% Works in Matlab 5.3, but with no transparency.
%
% :Example:
% ::
%
%    p(i) = imageCluster('cluster',region2struct(r(i)),'color',colors{i},'alpha',1, 'fwhm', 1.2, 'heightthresh', .3); 
%    view(135, 30); lighting gouraud; lightRestoreSingle; axis image; camlight right;
%
% ..
%    By Tor Wager, 10/3/2001, edited 7/2012.  1/18: removed verbose
%    printout
% ..

% ..
%    Set up arguments and default values
% ..
out = [];
clusters = [];
mycolor = 'y';
myalpha = 1;
heightthresh = .5;
fwhm = 1;
kernelsize = [5 5 5];

for i = 1:nargin
    if isstr(varargin{i})
        switch varargin{i}
        case 'cluster', cl = varargin{i+1};
        case 'getclusters', clusters = tor_ihb_getClusters;
        case 'getfigclusters', clusters = getappdata(gcf, 'clusters');
        case 'figure', h3dfig = figure; set(h3dfig,'Tag','myFig');
        case 'color',mycolor = varargin{i+1};
        case 'alpha',myalpha = varargin{i+1};
            
            case 'heightthresh', heightthresh = varargin{i+1};
            case 'fwhm', fwhm = varargin{i+1};
            case 'kernelsize', kernelsize = varargin{i+1};
                
        end
    end
end

% If you chose to get clusters, no imaging done - just returns all clusters.
if ~isempty(clusters), out = clusters; return, end
 
if ~exist('cl', 'var'), cl = []; disp('No clusters entered. Enter ''cluster'' followed by clusters.'); return, end

% do not do this [X,Y,Z] = meshgrid(x,y,z);

%------------------------------------------------------------------------------
% Prepare volume data and create isosurface patch object
%------------------------------------------------------------------------------
if isa(cl, 'region'), cl = region2struct(cl); end

% if a vector of clusters, concatenate
if length(cl) > 1
    cl = clusters2CLU(cl); 
end
    
    try
        %V = smooth3(cl.vTal);
        %V = smooth3(cl.vTal,'box',5);
        V = smooth3(cl.vTal,'gaussian', kernelsize, fwhm);

    catch
	% no vTal field, make it here.
	% disp(['No vTal field for cluster, attempting to create it.'])
	cl = get_cluster_volume(cl);
    %V = smooth3(cl.vTal);
	% V = smooth3(cl.vTal,'box',5);
    
    V = smooth3(cl.vTal,'gaussian', kernelsize, fwhm);
    
    %V = cl.vTal;
    end

% 	if isfield(cl,'hThreshold'), height_threshold = cl.hThreshold;
% 	elseif isfield(cl,'threshold'), height_threshold = cl.threshold;
% 	else height_threshold = min(cl.Z); height_threshold; %height_threshold = 1;,warning('No height threshold found in cluster.')
% 	end
    %height_threshold = floor(min(cl.Z) ./2); , height_threshold = height_threshold(1);
 
    cl.Z = double(cl.Z); % fix bug with single format
    height_threshold = (min(cl.Z) .* heightthresh); 
    height_threshold = height_threshold(1);
%height_threshold = 3

    % disp(['Range of cl is: ' num2str(min(cl.Z)) ',' num2str(max(cl.Z)) '. Threshold is :' num2str(height_threshold)])
    
 	set(0,'CurrentFigure', gcf)
	%set(gcf, 'CurrentAxes', gca);

	FV = isosurface(cl.xTal, cl.yTal, cl.zTal, V, height_threshold);
	hPatch = patch(FV);
	isonormals(cl.xTal, cl.yTal, cl.zTal, cl.vTal, hPatch);
	try
		set(hPatch, 'Tag', 'ihb_surfaceCluster', 'FaceColor', mycolor, 'EdgeColor', 'none', 'FaceAlpha', myalpha);
	catch
		set(hPatch, 'Tag', 'ihb_surfaceCluster', 'FaceColor', mycolor, 'EdgeColor', 'none');
	end

    %delete(hMsg);
    set(gcf, 'Visible', 'on');
    out = hPatch;
    
end % function
    
