function [out,cl] = imageCluster_block(varargin)
% :Usage:
% ::
%
%    [out,cl] = imageCluster_block(input arguments)
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
% :Outputs:
%
%   **out:**
%        Patch handle
%
%   **cl:**
%        cluster struct
%
% ..
%    Works in Matlab 5.3, but with no transparency.
%    By Tor Wager, 10/3/2001, last edit 10/18/01
% ..

% ..
%    Set up arguments and default values
% ..
clusters = [];
mycolor = 'y';
myalpha = 1;

for i = 1:nargin
    if isstr(varargin{i})
        switch varargin{i}
        case 'cluster', cl = varargin{i+1};
        case 'getclusters', clusters = tor_ihb_getClusters;
        case 'getfigclusters', clusters = getappdata(gcf, 'clusters');
        case 'figure', h3dfig = figure; set(h3dfig,'Tag','myFig');
        case 'color',mycolor = varargin{i+1};
        case 'alpha',myalpha = varargin{i+1};
        end
    end
end

% If you chose to get clusters, no imaging done - just returns all clusters.
if ~isempty(clusters), out = clusters;, return, end
 

% do not do this [X,Y,Z] = meshgrid(x,y,z);

%------------------------------------------------------------------------------
% Prepare volume data and create isosurface patch object
%------------------------------------------------------------------------------
    try
	    V = cl.vTal;
        %V = smooth3(cl.vTal,'gaussian',[5 5 5]);
        V = smooth3(cl.vTal);
    catch
	    % no vTal field, make it here.
	    disp(['No vTal field for cluster, attempting to create it.'])
	    cl = get_cluster_volume(cl);
	    V = cl.vTal;
        %V = smooth3(cl.vTal,'gaussian',[5 5 5]);
        V = smooth3(cl.vTal);
    end
 
	if isfield(cl,'hThreshold'), height_threshold = cl.hThreshold;
	elseif isfield(cl,'threshold'), height_threshold = cl.threshold;
	else height_threshold = 0;,warning('No height threshold found in cluster.')
	end
    
    %height_threshold = min(cl.Z);
    % to image the whole cluster.  Smoothing will increase the cluster size if you choose this!
    height_threshold = 0;   % set this to 0 because all above 0 should be imaged.
    
    % use with smoothing, to smooth out the clusters and keep them a reasonable size
    %height_threshold = .5 * mean(mean(mean(V)));
    
    
 	set(0,'CurrentFigure', gcf)
	%set(gcf, 'CurrentAxes', gca);
      
	FV = isosurface(cl.xTal, cl.yTal, cl.zTal, V, height_threshold);
	hPatch = patch(FV);
    IC = isocaps(cl.xTal, cl.yTal, cl.zTal, V, height_threshold);
    iPatch = patch(IC);
	isonormals(cl.xTal, cl.yTal, cl.zTal, cl.vTal, hPatch);
	try
		set(hPatch, 'Tag', 'ihb_surfaceCluster', 'FaceColor', mycolor, 'EdgeColor', 'none', 'FaceAlpha', myalpha);
        set(iPatch, 'Tag', 'ihb_surfaceCluster', 'FaceColor', mycolor, 'EdgeColor', 'none', 'FaceAlpha', myalpha);
	catch
		set(hPatch, 'Tag', 'ihb_surfaceCluster', 'FaceColor', mycolor, 'EdgeColor', 'none');
        set(iPatch, 'Tag', 'ihb_surfaceCluster', 'FaceColor', mycolor, 'EdgeColor', 'none');
	end

    %delete(hMsg);
    set(gcf, 'Visible', 'on');
    out = hPatch;
    
return
    
