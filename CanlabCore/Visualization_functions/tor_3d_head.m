% run this script AFTER using talairach space utility (tsu.m) to get clusters
% this copies the orignial tsu stuff in ihb_TalSpaceView3D.m, but in the workspace

%------------------------------------------------------------------------------
% Get cluster information
%------------------------------------------------------------------------------
clusters = getappdata(hFigMain, 'clusters');
indClusterToView = getappdata(hFigMain, 'indClusterToView');
cl = clusters(indClusterToView);
iv = cl.hThreshold;


%------------------------------------------------------------------------------
% Figure does not exist - create new one
%------------------------------------------------------------------------------
	hFig3D = figure(...
		'Visible', 'off',...
		'NumberTitle', 'off',...
		'MenuBar', 'none',...
		'Name', '3D view of cluster',...
        'Renderer', ihbdfl_renderer,...
		'Tag', 'ihb_TalSpace3D_fig'...
		);

	%------------------------------------------------------------------------------
    % Create and/or adjust light position
	%------------------------------------------------------------------------------
    hLight = findobj('Type', 'light', 'Tag', 'ihb_lightTalSpace3D');
    cameraPos = get(hAxis3D, 'CameraPosition');
    if isempty(hLight)
	    hLight = light('Tag', 'ihb_lightTalSpace3D', 'Position', cameraPos);
    else
        set(hLight, 'Position', cameraPos);  
    end
	%------------------------------------------------------------------------------
    % set callback to light follow the camera
	%------------------------------------------------------------------------------
	set(hFig3D, 'WindowButtonUpFcn', 'ihb_View3DButtonUpFcn');
	view([-150,32]);    % set initial camera position
	lighting phong;     % set lighting algorithm
	rotate3d on;        % enable 3D rotation

    %------------------------------------------------------------------------------
    % Smooth volume data
	%------------------------------------------------------------------------------
    if ihbdfl_bApplySmoothing == 1
        switch ihbdfl_smooth_filter
        case 'box'
            V = smooth3(cl.vTal, 'box', ihbdfl_smooth_size);
        case 'gaussian'
            V = smooth3(cl.vTal, 'gaussian', ihbdfl_smooth_size, ihbdfl_smooth_sd');
        end
    else
	    V = cl.vTal;
    end
    
    %------------------------------------------------------------------------------
    % Prepare volume data and create isosurface patch object
	%------------------------------------------------------------------------------
    % V is the volume of the cluster.
    
    FV = isosurface(cl.xTal, cl.yTal, cl.zTal, V, cl.hThreshold);
	hPatch = patch(FV);
	isonormals(cl.xTal, cl.yTal, cl.zTal, cl.vTal, hPatch);
	set(hPatch, 'Tag', 'ihb_surfaceCluster', 'FaceColor', ihbdfl_color_surface, 'EdgeColor', 'none', 'FaceAlpha', ihbdfl_transparency);
    set(hFig3D, 'Visible', 'on');
    