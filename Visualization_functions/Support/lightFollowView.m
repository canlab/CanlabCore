function lightFollowView
% this uses the current axis and sets the lightangle to follow the camera view. 
% to set, myLight = camlight(0,0); set(myLight,'Tag','myLight')
% call with:
% set(gcf, 'WindowButtonUpFcn', 'lightFollowView');


    %rotate3d('up'); % call default processing (3D rotation enabled)
    % Matlab 7.2 compat. : changed by Tor, 5/22/06
    rotate3d on
    
    %----------------------------------------------------------------------------------
    % Set light position into current camera position
    %----------------------------------------------------------------------------------

    set(0,'CurrentFigure', gcf)
    set(gcf, 'CurrentAxes', gca);
    hLight = findobj('Type', 'light', 'Tag', 'myLight');
    set(hLight, 'Position', get(gca, 'CameraPosition'));
    
return
    