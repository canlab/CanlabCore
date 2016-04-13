function bIsEmpty = tor_ihb_GetClusterSet()
% :Usage:
% ::
%
%    bIsEmpty = ihb_GetClusterSet
%
% Function to select set of cluster (SPM.mat) and select one cluster
% from this set
%
% Return value == 1 if cluster set is nonempty and one cluster selected
%                 0 otherwise
%
% ..
%   10.04.01    Sergey Pakhomov
% ..


bIsEmpty = 1; % Get clusters from SPM.mat
hFigMain = findobj('Type', 'figure', 'Tag', 'ihb_TalSpaceMain_fig');
ihb_HideShowFigureAll('off');
clusters = tor_ihb_GetClusters;
if isempty(clusters)
    bIsEmpty = 0;
    h = msgbox('No voxels above threshold', 'Warning!', 'warn', 'modal');
    waitfor(h);
    ihb_HideShowFigureAll('on');
    return;
end
setappdata(hFigMain, 'clusters', clusters);
%----------------------------------------------------------------------------------
% Select one of the clusters to view
%----------------------------------------------------------------------------------
[cl, indClusterToView] = ihb_selectClusters(clusters, 'single');

if isempty(cl)
    bIsEmpty = 0;
    h = msgbox('No clusters selected','Warning','warn', 'modal');  
    waitfor(h);
    ihb_HideShowFigureAll('on'); 
    return; 
end
%----------------------------------------------------------------------------------
% Set clusters UI combo box
%----------------------------------------------------------------------------------
hPopUp = findobj('Type', 'uicontrol', 'Style', 'popupmenu', 'Tag', 'ihb_ClusterPopUp');
set(hPopUp, 'String', {clusters(:).name}, 'Value', indClusterToView);
%----------------------------------------------------------------------------------
% Load 
%----------------------------------------------------------------------------------
ihb_PrepareClusterToView;
%----------------------------------------------------------------------------------
% Make main figure visible
%----------------------------------------------------------------------------------
ihb_HideShowFigureAll('on');
