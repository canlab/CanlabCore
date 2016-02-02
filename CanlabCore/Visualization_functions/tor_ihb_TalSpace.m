function tor_ihb_TalSpace()
% :Usage:
% ::
%
%    tor_ihb_TalSpace
%
% Main function to view SPM99 cluster's (contours & 3D) in Talariach space
%
%==================================================================================
%   Graphic objects and their 'Tags':
%==================================================================================
%   type        handle          Tag                     Description             
%----------------------------------------------------------------------------------
%   figure      hFigMain        ihb_TalSpaceMain_fig    Main figure    
%   figure      hFig3D          ihb_TalSpace3D_fig      3D view figure
%   figure      hFigLeg         ihb_TalSpaceLeg_fig     Legend figure
%
%   axes        hAxis3D         ihb_axesTalSpace3D      Axes for 3D view
%   light       hLight          ihb_lightTalSpace3D     Light object for 3D view
%
%   patch                       ihb_cntAxial            contours patch handles
%   patch                       ihb_cntFrontal              in 3D view window
%   patch                       ihb_cntSaggitalL
%   patch                       ihb_cntSaggitalR
% 
%   patch                       ihb_surfaceCluster      isosurface patch
%
%   uicontrol                   ihb_CursorCoord         coordinates under cursor
%   uicontrol                   ihb_ClusterLevelProb    cluster probability
%   uicontrol                   ihb_VoxelLevelProb      extr. probability
%   uicontrol                   ihb_VolumeInVox         cluster volume
%   uicontrol                   ihb_ClusterPopUp        select cluster combo box
%
%   uicontrol                   ihb_clr_Orig            pushbuttons for
%   uicontrol                   ihb_clr_Symm            color change in
%   uicontrol                   ihb_clr_Surface         legend window
%   uicontrol                   ihb_clr_Axial
%   uicontrol                   ihb_clr_Frontal
%   uicontrol                   ihb_clr_Saggital
%
%   line                        ihb_LineOriginal        lines for original clusters
%   line                        ihb_LineSymmetrical     lines for symmetrical clusters
%
%   uimenu                      ihb_VoxSizeMenu111
%   uimenu                      ihb_VoxSizeMenu222
%   uimenu                      ihb_VoxSizeMenu444
%
%   uimenu                      ihb_RendererMenuZbuffer
%   uimenu                      ihb_RendererMenuOpenGL
%
%   uimenu                      ihb_TransparencyMenu
%   uimenu                      ihb_TransparencyMenuOpaque
%   uimenu                      ihb_TransparencyMenuLow
%   uimenu                      ihb_TransparencyMenuMedium
%   uimenu                      ihb_TransparencyMenuHigh
%
%   uimenu                      ihb_Axes3DMenuClusterOnly
%   uimenu                      ihb_Axes3DMenuWholeBrain
% 
%   uimenu                      ihb_LineWidthMenuThin
%   uimenu                      ihb_LineWidthMenuNormal
%   uimenu                      ihb_LineWidthMenuThick
%
%==================================================================================
%  Application data (name used with setappdata the same as handle or variable name
%==================================================================================
% for hFigMain:
%
%   axes            hAxAxial            axial               
%   axes            hAxFront            frontal      
%   axes            hAxSagL             left saggital       
%   axes            hAxSagR             right saggital       
%   strucure array  clusters            array of selected clusters
%   integer         indClusterToView    index of currently viewed cluster in clusters array
%   bool            drawSymLS           == 1 if draw symmetrical on left saggital
%   bool            drawSymRS           == 1 if draw symmetrical on right saggital
%
%----------------------------------------------------------------------------------
%
% for hAxis (where hAxis is one of hAxAxial, hAxFront, hAxSagL, hAxSagR):
%   (see ihb_LoadSlice for detail)
%
%   strucure array  sliceInfo           array of Talarich slice information (ihbdfl_ax_info etc.)  
%   integer array   rngArray            array of local indexes in sliceInfo
%   integer         sliceIndex          index of slice in the local range array rngArray
%   double          sliceDist           slice distance from origin
%   string          axType              axis type: 'ax' 'fr' 'sl' 'sr' 
%   1x2 vector      xLim                limit values for X direction
%   1x2 vector      yLim                limit values for Y direction
%
%
%==================================================================================
% Default information in ihb_TalSpaceDfl.mat file (see ihb_ResetDefaults)
%==================================================================================
% To get default variable with name var use:
%
% load('ihb_TalSpaceDfl.mat', 'var');
%
% To save:
%
% save('ihb_TalSpaceDfl.mat', 'var', '-append');
%----------------------------------------------------------------------------------
% Figures tags for windows created by SPM during call spm_getSPM
%----------------------------------------------------------------------------------
% ihbdfl_spm_fig_Interactive
% ihbdfl_spm_fig_SelFileWin
% ihbdfl_spm_fig_ConMan
%----------------------------------------------------------------------------------
% Minimal size of clusters to draw
%----------------------------------------------------------------------------------
% ihbdfl_min_size_to_draw
%----------------------------------------------------------------------------------
% Sizes of voxels used in Talariach space
%----------------------------------------------------------------------------------
% ihbdfl_tal_x_vox
% ihbdfl_tal_y_vox
% ihbdfl_tal_z_vox
%----------------------------------------------------------------------------------
% Sizes of main window and axis for ihb_TalSpace
%----------------------------------------------------------------------------------
% ihbdfl_main_width
% ihbdfl_main_height
% ihbdfl_main_gap
% ihbdfl_main_x
% ihbdfl_main_y
% ihbdfl_main_z
% ihbdfl_main_info_h
%----------------------------------------------------------------------------------
% Talariach slices information
%----------------------------------------------------------------------------------
% ihbdfl_ax_info
% ihbdfl_fr_info
% ihbdfl_sl_info
% ihbdfl_sr_info
%----------------------------------------------------------------------------------
% Axis limits to draw slice
%----------------------------------------------------------------------------------
% ihbdfl_xLimD
% ihbdfl_yLimD
% ihbdfl_zLimD
% ihbdfl_xLimI
% ihbdfl_yLimI
% ihbdfl_zLimI
%----------------------------------------------------------------------------------
% 3D view Renderer ('zbuffer' or 'OpenGL')
% 3D view transparency (available only for OpenGL
%----------------------------------------------------------------------------------
% ihbdfl_renderer
% ihbdfl_transparency
%----------------------------------------------------------------------------------
% Axes for 3D view type. If ihbdfl_bAxesIsClusterOnly == 1 - axes have cluster range
%                                                        0 - axes have brain range
%----------------------------------------------------------------------------------
% ihbdfl_bAxesIsClusterOnly
%----------------------------------------------------------------------------------
% Colors defaults
%----------------------------------------------------------------------------------
% ihbdfl_color_clOrig
% ihbdfl_color_clSymm
% ihbdfl_color_surface
% ihbdfl_color_cntAx
% ihbdfl_color_cntFr
% ihbdfl_color_cntSag
%----------------------------------------------------------------------------------
% Line width
%----------------------------------------------------------------------------------
% ihbdfl_line_width
%==================================================================================
%
% ..
%   09.04.01    Sergey Pakhomov
%   01.08.01    last modified
% ..

% ..
%    Create Default information file if necessary
% ..
dflFileName = ihb_FileFolderName('dfl');
fid = fopen(dflFileName);
if fid == -1
    ihb_ResetDefaults;
else
    fclose(fid);
end
%==================================================================================
% Delete previous figures
%==================================================================================
ihb_ResetFigure('ihb_TalSpaceMain_fig');
ihb_ResetFigure('ihb_TalSpace3D_fig');
ihb_ResetFigure('ihb_TalSpaceLeg_fig');
%==================================================================================
% Create main figure
%==================================================================================
screenUnits = get(0, 'units');      % Save current screen units
set(0, 'units', 'pixels');          % Set pixels units for screen
scrSize = get(0,'ScreenSize');      % Get screen size in pixels
set(0, 'units', screenUnits);       % Restore current screen units

load(dflFileName, 'ihbdfl_main_width');
load(dflFileName, 'ihbdfl_main_height');
pos(3) = ihbdfl_main_width;  
pos(4) = ihbdfl_main_height;  
pos(1) = (scrSize(3) - pos(3))/2;
pos(2) = (scrSize(4) - pos(4))/2;

hFigMain = figure(...
    'Visible', 'off',...
    'NumberTitle', 'off',...
    'MenuBar', 'none',...
    'Name', 'Talairach slices',...
    'Tag', 'ihb_TalSpaceMain_fig',...
    'CloseRequestFcn', 'ihb_CloseTalSpaceCB',...
    'ResizeFcn', 'ihb_ResizeTalSpaceCB',...
    'units', 'pixels',...
    'Position', pos);
%==================================================================================
% Set menu
%==================================================================================
%----------------------------------------------------------------------------------
% File menu
%----------------------------------------------------------------------------------
hFileMenu = uimenu('Label', 'File');
uimenu(hFileMenu, 'Label', 'New clusters', 'Callback', 'ihb_GetClusterSet');
uimenu(hFileMenu, 'Label', 'Exit', 'Callback', 'close', 'Separator', 'on'); 
%----------------------------------------------------------------------------------
% 3D View menu
%----------------------------------------------------------------------------------
h3DViewMenu = uimenu('Label', '3D View', 'Callback', 'ihb_TalSpaceView3D');
%----------------------------------------------------------------------------------
% Color Legend menu
%----------------------------------------------------------------------------------
hLegendMenu = uimenu('Label', 'Colors', 'Callback', 'ihb_ColorsTalSpace');
%----------------------------------------------------------------------------------
% Option menu
%----------------------------------------------------------------------------------
hOptionMenu = uimenu('Label', 'Options');
	%------------------------------------------------------------------------------
	%     Voxel size submenu
	%------------------------------------------------------------------------------
    hVoxSizeMenu = uimenu(hOptionMenu, 'Label', 'Voxel size');
        uimenu(hVoxSizeMenu, 'Tag', 'ihb_VoxSizeMenu111', 'Label', '1x1x1', 'Callback', 'ihb_TalSpaceVoxSizeMenuCB(1)');
        uimenu(hVoxSizeMenu, 'Tag', 'ihb_VoxSizeMenu222', 'Label', '2x2x2', 'Callback', 'ihb_TalSpaceVoxSizeMenuCB(2)');
        uimenu(hVoxSizeMenu, 'Tag', 'ihb_VoxSizeMenu444', 'Label', '4x4x4', 'Callback', 'ihb_TalSpaceVoxSizeMenuCB(4)');
	%------------------------------------------------------------------------------
	%     Renderer submenu
	%------------------------------------------------------------------------------
    hRendererMenu  = uimenu(hOptionMenu, 'Label', '3D renderer');
        uimenu(hRendererMenu, 'Tag', 'ihb_RendererMenuZbuffer', 'Label', 'Zbuffer', 'Callback', 'ihb_TalSpaceRendererMenuCB(1)');
        uimenu(hRendererMenu, 'Tag', 'ihb_RendererMenuOpenGL', 'Label', 'OpenGL', 'Callback', 'ihb_TalSpaceRendererMenuCB(2)');
	%------------------------------------------------------------------------------
	%     Transparency submenu
	%------------------------------------------------------------------------------
    hTransparencyMenu  = uimenu(hOptionMenu, 'Tag', 'ihb_TransparencyMenu', 'Label', '3D Transparency');
        uimenu(hTransparencyMenu, 'Tag', 'ihb_TransparencyMenuOpaque', 'Label', 'Opaque', 'Callback', 'ihb_TalSpaceTransparencyMenuCB(1.0)');
        uimenu(hTransparencyMenu, 'Tag', 'ihb_TransparencyMenuLow', 'Label', 'Low', 'Callback', 'ihb_TalSpaceTransparencyMenuCB(0.8)');
        uimenu(hTransparencyMenu, 'Tag', 'ihb_TransparencyMenuMedium', 'Label', 'Medium', 'Callback', 'ihb_TalSpaceTransparencyMenuCB(0.5)');
        uimenu(hTransparencyMenu, 'Tag', 'ihb_TransparencyMenuHigh', 'Label', 'High', 'Callback', 'ihb_TalSpaceTransparencyMenuCB(0.2)');
	%------------------------------------------------------------------------------
	%     3D axes submenu
	%------------------------------------------------------------------------------
    hAxes3DMenu  = uimenu(hOptionMenu, 'Label', '3D view axes');
        uimenu(hAxes3DMenu, 'Tag', 'ihb_Axes3DMenuClusterOnly', 'Label', 'Cluster Only', 'Callback', 'ihb_TalSpaceAxes3DMenuCB(1)');
        uimenu(hAxes3DMenu, 'Tag', 'ihb_Axes3DMenuWholeBrain', 'Label', 'Whole Brain', 'Callback', 'ihb_TalSpaceAxes3DMenuCB(0)');
	%------------------------------------------------------------------------------
	%     Smoothing submenu
	%------------------------------------------------------------------------------
    % hSmoothMenu = uimenu(hOptionMenu, 'Label', 'Smoothing');
	%------------------------------------------------------------------------------
	%     Line width submenu
	%------------------------------------------------------------------------------
    hLineWidthMenu = uimenu(hOptionMenu, 'Label', 'Line width');
        uimenu(hLineWidthMenu, 'Tag', 'ihb_LineWidthMenuThin', 'Label', 'Thin', 'Callback', 'ihb_TalSpaceLineWidthMenuCB(1)');
        uimenu(hLineWidthMenu, 'Tag', 'ihb_LineWidthMenuNormal', 'Label', 'Normal', 'Callback', 'ihb_TalSpaceLineWidthMenuCB(2)');
        uimenu(hLineWidthMenu, 'Tag', 'ihb_LineWidthMenuThick', 'Label', 'Thick', 'Callback', 'ihb_TalSpaceLineWidthMenuCB(4)');
    hResetMenu = uimenu(hOptionMenu, 'Label', 'Reset Defaults', 'Callback', 'ihb_TalSpaceResetMenuCB');
%----------------------------------------------------------------------------------
% Help menu
%----------------------------------------------------------------------------------
hHelpMenu  = uimenu('Label', 'Help');
    uimenu(hHelpMenu, 'Label', 'Web site', 'Callback', 'ihb_TalSpaceHelp');
    uimenu(hHelpMenu, 'Label', 'About', 'Callback', 'ihb_AboutTalSpace');
%----------------------------------------------------------------------------------
ihb_TalSpaceSetMenuCheck;
%----------------------------------------------------------------------------------
% Set axis and information UI
%----------------------------------------------------------------------------------
ihb_SetAxis(hFigMain);
ihb_SetInfoUI(hFigMain);
%----------------------------------------------------------------------------------
% Set mouse move callback
%----------------------------------------------------------------------------------
set(hFigMain, 'WindowButtonMotionFcn', 'ihb_MotionTalSpaceCB');   
%----------------------------------------------------------------------------------
% Select cluster set and cluster to view    while
%----------------------------------------------------------------------------------
bIsEmpty = 0;
while bIsEmpty == 0
    bIsEmpty = tor_ihb_GetClusterSet;
    if bIsEmpty == 0
        h = msgbox('Cluster to view should be selected', 'Warning!', 'warn', 'modal');
        waitfor(h);
    end
end
