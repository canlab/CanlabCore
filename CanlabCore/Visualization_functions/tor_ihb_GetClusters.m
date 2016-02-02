function clusters = tor_ihb_GetClusters()
% :Usage:
% ::
%
%    clusters = ihb_getClusters
%
% Get cluster information (use [SPM,VOL,xX,xCon,xSDM] = spm_getSPM;)
%
% :Output:
%
%   **clusters:**
%        array of structs with fields:
%
% :Common to all clusters in the set:
%
%   **isSpmCluster:**
%        1 if got from SPM data; 0 - if constructed by 
%                     symmetrical or intersection
%
%   **title:**
%        title for comparison (string) (SPM.title)
%
%   **hThreshold:**
%        height threshold (SPM.u)
%
%   **voxSize:**
%        voxel dimensions {mm} - column vector (VOL.VOX)
%
% :Specific for every cluster:
%
%   **name:**
%        name of the cluster
%
%   **numVox:**
%        number of voxels in cluster
%
%   **Z:**
%        minimum of n Statistics {filtered on u and k} (1 x num_vox)
%
%   **XYZmm:**
%        location of voxels {mm} (3 x num_vox)
%
%   **pVoxelLev:**
%        corrected p for max value in cluster
%
%   **pClustLev:**
%        corrected p for given cluster (cluster level)
%
% :Talariach volume:
%
%   **xTal:**
%        x Talariach coordinates ready for contourslice & isosurface
%
%   **yTal:**
%        y Talariach coordinates ready for contourslice & isosurface
%
%   **zTal:**
%        z Talariach coordinates ready for contourslice & isosurface
%
%   **vTal:**
%        Talariach volume values ready for contourslice & isosurface
%
%   xMin, yMin, zMin, xMax, yMax, zMax - bounding box in mm for Talariach
%
% ..
%    10.04.01    Sergey Pakhomov
%    01.08.01    last modified
% ..

dflFileName = ihb_FileFolderName('dfl'); % Check existance of figures used by spm_getSPM
load(dflFileName, 'ihbdfl_spm_fig_Interactive');
load(dflFileName, 'ihbdfl_spm_fig_SelFileWin');
load(dflFileName, 'ihbdfl_spm_fig_ConMan');
hFigInteractive = findobj('Type', 'figure', 'Tag', ihbdfl_spm_fig_Interactive);
hFigSelFileWin = findobj('Type', 'figure', 'Tag', ihbdfl_spm_fig_SelFileWin);
hFigConMan = findobj('Type', 'figure', 'Tag', ihbdfl_spm_fig_ConMan);
%----------------------------------------------------------------------------------
% Get SPM information
%----------------------------------------------------------------------------------
[SPM,VOL,xX,xCon,xSDM] = spm_getSPM;
%----------------------------------------------------------------------------------
% Delete figures used by spm_getSPM if they were created during above call
%----------------------------------------------------------------------------------
if isempty(hFigInteractive) 
    delete(findobj('Type', 'figure', 'Tag', ihbdfl_spm_fig_Interactive)); 
end;
if isempty(hFigSelFileWin) 
    delete(findobj('Type', 'figure', 'Tag', ihbdfl_spm_fig_SelFileWin)); 
end;
if isempty(hFigConMan) 
    delete(findobj('Type', 'figure', 'Tag', ihbdfl_spm_fig_ConMan)); 
end;
%----------------------------------------------------------------------------------
if size(SPM.XYZ, 2) == 0
    clusters = [];
    h = warndlg('No voxels survive height threshold','Warning');
    delete(h);
    return;
end
A = spm_clusters(SPM.XYZ);
numClust = max(A);
clusters = [];
%----------------------------------------------------------------------------------
% Main cluster loop
%----------------------------------------------------------------------------------
load(dflFileName, 'ihbdfl_min_size_to_draw');
load(dflFileName, 'ihbdfl_tal_x_vox', 'ihbdfl_tal_y_vox', 'ihbdfl_tal_z_vox');
hWait = waitbar(0,'Preparing cluster data. Please wait...');
nTotalVox = size(A, 2);
nTotalProcessed = 0;

for curClust = 1:numClust
    %--------------------------------------------------------------------
    % Check current cluster size and skip if too small
    %--------------------------------------------------------------------
    a = find(A == curClust);
    if (size(a, 2) < ihbdfl_min_size_to_draw) | (size(a, 2) < SPM.k)
        nTotalProcessed = nTotalProcessed + size(a, 2);
        waitbar(nTotalProcessed/nTotalVox);
        continue; 
    end;
    cl = ihb_DefineCluster;
    %--------------------------------------------------------------------
    % General cluster data
    %--------------------------------------------------------------------
    cl.isSpmCluster = 1;
    cl.title = SPM.title;
    cl.hThreshold = SPM.u;
    cl.voxSize = VOL.VOX;
    %--------------------------------------------------------------------
    % Cluster specific data from SPM (MNI space)
    %--------------------------------------------------------------------
    cl.name = strcat(cl.title, '_', mat2str(size(a, 2)));
    cl.numVox = size(a, 2);   
    cl.Z = SPM.Z(a);
    cl.XYZmm = SPM.XYZmm(:,a);
    cl.pVoxelLev = spm_P(1, 0, max(cl.Z), SPM.df, SPM.STAT, VOL.R, SPM.n);
    cl.pClustLev = spm_P(1, cl.numVox/prod(VOL.FWHM), SPM.u, SPM.df, SPM.STAT, VOL.R, SPM.n);
    %====================================================================
    % Talariach volume related data
    %====================================================================
    clOut = tor_ihb_UpdateClusterTalVoxSize(cl, nTotalVox, nTotalProcessed); % tor changed only this line
    nTotalProcessed = nTotalProcessed + clOut.numVox;  
    clusters = [clusters, clOut];
end
close(hWait);
