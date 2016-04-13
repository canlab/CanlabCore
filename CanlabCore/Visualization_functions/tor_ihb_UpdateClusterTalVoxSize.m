function clOut = tor_ihb_UpdateClusterTalVoxSize(clIn, nTotalVox, nTotalProcessed)
% :Usage:
% ::
%
%    clOut = ihb_ChangeTalVoxSize(clIn)
%
% :Inputs:
%
%   **clIn:**
%        initial input cluster (structure see ihb_GetClusters)
%
%   **nTotalVox:**
%        total nmb of voxel in all clusters (used for waitbar)
%
%   **nTotalProcessed:**
%        total nmb of voxel in already processed clusters (used for waitbar)
%
% :Output:
%
%   **clOut:**
%        output cluster with new Talariach related fields
%
% It is assumed that waitbar already exists
%
% ..
%    24.04.01    Sergey Pakhomov
%
%    - Tor changed to leave in MNI space!!!
% ..

% ..
% Load current voxel size
% ..
dflFileName = ihb_FileFolderName('dfl');
load(dflFileName, 'ihbdfl_tal_x_vox', 'ihbdfl_tal_y_vox', 'ihbdfl_tal_z_vox');
%----------------------------------------------------------------------------------
% Copy input cluster data
%----------------------------------------------------------------------------------
clOut = ihb_DefineCluster;
clOut.isSpmCluster  = clIn.isSpmCluster;
clOut.title         = clIn.title;
clOut.hThreshold    = clIn.hThreshold;
clOut.voxSize       = clIn.voxSize;
clOut.name          = clIn.name;
clOut.numVox        = clIn.numVox;
clOut.Z             = clIn.Z;
clOut.XYZmm         = clIn.XYZmm;
clOut.pVoxelLev     = clIn.pVoxelLev;
clOut.pClustLev     = clIn.pClustLev;
%----------------------------------------------------------------------------------
% Prepare MNI volume data xMNI, yMNI, zMNI, and volMNI
%----------------------------------------------------------------------------------
rMin = min(clOut.XYZmm');
rMax = max(clOut.XYZmm');
voxMNI = clOut.voxSize;
[xMNI, yMNI, zMNI] = meshgrid(...
    rMin(1) - voxMNI(1) : voxMNI(1) : rMax(1) + voxMNI(1),...
    rMin(2) - voxMNI(2) : voxMNI(2) : rMax(2) + voxMNI(2),...
    rMin(3) - voxMNI(3) : voxMNI(3) : rMax(3) + voxMNI(3));
    
volMNI = zeros(size(xMNI));

[x y z] = size(xMNI);

for curVox = 1:clOut.numVox
    tmp = clOut.XYZmm(:, curVox);       
    ind1 = (tmp(1) - rMin(1))/voxMNI(1) + 1;
    ind2 = (tmp(2) - rMin(2))/voxMNI(2) + 1;
    ind3 = (tmp(3) - rMin(3))/voxMNI(3) + 1;
    index = x*y*ind3 + x*ind1 + ind2 + 1;        
    volMNI(index) = clOut.Z(curVox);
    nTotalProcessed = nTotalProcessed + 1;
    waitbar(nTotalProcessed/nTotalVox);
end
%----------------------------------------------------------------------------------
% Transform to Talariach volume data
%----------------------------------------------------------------------------------
sizeMNI = size(volMNI);
sizeMNIprod = prod(sizeMNI);
xMNIreshaped = reshape(xMNI, [1 sizeMNIprod]);
yMNIreshaped = reshape(yMNI, [1 sizeMNIprod]);
zMNIreshaped = reshape(zMNI, [1 sizeMNIprod]);
xyzReshaped = [xMNIreshaped; yMNIreshaped; zMNIreshaped]; 

%xyzTal = ihb_Mni2Tal(xyzReshaped);
% tor eliminated this transformation line.
xyzTal = xyzReshaped;

xx = reshape(xyzTal(1, :), sizeMNI);
yy = reshape(xyzTal(2, :), sizeMNI);
zz = reshape(xyzTal(3, :), sizeMNI);

clOut.xMin = ihbdfl_tal_x_vox*floor(min(xyzTal(1, :))/ihbdfl_tal_x_vox);
clOut.yMin = ihbdfl_tal_y_vox*floor(min(xyzTal(2, :))/ihbdfl_tal_y_vox);
clOut.zMin = ihbdfl_tal_z_vox*floor(min(xyzTal(3, :))/ihbdfl_tal_z_vox);

clOut.xMax = ihbdfl_tal_x_vox*ceil(max(xyzTal(1, :))/ihbdfl_tal_x_vox);
clOut.yMax = ihbdfl_tal_y_vox*ceil(max(xyzTal(2, :))/ihbdfl_tal_y_vox);
clOut.zMax = ihbdfl_tal_z_vox*ceil(max(xyzTal(3, :))/ihbdfl_tal_z_vox);

[clOut.xTal, clOut.yTal, clOut.zTal] = meshgrid(clOut.xMin:ihbdfl_tal_x_vox:clOut.xMax,...
                                                clOut.yMin:ihbdfl_tal_y_vox:clOut.yMax,...
                                                clOut.zMin:ihbdfl_tal_z_vox:clOut.zMax);
volTal = interp3(xx, yy, zz, volMNI, clOut.xTal, clOut.yTal, clOut.zTal, 'linear');
index = find(isnan(volTal(:,:,:)));
volTal(index) = 0;
clOut.vTal = volTal;
