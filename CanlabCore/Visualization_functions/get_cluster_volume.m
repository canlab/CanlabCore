function clOut = get_cluster_volume(clIn)
% :Usage:
% ::
%
%    clOut = get_cluster_volume(clIn)
%
% Adapted from:
%
% FORMAT clusters = ihb_getClusters
%
% Get cluster information (use [SPM,VOL,xX,xCon,xSDM] = spm_getSPM;)
%
% ..
%    01.03.01    Sergey Pakhomov
%    01.08.01    last modified
%
%    Only change is to eliminate stuff in beginning - just update the cluster.
%    Tor Wager, 10/18/01
% ..

clOut = clIn;

%----------------------------------------------------------------------------------
% Define default voxel size for canonical image
% This is OK to leave at 1's, because the head plotting utility adjusts
% to millimeters.
%----------------------------------------------------------------------------------
ihbdfl_tal_x_vox = 1;
ihbdfl_tal_y_vox = 1;
ihbdfl_tal_z_vox = 1;



%----------------------------------------------------------------------------------
% Prepare MNI volume data xMNI, yMNI, zMNI, and volMNI
%----------------------------------------------------------------------------------
rMin = min(clOut.XYZmm');
rMax = max(clOut.XYZmm');
voxMNI = clOut.voxSize;

% make all voxMNI positive: Correct bug with negative x voxel size
voxMNI = abs(voxMNI);

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
    index = round(x*y*ind3 + x*ind1 + ind2 + 1);  % tor changed to avoid rounding error sometimes - 11/13        
    volMNI(index) = clOut.Z(curVox);
    
    % waitbar stuff too much of a pain.
    %nTotalProcessed = nTotalProcessed + 1;
    %waitbar(nTotalProcessed/nTotalVox);
end

sizeMNI = size(volMNI);
sizeMNIprod = prod(sizeMNI);
xMNIreshaped = reshape(xMNI, [1 sizeMNIprod]);
yMNIreshaped = reshape(yMNI, [1 sizeMNIprod]);
zMNIreshaped = reshape(zMNI, [1 sizeMNIprod]);
xyzReshaped = [xMNIreshaped; yMNIreshaped; zMNIreshaped]; 

% tor eliminated this transformation line to keep in MNI space
%xyzTal = ihb_Mni2Tal(xyzReshaped);

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
