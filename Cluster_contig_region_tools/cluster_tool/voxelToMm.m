function XYZmm = voxelToMm(XYZ,M)
% Usage:
% XYZmm = voxelToMm(XYZ,M)
% 
% XYZ is 3 vector point list, with XYZ assumed to be in ROWS if you input a
% list with 3 points.
% m is SPM mat - 4 x 4 affine transform (usually .M or .mat in a structure)
% (what's stored in the .mat file)
%
% essentially a copy of voxel2mm, exists for symmetry with mmToVoxel.
%
% Example:
% XYZmm = voxelToMm([x y z],volInfo.mat);

flip=0;
if isempty(XYZ), XYZmm = [];, return, end
if size(XYZ,1)~=3
    if size(XYZ,2)~=3,error('XYZ matrix must have 3 elements in one of its dimensions!')
    else XYZ=XYZ';flip=1;
    end
end

XYZ(4,:) = 1;	
XYZmm = M*XYZ;
XYZmm = XYZmm(1:3,:);
if flip,XYZmm=XYZmm';end