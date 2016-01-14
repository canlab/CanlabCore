function CLU = transform_coordinates(CLU,mat)
% Transforms XYZ voxel coordinates in CLU (clusters or CLU)
% to new voxel coordinates, given mm coordinates in CLU
% and a mat file describing the transformation, as in SPM99
%
% :Usage:
% ::
%
%     CLU = transform_coordinates(CLU,mat)
%
% This preserves the order of the voxels, but is slower
% and gives UNIQUE XYZ voxels given XYZmm.
% see mm2voxel.m
%
% ..
%    tor wager
% ..

V.M = mat;

for i = 1:length(CLU)
    
    CLU(i).XYZ = mm2voxel(CLU(i).XYZmm,V)';
    CLU(i).XYZmm = voxel2mm(CLU(i).XYZ,mat);
    
end

return

    
