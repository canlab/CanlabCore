function XYZmm = voxel2mm(XYZ,m)
% :Usage:
% ::
%
%     function XYZmm = voxel2mm(XYZ,m)
%
% :Inputs:
%
%   **XYZ:**
%        is 3 vector point list (3 rows, n columns)
%
%   **m:**
%        is SPM mat - 4 x 4 affine transform
%        (what's stored in the .mat file)
%
% :Example:
% ::
%
%    XYZmm = voxel2mm([x y z]',V.mat);
%
% ..
%    Verified that this works 10/27/01.
%    Tor Wager, 10/27/01
%    3/8/26 - cosmetic changes to code/explanation only
% ..

if isempty(XYZ), XYZmm = []; return, end

% add one for the origin -- constant shift (offset from edge) in mat
% x_world = M * [x_i y_i z_i 1]';  
% -------------------------------------------------------------------
XYZ(4, :) = 1;	

XYZmm = m * XYZ;

XYZmm = XYZmm(1:3, :);

return
