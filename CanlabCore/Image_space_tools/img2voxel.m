function [XYZ,XYZmm,val,V] = img2voxel(P,varargin)
%
% given a mask or filtered image file name, 
% returns XYZ coordinates in voxels and mm
% of nonzero, non-NaN voxels
%
% and img values at these coordinates in val
%
% Tor Wager 02/04/02
% Modified 2/15/06 to take raw data as well as image name
% 

XYZ = [];
XYZmm = [];
val = [];
V = [];

% needed for raw data input only
if length(varargin) > 0, vmat = varargin{1}; V.mat = vmat;, end

if isstr(P)
% it's an image name

    V = spm_vol(P);

    if length(V) > 1, error('Input only 1 image at once!'), end

    try
        vol = spm_read_vols(V);
    catch
        disp(['filename is: ' V.fname])
        disp('Cannot read!')
        return
    end

else
    % it's loaded image data
    vol = P;
end


% mask out zeros
% -------------------------------------------------------------------
vol = double(vol);
vol(vol == 0) = NaN;


% get XYZ
% -------------------------------------------------------------------
indx = find(~isnan(vol));
[x,y,z] = ind2sub(size(vol),indx);

XYZ = [x y z]';
XYZ = XYZ(1:3,:);


% get XYZmm
% (add one to multiply by the constant shift (offset from edge) in mat)
% -------------------------------------------------------------------
XYZ(4,:) = 1;	

XYZmm = V.mat * XYZ;

XYZmm = XYZmm(1:3,:);


% get val
% -------------------------------------------------------------------
val = vol(indx);

return