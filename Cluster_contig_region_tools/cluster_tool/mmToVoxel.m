function varargout = mmToVoxel(XYZmm,M,varargin)
% Usage:
% XYZ = mmToVoxel(XYZmm,M)
% [XYZ,point_index] = mmToVoxel(XYZmm,M,'valid')
% 
% XYZ is 3 vector point list, with XYZ assumed to be in ROWS if you input a
% list with 3 points.
% m is SPM mat - 4 x 4 affine transform (usually .M or .mat in a structure)
% (what's stored in the .mat file)
% 
% use of the 'valid' text flag will force all output XYZ values to be real
% positive integers. This is done by rounding all values and eliminating
% all points with an X, Y, or Z value <1. Only unique points will be
% reported.
% 
% the point_index output variable is a list of indices that can be used
% in a corresponding matrix or vector (e.g. a vector of Z scores 
% corresponding to each voxel) for removing non-unique points.
% For example, the following code will take the mean of the Z scores of
% non-unique points assigned to the same coordinate by using the 'valid'
% switch, and ensure that each point in the Z vector still corresponds to
% the correct coordinate:
% [XYZ,point_index]=mmToVoxel(XYZmm,M,'valid');
% for m=1:max(point_index)
%     ind=find(point_index(:)==point_index(m));
%     Z(ind)=mean(Z(ind));
% end
% Z=Z(unique(point_index));
% 
% the above is primarily useful when converting a list of XYZmm points that
% did not originate in the same space defined by the affine matrix.
%
% intended to replace mm2voxel, which is not symmetric with voxel2mm and
% requires a structure with a .M or .mat field in it to be passed in,
% rather than looking for the affine matrix to be passed in directly.
% 
% Example:
% XYZ = mmToVoxel([x y z],volInfo.mat);

valid=0;
if ~isempty(varargin)
    for k=1:length(varargin)
        if strcmp(varargin{k},'valid')
            valid=1;
        end
    end
end

flip=0;
if isempty(XYZmm), XYZ = [];, return, end
if size(XYZmm,1)~=3
    if size(XYZmm,2)~=3,error('XYZ matrix must have 3 elements in one of its dimensions!')
    else XYZmm=XYZmm';flip=1;
    end
end

XYZmm(4,:) = 1;	
XYZ=(M^-1)*XYZmm;
XYZ=XYZ(1:3,:);
if valid
    XYZ=round(XYZ);
    [row,col]=find(XYZ<1);
    XYZ(:,col)=[];
    [l,m,n]=unique(XYZ','rows');
    XYZ=XYZ(:,unique(m));
    varargout{2}=m;
end
if flip,XYZ=XYZ';end
varargout{1}=XYZ;