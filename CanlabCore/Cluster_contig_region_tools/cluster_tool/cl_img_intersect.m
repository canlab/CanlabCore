function [clout]=cl_img_intersect(clin,img_path)

% Usage:
% clout=cl_img_intersect(clin1,img)
% 
% img may be a pathname or a volInfo structure. The function does not check
% the integrity of a volInfo structure but will attempt to operate on any
% structure that is passed to it.
% 
% clout is a cl structure containing clusters of only those voxels that
% exist in both cl structure clin and are in-mask in the image specified by
% img_path.
% 
% If clin and the image are not in the same image space, clin will be
% converted to the image space by finding the nearest mm point in the
% image space to each mm point in clin. NOTE: ***** This may be undiserable
% if the voxel sizes in the image are substantially smaller than those in
% clin2--eg., a 4x4x4 voxel will be translated into a single 2x2x2 voxel,
% resulting in 3 'empty' 2x2x2 voxels within the space covered by the 4x4x4
% voxel.
% 
% The cl.Z values and the threshold field are taken from clin1, and the
% values from the image are dropped. 
% 
% If multiple points in clin are mapped to the same point in the image, the
% mean of all Z values for points in clin mapped to that point in the image
% is used.
% 

if isstruct(img_path)
    volInfo=img_path;
    img_path=volInfo.fname;
else
    volInfo=iimg_read_img(img_path,1);
end

if size(volInfo.xyzlist,1)~=3
    volInfo.xyzlist=volInfo.xyzlist';
end

XYZ=[];
Z=[];
for k=1:length(clin)
    if size(clin(k).XYZmm,1)~=3
        clin(k).XYZmm=clin(k).XYZmm';
    end
    if size(clin(k).XYZ,1)~=3
        clin(k).XYZ=clin(k).XYZ';
    end
    if clin(k).voxSize(:)~=diag(volInfo.mat(1:3,1:3))
        clin(k).voxSize=diag(volInfo.mat(1:3,1:3));
        [clin(k).XYZ,p_ind]=mmToVoxel(clin(k).XYZmm,volInfo.mat,'valid');
        for m=1:length(p_ind)
            ind=find(p_ind(:)==p_ind(m));
            clin(k).Z(ind)=mean(clin(k).Z(ind));
        end
        clin(k).Z=clin(k).Z(unique(p_ind));
        clin(k).XYZmm=voxelToMm(clin(k).XYZ,volInfo.mat);
        clin(k).M=volInfo.mat;
        clin(k).dim=volInfo.dim;
        clin(k).numVox=size(clin(k).XYZmm,2);
    end
    XYZ(:,end+1:end+clin(k).numVox)=clin(k).XYZ;
    Z(end+1:end+clin(k).numVox)=clin(k).Z;
end

XYZ_out=[];
Zout=[];

for k=1:size(XYZ,2)
    a=find(volInfo.xyzlist(1,:)==XYZ(1,k));
    if ~isempty(a)
        b=find(volInfo.xyzlist(2,a)==XYZ(2,k));
        if ~isempty(b)
            c=find(volInfo.xyzlist(3,a(b))==XYZ(3,k));
            if ~isempty(c)
                XYZ_out(:,end+1)=volInfo.xyzlist(:,a(b(c)));
                Zout(:,end+1)=Z(k);
            end
        end
    end
end
if isempty(XYZ_out);
    clout=struct('XYZ',{},'XYZmm',{},'name',{},'numVox',{},'M',{},'voxSize',{},'dim',{},'title',{},'threshold',{},'Z',{},'mm_center',{});
else
    cl_ind=spm_clusters(XYZ_out);
    for k=1:max(cl_ind)
        clout(k).XYZ=XYZ_out(:,cl_ind==k);
        clout(k).XYZmm=voxelToMm(clout(k).XYZ,clin(1).M);
        clout(k).name='';
        clout(k).numVox=length(find(cl_ind==k));
        clout(k).M=clin(1).M;
        clout(k).voxSize=clin(1).voxSize;
        clout(k).dim=clin(1).dim;
        clout(k).title=['Cluster of ' num2str(clout(k).numVox) ' voxels from the intersection of a cluster and ' img_path ', created by cl_img_intersect.m.'];
        clout(k).threshold=clin(1).threshold;
        clout(k).Z=Zout(cl_ind==k);
        clout(k).mm_center=center_of_mass(clout(k).XYZmm,clout(k).Z);
    end
end