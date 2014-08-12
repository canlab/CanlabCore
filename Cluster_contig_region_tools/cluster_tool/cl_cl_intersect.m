function [clout]=cl_cl_intersect(clin1,clin2,varargin)

% Usage:
% clout=cl_cl_intersect(clin1,clin2)
% clout=cl_cl_intersect(clin1,clin2,'switch')
% 
% clout is a cl structure containing clusters of only those voxels that
% exist in both cl structures clin1 and clin2.
% 
% If the clin1 and clin2 are not in the same image space, clin2 will be
% converted to the space of clin1 by finding the nearest mm point in the
% space of clin1 to each mm point in clin2. NOTE: ***** This may be
% undiserable if the voxel sizes in clin1 are substantially smaller than
% those in clin2--eg., a 4x4x4 voxel will be translated into a single 2x2x2
% voxel, resulting in 3 'empty' 2x2x2 voxels in the space covered by the
% 4x4x4 voxel. For this reason, it is recommended that you place the lowest
% resolution cl structure in the first input.
% 
% The cl.Z values and the threshold field are taken from clin1, and the
% values from clin2 are dropped. If the string 'switch' is included as in
% input, the Z values from clin2 will be used instead.
% 

sw=0;
if ~isempty(varargin)
    for k=1:length(varargin)
        if strcmp(varargin{k},'switch')
            sw=1;
        end
    end
end

XYZmm1=[];
XYZmm2=[];
Z=[];
for k=1:length(clin1)
    if size(clin1(k).XYZmm,1)==3
        XYZmm1(:,end+1:end+clin1(k).numVox)=clin1(k).XYZmm;
    else
        XYZmm1(:,end+1:end+clin1(k).numVox)=clin1(k).XYZmm';
    end
    if ~sw
        Z(end+1:end+clin1(k).numVox)=clin1(k).Z;
    end
end

for k=1:length(clin2)
    if size(clin2(k).XYZmm,1)==3
        XYZmm2(:,end+1:end+clin2(k).numVox)=clin2(k).XYZmm;
    else
        XYZmm2(:,end+1:end+clin2(k).numVox)=clin2(k).XYZmm';
    end
    if sw
        Z(end+1:end+clin2(k).numVox)=clin2(k).Z;
    end
end

if clin1(1).voxSize(:)~=clin2(1).voxSize(:)
    [XYZmm2,p_ind]=mmToVoxel(XYZmm2,clin1(1).M,'valid');
    if sw
        for m=1:length(p_ind)
            ind=find(p_ind(:)==p_ind(m));
            Z(ind)=mean(Z(ind));
        end
        Z=Z(unique(p_ind));
    end
    XYZmm2=voxelToMm(XYZmm2,clin1(1).M);
end

XYZmm_out=[];
Zout=[];

for k=1:size(XYZmm1,2)
    a=find(XYZmm2(1,:)==XYZmm1(1,k));
    if ~isempty(a)
        b=find(XYZmm2(2,a)==XYZmm1(2,k));
        if ~isempty(b)
            c=find(XYZmm2(3,a(b))==XYZmm1(3,k));
            if ~isempty(c)
                XYZmm_out(:,end+1)=XYZmm2(:,a(b(c)));
                Zout(:,end+1)=Z(a(b(c)));
            end
        end
    end
end
if isempty(XYZmm_out)
	clout=struct('XYZ',{},'XYZmm',{},'name',{},'numVox',{},'M',{},'voxSize',{},'dim',{},'title',{},'threshold',{},'Z',{},'mm_center',{});
else
    XYZout=mmToVoxel(XYZmm_out,clin1(1).M,'valid');
    cl_ind=spm_clusters(XYZout);
    for k=1:max(cl_ind)
        clout(k).XYZ=XYZout(:,cl_ind==k);
        clout(k).XYZmm=voxelToMm(clout(k).XYZ,clin1(1).M);
        clout(k).name='';
        clout(k).numVox=length(find(cl_ind==k));
        clout(k).M=clin1(1).M;
        clout(k).voxSize=clin1(1).voxSize;
        clout(k).dim=clin1(1).dim;
        clout(k).title=['Cluster of ' num2str(clout(k).numVox) ' voxels from the intersection of two clusters, created by cl_cl_intersect.m.'];
        clout(k).threshold=clin1(1).threshold;
        clout(k).Z=Zout(cl_ind==k);
        clout(k).mm_center=center_of_mass(clout(k).XYZmm,clout(k).Z);
    end
end
        
