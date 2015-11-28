function cl = cluster_interp(cl,varargin)
% Interpolates voxels and mm in a clusters structure (cl) to match the
% image dimensions and voxel sizes of a target mask image.
%
% :Usage:
% ::
%
%    function cl = cluster_interp(cl,maskimg,[keep sep clusters flag])
%
% :Example:
% ::
%
%    cl = cluster_interp(cl,maskimg,1);
%
%    %default mask 2 x 2 x 2, SPM2 default:
%    cl = cluster_interp(cl,[],1);
%

%maskimg = which('scalped_avg152T1_graymatter_smoothed.img');
maskimg = which('scalped_avg152T1.img');
if length(varargin) > 0, maskimg = varargin{1};,end
if isempty(maskimg),maskimg = which('scalped_avg152T1.img');, end


dosep = 0;if length(varargin) > 1, dosep = varargin{2};,end
%dotext = 1; if dosep, dotext = 0;,end

% recursive call if clusters are to be kept separate
if dosep, 
    fprintf(1,'\n---------------------------------------\n');
    fprintf(1,'cluster_interp.m: creating separate masks for each cluster to preserve separation\n');
    fprintf(1,'\n---------------------------------------\n');
    for i = 1:length(cl)
        fprintf(1,'\n---------------------------------------\n');
        fprintf(1,'Cluster %3.0f ',i);
        fprintf(1,'\n---------------------------------------\n');
        cltmp = cluster_interp(cl(i),maskimg);
        
        if i == 1, cl2 = cltmp;, else, cl2 = merge_clusters(cl2,cltmp);, end
    end
    fprintf(1,'\n---------------------------------------\n');
    cl = cl2;
    
    return
end

% make mask of clusters

Vcl.mat = cl(1).M;
Vcl.dim = [91 109 91 2];
[m,Vcl] = clusters2mask(cl,Vcl);



% reslice mask

[tmp,newname] = reslice_imgs(maskimg,'clustermask.img');
if newname(1) == filesep, newname = newname(2:end);,end

% re-threshold mask to avoid larger clusters
Vcl = spm_vol(newname); v = spm_read_vols(Vcl);
v(v < .4) = 0;
spm_write_vol(Vcl,v);

% get clusters back
cl2 = mask2clusters(newname);


% put coords, etc. back in original cluster structure
if length(cl) == length(cl2)
    
for i = 1:length(cl2)
    cl(i).numVox = cl2(i).numVox;
    cl(i).voxSize = cl2(i).voxSize;
    cl(i).M = cl2(i).M;
    cl(i).Z = cl2(i).Z;
    cl(i).XYZmm = cl2(i).XYZmm;    
    cl(i).XYZ = cl2(i).XYZ;
    cl(i).center = cl2(i).center;
    cl(i).mm_center = cl2(i).mm_center;
end

else
    disp('Number of clusters has changed; cannot save old info.');
    cl = cl2;
end


%disp('Removing clustermask.img')
!rm clustermask.img
!rm clustermask.hdr

return
