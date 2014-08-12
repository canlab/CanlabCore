function [clusters,CLU] = sphere_mask(P,XYZmm,r,outname,varargin)
% [clusters,maskCLU] = sphere_mask(fname,XYZmm,r,outname,[maskname],[overlay])
%
% Creates mask images and clusters for a set of spheres defined around
% coordinates you specify.  Used for creating regions of interest (ROIs).
% Spheres may be masked with an anatomical mask file.
%
%
% P       is input image name with correct dimensions and vox sizes for your
% study
% XYZmm   is mm coordinates (row vector) for sphere center
% r       is radius in mm
% outname is string for output mask name, e.g., 'sphere_mask.img'
% [maskname]    is optional mask .img containing additional constraints
%       (e.g., gray matter mask, etc.), can be in different dimensions
% see also mask2clusters, montage_clusters
%
% tor wager, Fall 2003.
% modified March 2004 to provide additional anatomical masking based on
% anatomical ROI (e.g., a gray matter region from the ICBM template)
%
% Example:
% tmp = sphere_mask(EXPT.SNPM.P{1}(1,:),[25.6 -45.5 77.0],10,'test.img','ICBM_area74.img');
% icbm_localize(tmp)
%
% M =which('ICBM_brainonly_1mm_seg1.img')
% cl = sphere_mask(d(1).name,[35 -57 54; -23 -59 56;13 -63 62; -7 -77 50; 25 -79 30; -23 -81 18;41 -7 46;-37 -7 46; 35 -1 30; 29 -3 60; -29 -13 46; -5 -1 56;-37 29 30; 53 21 32;51 11 -4;43 -69 14;-45 -71 12;29 -83 4;11 -89 -4;-11 -79 2; -1 -95 -14; 17 -99 -8;31 -77 -20;-23 -81 -20;41 -71 -20;-43 -79 -12],8,'tmp.img',M);
%
% Matlab 6.5/OSX bug gives seg fault or something if mask is too big.

CLU = [];

% call this function recursively, saving mask CLU after first time
if size(XYZmm,1) > 1
    disp('Multiple ROI mode')
    cl = [];
    for i = 1:size(XYZmm,1)
        if length(varargin) > 0, [out,varargin{1}] = sphere_mask(P,XYZmm(i,:),r,outname,varargin{1});,
        else, [out] = sphere_mask(P,XYZmm(i,:),r,outname);
        end
        cl = [cl,out];
    end
    clusters = cl;
    
    if length(varargin) > 1, ovl = varargin{2};, else, ovl = [];, end
    %montage_clusters(which(P),clusters,{'r'});
    cluster_orthviews(clusters,{[1 0 0]},'overlay',ovl)
    
    V = spm_vol(P); V.fname = outname;
    [m] = clusters2mask(clusters,V.dim(1:3));
    spm_write_vol(V,m);
    
else
    
P = which(P);    
if isempty(P), P = which('scalped_avg152T1_graymatter_smoothed.img');,end

if isempty(P), disp(['Cannot find image with vox sizes for your study.']);,end

V = spm_vol(P);


if length(varargin) > 0
    m = varargin{1};
    if isstr(m), m = which(m);, end
    if isempty(m), error('Cannot find mask image.');, end
end

V.M = V.mat;
XYZ = mm2voxel(XYZmm,V);
rv = abs(r ./ diag(V.M(1:3,1:3))');   % voxel coords


% build box (list of XYZ voxel coords)

lim = round([XYZ - rv; XYZ + rv]);
diffs = diff(lim);

xtmp = prod([diffs(2)+1 diffs(3)+1]);
ztmp = prod([diffs(1)+1 diffs(2)+1]);

x = repmat((lim(1,1):lim(2,1))',xtmp,1);

y = []; for i=1:diffs(2)+1, 
    ytmp = repmat(lim(1,2)+i-1,diffs(1)+1,1); y = [y;ytmp];,
end
y = repmat(y,diffs(3)+1,1);

ztmp = repmat(lim(1,3),ztmp,1);
z = [];
for i = 1:diffs(3)+1
    z = [z; ztmp+i-1];
end

xyz2 = [x y z];

xyz2mm = voxel2mm(xyz2',V.mat)';
xyzmm = repmat(XYZmm,size(xyz2,1),1);

if isempty(xyzmm), error('No coordinates in mask.'), end
    
d = sum((xyz2mm - xyzmm) .^ 2,2) .^ .5;


if length(varargin) > 0
    % mask - determine in-mask voxels
    fprintf(1,'Applying mask... ')
    %m = varargin{1}; 
    if isstr(m), [dum,CLU] = mask2clusters(m);  % mask name,
    else
        CLU = m;    % input CLU structure
    end
    [t1,t2,t3]=intersect(round(xyz2mm),round(CLU.XYZmm'),'rows');
    t4 = zeros(size(xyz2mm,1),1); t4(t2) = 1;
else    
    % No mask, all voxels are OK
    t4 = ones(size(xyz2mm,1),1);
end

wh = find(d <= r & t4);



xyz2 = xyz2(wh,:);
xyz2mm = xyz2mm(wh,:);

clusters = xyz2clusters(xyz2mm,P);

% mask = voxel2mask(xyz2,V.dim(1:3));
% 
% if length(size(mask)) > 3
%     mask = mask(:,:,:,1);
% end
% 
% V.fname = outname;
% spm_write_vol(V,mask);
% fprintf(1,['Written ' V.fname])
% 
% clusters = mask2clusters(outname);
% if isempty(clusters)
%     fprintf(1,'Empty clusters!')
% else
%     %cluster_orthviews(clusters,{[1 0 0]})
%     %montage_clusters(P,clusters,{'r'});
% end

end


return
