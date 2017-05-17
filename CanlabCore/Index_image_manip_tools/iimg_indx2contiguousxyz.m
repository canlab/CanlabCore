function cl = iimg_indx2contiguousxyz(dat,volInfo,remove_mean_flag)
% Take in index image data vector (in-mask values only) and a volume info structure with xyzlist
% and return a cl structure whose XYZ values list contiguous sets of voxels
% ("blobs")
%
% :Usage:
% ::
%
%     cl = iimg_indx2contiguousxyz(dat,volInfo,[remove_mean_flag])
%
% If a 3rd arg is entered, means of each blob are subtracted
% This is to facilitate randomizing blob centers in
% meta_stochastic_activation_blobs.m
%
% ..
%    Tor Wager, June 06
% ..

if nargin > 2, docenter = 1; else docenter = 0; end

%n = size(dat,1);

wh = find(dat);
xyz = volInfo.xyzlist(wh,:)';

cl = [];
if isempty(xyz), return, end

clusterid = spm_clusters(xyz);

for i = 1:max(clusterid)
    cl(i).XYZ = xyz(:,clusterid==i);
    
    if docenter
        m = mean(cl(i).XYZ,2);
        cl(i).XYZ(1,:) = cl(i).XYZ(1,:) - m(1);
        cl(i).XYZ(2,:) = cl(i).XYZ(2,:) - m(2);
        cl(i).XYZ(3,:) = cl(i).XYZ(3,:) - m(3);
        
        cl(i).XYZ = round(cl(i).XYZ);
    end
    
end

return
