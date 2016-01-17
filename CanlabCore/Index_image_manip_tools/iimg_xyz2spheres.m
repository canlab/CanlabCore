function indx = iimg_xyz2spheres(xyz, mask_xyzlist, r)
% :Usage:
% ::
%
%     indx = iimg_xyz2spheres(xyz, mask_xyzlist, r)
%
% :Inputs:
%
%   **mask_xyzlist:**
%        is a list of voxel coords of all in-mask voxels
%
%   **xyz:**
%        is a list of voxels to be convolved with spheres
%
%   **r:**
%        is sphere radius (voxels)
%
% more generally, finds mask_xyzlist entries that are within r units of xyz
% i.e., to find database points within r mm of a cluster xyz list:
% indx = iimg_xyz2spheres(clusterxyz,databasexyz,r)
%
% indx is length mask_xyzlist and contains all in-mask, in-sphere voxels
% for all xyz coordinates.
%
% ..
%    See programmers' notes testing this function: Tested OK, 7/28/13, tor
% ..

    num_mask_voxels = size(mask_xyzlist,1);
    indx = zeros(num_mask_voxels,1);
    
    if isempty(xyz), return; end

    wh_orig = (1:num_mask_voxels)';        % original voxel index


    % eliminate mask_xyzlist points that are not even close (outside box)
    cmax = max(xyz,[],1);
    cmin = min(xyz,[],1);
    d = [(mask_xyzlist - repmat(cmax,num_mask_voxels,1)) (repmat(cmin,num_mask_voxels,1) - mask_xyzlist)];
    wh = any(d > r,2);

    mask_xyzlist(wh,:) = [];
    wh_orig(wh) = [];

    r2 = r.^2;

    % for each voxel, find mask_xyzlist coords within r
    % add ones to indx
    for i = 1:size(xyz,1)
        sphere_center = xyz(i,:);
        xyz_candidates = mask_xyzlist;
        wh_candidates = wh_orig;

        % find list xyz in box around coord
        wh = abs(xyz_candidates(:,1) - sphere_center(1)) > r;
        xyz_candidates(wh,:) = [];
        wh_candidates(wh) = [];

        wh = abs(xyz_candidates(:,2) - sphere_center(2)) > r;
        xyz_candidates(wh,:) = [];
        wh_candidates(wh) = [];

        wh = abs(xyz_candidates(:,3) - sphere_center(3)) > r;
        xyz_candidates(wh,:) = [];
        wh_candidates(wh) = [];

        num_candidates = size(xyz_candidates,1);
        d = sum((repmat(sphere_center,num_candidates,1) - xyz_candidates ).^2, 2);

        wh = wh_candidates(d <= r2);

        indx(wh) = 1;
    end
end


% Programmers' notes: Test

% img = which('spm2_hipp.img')
% dat = fmri_data(img);
% dat = remove_empty(dat);
% 
% % Target
% myxyz1 = dat.volInfo.xyzlist(~dat.removed_voxels , :);
% XYZc = voxel2mm(myxyz1', dat.volInfo.mat)';
% 
% r = region(dat);
% cluster_orthviews(r, {[1 0 0]}, 'solid');
% 
% dat2 = preprocess(dat, 'smooth', 4);
% dat2 = threshold(dat2, [.1 Inf], 'raw-between');
% myxyz = dat.volInfo.xyzlist(~dat2.removed_voxels , :);
% myxyz = voxel2mm(myxyz', dat2.volInfo.mat)';
% 
% indx2 = iimg_xyz2spheres(XYZc, myxyz, mydist);
% dat2.dat(~indx2) = 0;
% dat2 = remove_empty(dat2);
% 
% r2 = region(dat2);
% cluster_orthviews(r2, {[0 1 0]}, 'add');
% 
% %%
% mydist = 2;
% 
% dat2 = preprocess(dat, 'smooth', 4);
% dat2 = threshold(dat2, [.1 Inf], 'raw-between');
% myxyz = dat.volInfo.xyzlist(~dat2.removed_voxels , :);
% myxyz = voxel2mm(myxyz', dat2.volInfo.mat)';
% 
% indx2 = iimg_xyz2spheres(XYZc, myxyz, mydist);
% dat2.dat(~indx2) = 0;
% dat2 = remove_empty(dat2);
% 
% r2 = region(dat2);
% cluster_orthviews(r, {[1 0 0]}, 'solid');
% cluster_orthviews(r2, {[0 0 1]}, 'add');

