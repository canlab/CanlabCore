function [data, XYZvoxSphere, XYZmmSphere] = iimg_sphere_timeseries(images, XYZmmCenter, radius)
% :Usage:
% ::
%
%     function [data, XYZvoxSphere, XYZmmSphere] = iimg_sphere_timeseries(images, XYZmm, radius)
%
% :Inputs:
%
%   **images:**
%        list of image files
%
%   **XYZmm:**
%        [3 x n] array of mm coords
%
%   **radiu:**
%        radius in mm of sphere to generate
%
% :Outputs:
%
%   **data:**
%        voxel data
%
%   **XYZvoxSphere:**
%        voxel

    Vimages = spm_vol(images);
    [XYZvox(1,:) XYZvox(2,:) XYZvox(3,:)] = ind2sub(Vimages(1).dim(1:3), 1:prod(Vimages(1).dim(1:3)));
    XYZmm = Vimages(1).mat(1:3, :)*[XYZvox; ones(1, size(XYZvox, 2))];
    
    dist = [XYZmm(1,:) - XYZmmCenter(1); XYZmm(2,:) - XYZmmCenter(2); XYZmm(3,:) - XYZmmCenter(3)];
    whVoxelsInSphere = find(sum(dist.^2) <= radius^2);
    
    XYZmmSphere = unique(XYZmm(:,whVoxelsInSphere)', 'rows')';
    XYZvoxSphere = unique(XYZvox(:,whVoxelsInSphere)', 'rows')';
    
    data = spm_get_data(Vimages, XYZvoxSphere);
end
