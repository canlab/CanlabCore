function [xyz_out, indx, min_distance] = select_coordinates_near_regions(region_obj, xyz, cutoff_in_mm) 
% xyz_out = select_coordinates_near_regions(region_obj, xyz, cutoff_in_mm) 
%
% xyz = r x 2 matrix of XYZ mm coordinates, e.g., DB.xyz from meta-analysis
% cutoff_in_mm: xyz coordinates within this distance from any coordinate in the region object will be saved 

xyz1 = cat(2, region_obj.XYZmm)'; % region object

% xyz1 : q x 3 matrix of XYZ mm coordinates
% xyz : 3 x r matrix of XYZ mm coordinates

D = dist(xyz1, xyz');

min_distance = min(D, [], 1)';

indx = min_distance < cutoff_in_mm;

xyz_out = xyz(indx, :);

end % function

