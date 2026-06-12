function [xyz_out, indx, min_distance] = select_coordinates_near_regions(region_obj, xyz, cutoff_in_mm)
% select_coordinates_near_regions Filter coordinates by minimum distance to a region object.
%
% Take an r x 3 matrix of MNI mm coordinates (e.g., a meta-analysis
% coordinate database) and return only those coordinates that lie within
% cutoff_in_mm of any voxel in the supplied region object.
%
% :Usage:
% ::
%
%     [xyz_out, indx, min_distance] = select_coordinates_near_regions(region_obj, xyz, cutoff_in_mm)
%
% :Inputs:
%
%   **region_obj:**
%        A region-class object array.
%
%   **xyz:**
%        r x 3 matrix of XYZ mm coordinates (e.g., DB.xyz from a
%        meta-analysis).
%
%   **cutoff_in_mm:**
%        Distance threshold (in mm). Coordinates within this distance
%        from any coordinate in the region object will be retained.
%
% :Outputs:
%
%   **xyz_out:**
%        Subset of xyz containing only coordinates within cutoff_in_mm
%        of the region object.
%
%   **indx:**
%        Logical r x 1 vector indicating which rows of xyz were kept.
%
%   **min_distance:**
%        r x 1 vector of minimum distances from each input coordinate to
%        the region object (in mm).
%
% :See also:
%   - region
%   - dist

xyz1 = cat(2, region_obj.XYZmm)'; % region object

% xyz1 : q x 3 matrix of XYZ mm coordinates
% xyz : 3 x r matrix of XYZ mm coordinates

D = dist(xyz1, xyz');

min_distance = min(D, [], 1)';

indx = min_distance < cutoff_in_mm;

xyz_out = xyz(indx, :);

end % function

