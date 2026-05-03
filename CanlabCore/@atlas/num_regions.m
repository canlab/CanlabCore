function [n_regions, n_regions_with_data, missing_regions] = num_regions(obj)
% num_regions Count number of regions in atlas object, even with incomplete data.
%
% Returns the total number of regions defined by the atlas (taking the
% maximum across probability_maps columns, labels, and unique values in
% obj.dat), the number of regions that have at least one assigned voxel,
% and the index values of any regions with no assigned voxels.
%
% Missing regions can cause problems for some atlas functions; see
% atlas.check_properties.
%
% :Usage:
% ::
%
%     [n_regions, n_regions_with_data, missing_regions] = num_regions(obj)
%
% :Inputs:
%
%   **obj:**
%        An atlas-class object.
%
% :Outputs:
%
%   **n_regions:**
%        Total number of regions defined by the atlas.
%
%   **n_regions_with_data:**
%        Number of regions with at least one assigned voxel.
%
%   **missing_regions:**
%        Index values of regions with no assigned voxels.
%
% :See also:
%   - check_properties
%   - get_region_volumes

n1 = size(obj.probability_maps, 2);
n2 = length(obj.labels);

u = unique(obj.dat(obj.dat ~= 0));
n3 = length(u); % will be large if interpolation issue, needs to be integers

n_regions = max([n1 n2 n3]);

if nargout == 1
    return
end

% n_regions_with_data
n4 = 0;
if ~isempty(obj.probability_maps), n4 = sum(sum(obj.probability_maps) > 0); end

n_regions_with_data = max([n3 n4]);

if nargout == 2
    return
end

% which regions are missing if any
missing_regions = setxor(u, [1:n_regions]');

end % function