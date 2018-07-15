function [n_regions, n_regions_with_data, missing_regions] = num_regions(obj)
% Count number of regions in atlas object, even with incomplete data

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