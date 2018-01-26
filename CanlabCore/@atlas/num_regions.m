function n_regions = num_regions(obj)
% Count number of regions in atlas object, even with incomplete data

n1 = size(obj.probability_maps, 2);
n2 = length(obj.labels);
n3 = length(unique(obj.dat(obj.dat ~= 0)));

n_regions = max([n1 n2 n3]);

end