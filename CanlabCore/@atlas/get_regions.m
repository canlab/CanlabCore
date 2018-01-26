function dat_out = get_regions(dat, varargin)
% Given an atlas object dat, select voxels corresponding to a subset of
% regions and return those in an output atlas object, dat_out
%
% dat_out = get_regions(dat, varargin)


% Find which names match
wh = ~cellfun(@isempty, strfind(names, string_to_find));

% Index values to look for
indx = find(wh)';

output_names = names(indx);


[subregions, voxel_matches] = select_voxels_by_value(dat, [5 7 9]);
