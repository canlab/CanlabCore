function atlas_out = split_atlas_by_hemisphere(atlas_obj)
% Divide regions that are bilateral into separate left- and right-hemisphere regions
%
% Take an atlas object whose labeled regions are bilateral and divide each
% region into separate left- and right-hemisphere regions.
% - Uses x-coordinates. x-coords exactly zero will be included in R hem regions
%
% July 2018, Tor Wager
%
% "split_atlas" method use cases:
% split_atlas_by_hemisphere: We have a defined set of bilateral regions that 
% we want to "hard-split" into left and right. Multiple discontiguous regions 
% with the same label will be kept together. 
%
% split_atlas_into_contiguous_regions: We want to (1) keep contiguous blobs together 
% that may cross the midline, and (2) separate contiguous blobs with the
% same label into separate labeled regions.


n_regions = num_regions(atlas_obj);

atlas_out = atlas_obj;

atlas_out.labels = {};
atlas_out.label_descriptions = {};
atlas_out.probability_maps = []; 
atlas_out.dat = zeros(size(atlas_out.dat)); 

fprintf('Splitting %d regions: 000', n_regions)

for i = 1:n_regions
    
    fprintf('%3.0f', i);
    
    subatlas = select_atlas_subset(atlas_obj, i);   % a single label
    
    wh = logical(subatlas.dat);                      % which voxels in region
    
    xyz = subatlas.volInfo.xyzlist(wh, :);           % x-coordinates in voxels
    
    XYZmm = voxel2mm(xyz',subatlas.volInfo.mat);
    
    isleft = XYZmm(1, :) < 0;
    
    % Data field
    newdat = zeros(sum(wh), 1);
    newdat(isleft) = 1;
    newdat(~isleft) = 2;
    
    subatlas.dat(wh) = newdat;
    
    % Probability map field
    if ~isempty(subatlas.probability_maps)
        
        pmaps = repmat(subatlas.probability_maps(wh, :), 1, 2);  % separate into L and R. should always contain 1 map
        pmaps(~isleft, 1) = 0;
        pmaps(isleft, 2) = 0;
        
        subatlas.probability_maps = repmat(subatlas.probability_maps, 1, 2);
        subatlas.probability_maps(wh, :) = pmaps;
        
    end
    
    subatlas.labels = get_new_labels(subatlas.labels{1});
    
    subatlas.label_descriptions = repmat(subatlas.label_descriptions, 2, 1);
    
    atlas_out = merge_atlases(atlas_out, subatlas);
    
    
end % region

%atlas_out.dat(atlas_out.dat == 9999) = 0;

atlas_out.references = unique(atlas_out.references, 'rows');

disp('Done.');


end % function





function new_labels = get_new_labels(old_label)

labelstr = '_L';
new_labels{1} = [old_label labelstr];

labelstr = '_R';
new_labels{2} = [old_label labelstr];

end % subfunction





