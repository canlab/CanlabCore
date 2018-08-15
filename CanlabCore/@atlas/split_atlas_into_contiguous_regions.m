function atlas_out = split_atlas_into_contiguous_regions(atlas_obj)
% Divide regions with multiple contiguous blobs into separate labeled regions for each blob
%
% Take an atlas object whose labeled regions contain multiple contiguous blobs and divide each
% contiguous blob into a separate labeled region.
% - Note: This version eliminates probability maps - handling them not implemented yet
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

[n_regions, n_regions_with_data, missing_regions] = num_regions(atlas_obj);

atlas_out = atlas_obj;

atlas_out.labels = {};
atlas_out.label_descriptions = {};

atlas_out.probability_maps = []; % Eliminate probability maps - handling them not implemented yet

atlas_out.dat = zeros(size(atlas_out.dat)); % avoid 0 for resampling

fprintf('Splitting %d regions: 000', n_regions)

for i = 1:n_regions
    
    fprintf('%3.0f', i);
    
    subatlas = select_atlas_subset(atlas_obj, i);
    
    %     wh = logical(subatlas.dat);                      % which voxels in region
    %
    %     xyz = subatlas.volInfo.xyzlist(wh, :);           % x-coordinates in voxels
    %
    %     XYZmm = voxel2mm(xyz',subatlas.volInfo.mat);
    
    % parse into contiguous regions. split if there are multiple.
    r = atlas2region(subatlas);
    r = reparse_continguous(r);
    
    % clean up left-out voxels
    wh_omit = cat(1, r.numVox) < 3;
    r(wh_omit) = [];
    
    % newr is atlas object for this original region, with contig regions
    % divided
    
    new_subatlas = region2atlas(r, subatlas);
    
    new_subatlas.probability_maps = []; % Eliminate probability maps - handling them not implemented yet

    new_subatlas.labels = get_new_labels(r, subatlas.labels{1});
    
    atlas_out = merge_atlases(atlas_out, new_subatlas);
    
    
end % region

%atlas_out.dat(atlas_out.dat == 9999) = 0;

atlas_out.references = unique(atlas_out.references, 'rows');

disp('Done.');


end % function





function new_labels = get_new_labels(r, old_label)
% later: could merge all L, all R, all M

k = length(r);

for j = 1:k
    
    voxsign = sign(r(j).XYZmm(1, :)); % sign of each x coord. - = left
    modal_x = mode(voxsign);
    
    switch modal_x
        case -1
            labelstr = '_L';
        case 1
            labelstr = '_R';
        case 0
            labelstr = '_M';
    end
    
    % proportional asymmetry. 1 is vary lateralized. 0 is balanced
    % across L/R.  If symmetrical, then label _M
    prop_asym = abs(sum(voxsign == -1) - sum(voxsign == 1)) ./ length(voxsign);
    
    if prop_asym < .5, labelstr = '_M'; end
    
    new_labels{j} = [old_label labelstr];
    
    
end % j

% add left/right/mid

%     % split region into left, right, and midline within 10 mm
%     isleft = XYZmm(1, :) < 0;
%     isright = XYZmm(1, :) > 0;   % what to do with midline?
%     ismid = XYZmm == 0;


end % subfunction





