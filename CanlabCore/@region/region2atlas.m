function atlas_obj = region2atlas(r, reference_image_name)
% region2atlas Transform a region object into an atlas object.
%
% Build an atlas-class object from a region object, optionally
% resampling to the space of a reference image. Each region becomes a
% labeled parcel; voxels are coded with an integer index into the
% region array (or with the region's own .dat values, if available).
% Region shorttitles populate atlas labels and titles populate label
% descriptions.
%
% :Usage:
% ::
%
%     atlas_obj = region2atlas(r, [reference_image_name])
%
% :Inputs:
%
%   **r:**
%        A region-class object array.
%
% :Optional Inputs:
%
%   **reference_image_name:**
%        Filename string or image_vector object specifying the space to
%        resample the resulting atlas to. If omitted, the atlas is
%        returned in the space of r.
%
% :Outputs:
%
%   **atlas_obj:**
%        An atlas-class object with .dat coded by region index, .labels
%        from r.shorttitle, and .label_descriptions from r.title.
%
% :Examples:
% ::
%
%     % see parcellate_pain_predictive_regions.m
%     lindquist2015_pos_region_obj = region2fmri_data(r, obj);
%     % back to regions to test:
%     rtest = region(lindquist2015_pos_region_obj, 'unique_mask_values');
%     orthviews(rtest, 'unique')
%
% :See also:
%   - region2fmri_data
%   - region2imagevec
%   - atlas2region
%
% ..
%    Tor Wager: edited July 16, 2018
% ..

has_reference_image = nargin > 1;

if isempty(r)
    disp('Region object empty!!')
    return
end

ivec = region2imagevec(r);  % whole mask, to include all regions 
ivec = replace_empty(ivec);
ivec.dat = zeros(size(ivec.dat));

for i = 1:length(r)
    
    ivec_tmp = region2imagevec(r(i)); % For some reason this operation appends data to ivec_tmp.
    ivec_tmp = replace_empty(ivec_tmp);
    ivec_tmp.dat = [];

    whvox = ivec.volInfo.image_indx(ivec.volInfo.wh_inmask) & ivec_tmp.volInfo.image_indx(ivec.volInfo.wh_inmask);
    
    %ivec.dat(whvox) = ivec.dat(whvox) + i * double(ivec_tmp.dat ~= 0); % code with region number
    
    if isempty(ivec_tmp.dat) % Currently this will always be true. A future update should make it so we can port data over from region to atlas.
        ivec.dat(whvox) = i; % code with region number
    else
        ivec.dat(whvox) = ivec_tmp.dat; 
    end
    
    
end

% Convert to minimal atlas object
atlas_obj = atlas;
atlas_obj.volInfo = ivec.volInfo;
atlas_obj.dat = int32(round(ivec.dat));
atlas_obj.image_names = [];
atlas_obj.fullpath = [];
atlas_obj.files_exist = false;
atlas_obj.history = {'Created from region object.'};

atlas_obj.probability_maps = [];
atlas_obj.labels = {r.shorttitle};
atlas_obj.label_descriptions = {r.title}';

if has_reference_image
    
    if ischar(reference_image_name)
        
        reference_image_name = check_valid_imagename(reference_image_name, 1);
        
        reference_obj = fmri_data(reference_image_name, 'noverbose');
        
    elseif isa(reference_image_name, 'image_vector')
        
        reference_obj = reference_image_name; % it's already an object
        
    end
    
    atlas_obj = resample_space(atlas_obj, reference_obj);
    atlas_obj.dat = int32(round(atlas_obj.dat));

end


end % main function

