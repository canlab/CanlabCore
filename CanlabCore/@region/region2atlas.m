function atlas_obj = region2atlas(r, reference_image_name)
% Transform region object r into atlas object atlas_obj, given reference
% image name in the space to resample to
%
% objout = region2atlas(r, [reference_image_name to sample to])
%
% Examples:
% see parcellate_pain_predictive_regions.m
% lindquist2015_pos_region_obj = region2fmri_data(r, obj);
% back to regions to test:
% rtest = region(lindquist2015_pos_region_obj, 'unique_mask_values');
% orthviews(rtest, 'unique')

% Tor Wager: edited July 16, 2018

has_reference_image = nargin > 1;

if isempty(r)
    disp('Region object empty!!')
    return
end

ivec = region2imagevec(r);  % whole mask, to include all regions 
ivec = replace_empty(ivec);
ivec.dat = zeros(size(ivec.dat));

for i = 1:length(r)
    
    ivec_tmp = region2imagevec(r(i));
    
    %ivec_tmp = replace_empty(ivec_tmp);
    
    whvox = ivec.volInfo.image_indx(ivec.volInfo.wh_inmask) & ivec_tmp.volInfo.image_indx(ivec.volInfo.wh_inmask);
    
    %ivec.dat(whvox) = ivec.dat(whvox) + i * double(ivec_tmp.dat ~= 0); % code with region number
    
    ivec.dat(whvox) = i; % code with region number
    
    
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

