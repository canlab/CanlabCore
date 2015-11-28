function mask = convert2mask(stats_image_obj)
% Converts each image in a statistic_image object into a mask object, based
% on significant voxels in the .sig field.
%
%
% :Example:
% ::
%
%    cl = region(convert2mask(timg), group)

mask = fmri_mask_image(stats_image_obj); % Copy into a mask image

% set dat to sig, and copy over other key fields
mask.dat = double(stats_image_obj.sig);
mask.removed_images = stats_image_obj.removed_images;
mask.removed_voxels = stats_image_obj.removed_voxels;

end


%{  
PREVIOUS CODE, NOT WORKING -- YONI REPLACED W/ ABOVE ON 8/5/2015

for i = 1:size(stats_image_obj.dat, 2)
    % for each image
    
    sigi = stats_image_obj.sig(:, i);
    
    if isempty(sigi)
        stats_image_obj.sig(:, i) = true(size(stats_image_obj.dat(:, i)));
    end
    
    % Copy into a mask image
    mask = fmri_mask_image;
    
    mask.volInfo = stats_image_obj.volInfo;
    
    mask.volInfo.image_indx(mask.volInfo.wh_inmask(~sigi)) = 0;
    mask.volInfo.wh_inmask(~sigi) = [];
    
    mask.volInfo.xyzlist(~sigi, :) = [];
    mask.volInfo.cluster(~sigi) = [];
    
    mask.volInfo.n_inmask = length(mask.volInfo.wh_inmask);
    
    mask.dat = zeros(mask.volInfo.nvox, 1, 'single');
    
    mask.dat(stats_image_obj.volInfo.wh_inmask(sigi)) = 1;

    if ~isempty(stats_image_obj.image_names)
        mask.image_names = stats_image_obj.image_names(i, :);
    end
    
    mask.history = stats_image_obj.history;
    
    mask.history{end + 1} = 'Converted to mask_image object.';
    
    varargout{i} = mask;
    
end

end % function
%}
