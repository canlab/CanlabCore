function varargout = convert2mask(stats_image_obj)
% Converts each image in a statistic_image object into a mask object, based
% on significant voxels in the .sig field.
%
% [mask1, mask2, etc...] = convert2mask(stats_image_obj)
%
% Examples
% cl = region(convert2mask(timg), group)

% enforce logical
stats_image_obj.sig = logical(stats_image_obj.sig);

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