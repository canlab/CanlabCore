function [values, components] = extract_gray_white_csf(obj)
% [values, components] = extract_gray_white_csf(obj)
%
% Extracts mean values (values) and top 5 component scores (components)
% from each of gray, white, and CSF masks.
% Images must be in standard MNI space for this to apply.
%
% obj = an image_vector (e.g., fmri_data) object
%
% Tor Wager, July 21, 2015

numcomps = 5;

masks = {'gray_matter_mask.img' 'canonical_white_matter.img' 'canonical_ventricles.img'};

values = NaN * zeros(size(obj.dat, 2), length(masks));
components = cell(1, length(masks));

for i = 1:length(masks)
    
    maskname = which(masks{i});
    
    if isempty(maskname)
       fprintf('Image %s cannot be found on path.\n', masks{i});
       error('Exiting');
    else
        fprintf('Extracting from %s.\n', masks{i});
    end
    
    masked_obj = apply_mask(obj, maskname);

    
    % get means
    masked_obj.dat(masked_obj.dat == 0) = NaN;
    
    values(:, i) = nanmean(masked_obj.dat, 1)';
    
       
    % get components
    if nargout > 1
        
        % NaNs will mess this up - remove voxel-wise
        [wasnan, dataforpca] = nanremove(masked_obj.dat);
        if any(wasnan), fprintf('Removing %3.0f voxels with one or more NaNs\n', sum(wasnan)); end
        
        [~, components{i}] = pca(dataforpca', 'Economy', true, 'NumComponents', numcomps);
        
    end
    
    
    
end % end masks


end % function

