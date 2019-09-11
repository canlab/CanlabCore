function [values, components, full_data_objects, l2norms] = extract_gray_white_csf(obj,varargin)
% Extracts mean values (values) and top 5 component scores (components)
% from each of gray, white, and CSF masks.
%
% - Images must be in standard MNI space for this to apply.
% - Uses canonical masks in CANlab tools:'gray_matter_mask.img' 'canonical_white_matter.img' 'canonical_ventricles.img' 
%
% This currently uses the images
%       'gray_matter_mask.img' 'canonical_white_matter.img' 'canonical_ventricles.img'
%       These images are based on the SPM8 a priori tissue probability
%       maps, but they have been cleaned up and made symmetrical and/or eroded
%       so that the white and CSF compartments are unlikely to contain very
%       much gray matter.  The gray compartment is currently more
%       inclusive. The potential value of this is that signal in the CSF/white
%       compartments may be removed from images prior to/during analysis
%
% :Usage:
% ::
%
%     [values, components] = extract_gray_white_csf(obj,options)
%
% :Inputs:
%
%   **obj:**
%        an image_vector (e.g., fmri_data) object
%
% :Options:
%
%   **'eval':**
%       A function handle to use for computing summary statistics of each
%       tissue class. Must accept exactly 1 argument and handle 'nan' 
%       gracefully. e.g. '@(x1)(nanvar(x1))'
%
% :Outputs:
% 
%   **values:**
%        mean gray matter, white, CSF. If 'eval' option is passed,
%        specified function will be used in place of mean.
%
%   **components:**
%        first 5 components from each tissue class, observation x 5
%
%   **full_data_objects:**
%        Masked data objects for {gray white CSF}
%
%   **l2norms:**
%        Length-adjusted L2 norms (divided by sqrt(nvox))
%
%   **masks:**
%        optional input for using different masks; it should provide gray,
%        white and ventricle images in order. 
% ..
%    Tor Wager, July 21, 2015
% ..

% Programmers' notes:
% Jan 2017:  Issue with vector lengths if obj has removed images, fixed (Tor)
%
% Feb 2017: change to mask. Gray matter is now sparse, old gray_matter_mask.img thresholded at .5.
%
% Sept 2017: added option for custom function evaluation in place of mean.

fxn = @(x1)(nanmean(x1,1));
masks = {'gray_matter_mask_sparse.img' 'canonical_white_matter.img' 'canonical_ventricles.img'};

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'eval'
                fxn = varargin{i+1};
            case 'masks'
                masks = varargin{i+1};
        end
    end
end

numcomps = 5;

obj = remove_empty(obj);  % return only non-empty values
nimgs = size(obj.dat, 2); % - sum(obj.removed_images);

values = NaN * zeros(nimgs, length(masks));
components = cell(1, length(masks));
full_data_objects = [];

for i = 1:length(masks)
    
    maskname = which(masks{i});
    
    if isempty(maskname)
       fprintf('Image %s cannot be found on path.\n', masks{i});
       error('Exiting');
    else
        fprintf('Extracting from %s.\n', masks{i});
    end
    
    masked_obj = apply_mask(obj, maskname);

    % get all values, if requested
    if nargout > 2
        
        full_data_objects{i} = remove_empty(masked_obj);

    end
           
    % get means
    masked_obj.dat(masked_obj.dat == 0) = NaN;
    
    myvalues = fxn(masked_obj.dat);
    myvalues = myvalues(:);
    
    % may need to insert omitted - no, return in reduced space
%     if length(masked_obj.removed_images) > 1
%         myvalues = naninsert(masked_obj.removed_images, myvalues);
%     end
    
    values(:, i) = myvalues;
    

    % get components
    if nargout > 1
        
        % NaNs will mess this up - remove voxel-wise
        [wasnan, dataforpca] = nanremove(masked_obj.dat);
        if any(wasnan), fprintf('Removing %3.0f voxels with one or more NaNs\n', sum(wasnan)); end
        
        [~, components{i}] = pca(dataforpca', 'Economy', true, 'NumComponents', numcomps);
        
        % may need to insert omitted
%         if length(masked_obj.removed_images) > 1
%             components{i} = naninsert(masked_obj.removed_images, components{i});
%         end
        
    end
    
    if nargout > 3

        l2norms(:, i) = getnorms(masked_obj);

    end
    
end % end masks


end % function



function n = getnorms(masked_obj)

% Vector L2 norm / sqrt(length) of vector

% divide by this value to normalize image

normfun = @(x) sum(x .^ 2) .^ .5;

%masked_obj = remove_empty(masked_obj);

x = masked_obj.dat;

%nv = size(x, 1); 

for i = 1:size(x, 2)

    % remove nans, 0s
    xx = x(:, i);
    xx(xx == 0) = [];
    xx(isnan(xx)) = [];

    % divide by sqrt(length) so number of elements will not change scaling

    n(i) = normfun(xx) ./ sqrt(length(xx)); 


end


end % function


