function obj = resample_space(obj, sampleto, varargin)
% Resample the images in an fmri_data object (obj) to the space of another
% image (sampleto; e.g., a mask image). Works for all image_vector objects.
%
% obj = resample_space(obj, sampleto, [sampling method])
%
% Sampleto may be one of these:
% 1) a volInfo structure (the image does not have to exist on the path)
% 2) an image_vector, fmri_data, fmri_mask_image object
% 3) a string with the name of an image
%
% Can enter resampling method as optional input. Takes any input to
% interp3:
%       'nearest' - nearest neighbor interpolation
%       'linear'  - linear interpolation (default)
%       'spline'  - spline interpolation
%       'cubic'   - cubic interpolation as long as the data is uniformly
%                   spaced, otherwise the same as 'spline'
%
% Examples:
% label_mask = fmri_data(which('atlas_labels_combined.img'));
% label_mask = resample_space(label_mask, ivec, 'nearest') % resamples and masks label image

% programmers' notes:
% 1/27/2012 Tor edited to handle .mask field in fmri_data and .sig field in
% statistic_image.  Was causing errors otherwise...
%           Also changed automatic behavior to reparse contig voxels with
%           'nonempty' in output obj

n_imgs = size(obj.dat, 2);


Vto = sampleto.volInfo;
SPACEto = define_sampling_space(Vto, 1);

Vfrom = obj.volInfo;
SPACEfrom = define_sampling_space(Vfrom, 1);

obj_out = obj;
obj_out.dat = [];
obj_out.volInfo = Vto;

obj = replace_empty(obj);  % to make sure vox line up

for i = 1:n_imgs
    
    voldata = iimg_reconstruct_vols(obj.dat(:, i), obj.volInfo);
    
    resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, varargin{:});
    
    resampled_dat = resampled_dat(:);
    
    obj_out.dat(:, i) = resampled_dat(Vto.wh_inmask);
    
end

% in case of NaN values
obj_out.dat(isnan(obj_out.dat)) = 0;

if isa(obj_out, 'statistic_image')
   % Rebuild fields specific to statistic_images
   
   obj_out = replace_empty(obj_out);
   k = size(obj_out.dat, 2);
   
   for i = 1:k
       % this may break if nvox (total in image) is different for 2
       % images...
       p = ones(obj.volInfo.nvox, k);
       p(obj.volInfo.wh_inmask, i) = obj.p(:, i);
       
       ste = Inf .* ones(obj.volInfo.nvox, k);
       ste(obj.volInfo.wh_inmask, i) = obj.ste(:, i);
       
       sig = zeros(obj.volInfo.nvox, k);
       sig(obj.volInfo.wh_inmask, i) = obj.sig(:, i);
   
   end
   
   obj_out.p = p(Vto.wh_inmask, :);
   obj_out.ste = ste(Vto.wh_inmask, :);
   obj_out.sig = sig(Vto.wh_inmask, :);
   
end


if size(obj_out.dat, 1) == sum(obj_out.volInfo.image_indx)
    % this should always/almost always be true - assign missing/removed vox
    obj_out.removed_voxels = ~obj_out.volInfo.image_indx;
    obj_out.removed_voxels = obj_out.removed_voxels(obj_out.volInfo.wh_inmask);
    
else
    obj_out.removed_voxels = false;
end

% add clusters if needed
if obj_out.volInfo(1).n_inmask < 50000
    obj_out.volInfo(1).cluster = spm_clusters(obj_out.volInfo(1).xyzlist')';
else
    obj_out.volInfo(1).cluster = ones(obj_out.volInfo(1).n_inmask, 1);
end

% No longer need to remove - tor 5/27/15
% if isa(obj_out, 'statistic_image') && ~isempty(obj_out.sig)
%     disp('resample_space: removing threshold information from statistic_image')
%     obj_out.sig = [];
% end

obj = obj_out;

% This stuff below added 1/27/13 by tor

% re-parse clusters
obj = reparse_contiguous(obj, 'nonempty');

if isa(obj, 'fmri_data')
    % fmri_data has this field, but other image_vector objects do not.
    obj.mask = resample_space(obj.mask, sampleto);
end

% if isa(obj, 'statistic_image')
%     % statistic_image has this field, but other image_vector objects do not.
%     obj.sig = ones(size(obj.dat));
%     disp('.sig field reset. Re-threshold if necessary.');
% end

obj.history{end+1} = sprintf('Resampled data to space of %s', sampleto.volInfo.fname);


end

% %
% % % ---------------------------------------------------------
% % % Define image name to sample to and volInfo from input
% % % ---------------------------------------------------------
% % % Sampleto may be one of these:
% % % 1) a volInfo structure (the image does not have to exist on the path)
% % % 2) an image_vector, fmri_data, fmri_mask_image object
% % % 3) a string with the name of an image
% %
% % [image_name_to_sample_to, volInfo_to, imgobj_to] = define_sample_to_image(sampleto);
% %
% % fprintf('Resampling data from %3.0f images to space of %s\n', n_imgs, image_name_to_sample_to);
% %
% % % ---------------------------------------------------------
% % % If the images corresponding to data we're resampling don't exist,
% % % we must write them to disk
% % % ---------------------------------------------------------
% % obj = check_image_filenames(obj);
% %
% % if ~all(obj.files_exist)
% %     disp('Writing temporary image file for resampling. This will be removed.')
% %     fullpath_orig = obj.fullpath; % save the name for later
% %
% %     str = 'qwertyuiopasdfghjkl';
% %     str = str(randperm(length(str)));
% %     name = sprintf('tmp_image%s.img', str);
% %     obj.fullpath = fullfile(pwd, name);
% %     write(obj);
% % end
% %
% % newdat = zeros(volInfo_to.n_inmask, n_imgs, 'single');
% %
% % for i = 1:n_imgs
% %
% %     if size(obj.fullpath, 1) < i
% %         error('%3.0f images in obj.dat, but only %3.0f image names in obj.fullpath.', n_imgs, size(obj.fullpath, 1));
% %     end
% %
% %     % resample to new space and re-extract vector
% %     loadImg = deblank(obj.fullpath(i, :)); % ',' num2str(i)]); %**may be issues with 4-D vs. 3-D files
% %
% %     imgData = scn_map_image(loadImg, volInfo_to);
% %     newdat(:, i) = imgData(volInfo_to.wh_inmask);
% %
% % end
% %
% % obj.dat = newdat;
% % obj.volInfo = volInfo_to;
% %
% % if isa(obj, 'fmri_data')
% %     % fmri_data has this field, but other image_vector objects do not.
% %     obj.mask = resample_to_image_space(obj.mask, imgobj_to);
% % end
% %
% % if isa(obj, 'statistic_image')
% %     % statistic_image has this field, but other image_vector objects do not.
% %     obj.sig = ones(size(obj.dat));
% %     disp('.sig field reset. Re-threshold if necessary.');
% % end
% %
% % obj.history{end+1} = sprintf('Resampled data to space of %s', image_name_to_sample_to);
% %
% % % Clean up: remove temporary images
% % if ~all(obj.files_exist)
% %     str = sprintf('!rm %s*', obj.fullpath(1:end-4));
% %     disp('Removing temporary files:')
% %     disp(str)
% %     eval(str)
% %     obj.fullpath = fullpath_orig;
% % end
% %
% % % removed voxels must be updated due to resampling
% % obj.removed_voxels = false(obj.volInfo.n_inmask, 1);
% %
% % end % function
% %
% % %
% % %
% % %
% % %
% % %
% %
% %
% % function [image_name_to_sample_to, volInfo_to, imgobj_to] = define_sample_to_image(sampleto)
% %
% % switch class(sampleto)
% %     case 'char'
% %         image_name_to_sample_to = sampleto;
% %
% %         volInfo_to = iimg_read_img(sampleto, 2, 1, 1); % read data from file, first volume only
% %
% %         imgobj_to = fmri_mask_image(image_name_to_sample_to);
% %
% %     case {'fmri_mask_image'}
% %
% %         if ~isempty(sampleto.space_defining_image_name)
% %             image_name_to_sample_to = sampleto.space_defining_image_name;
% %         else
% %             image_name_to_sample_to = sampleto.volInfo.fname;
% %         end
% %
% %         volInfo_to = sampleto.volInfo;
% %         imgobj_to = sampleto;
% %
% %     case {'image_vector', 'fmri_data', 'statistic_image'}
% %         image_name_to_sample_to = sampleto.volInfo.fname;
% %         volInfo_to = sampleto.volInfo;
% %         imgobj_to = sampleto;
% %
% %     case 'struct'
% %         % assume its a volInfo structure
% %         image_name_to_sample_to = sampleto.fname;
% %         volInfo_to = sampleto;
% %
% %         if exist(image_name_to_sample_to, 'file')
% %             imgobj_to = fmri_mask_image(image_name_to_sample_to);
% %         else
% %             fprintf('Cannot resample to volInfo input struct because volInfo.fname image\n');
% %             fprintf('%s\n cannot be found.\n', image_name_to_sample_to);
% %         end
% %
% %     otherwise
% %         disp('fmri_mask_image.resample_to_image_space: illegal sampleto input.');
% %         disp('Must be volInfo structure or image_vector, fmri_data, fmri_mask_image object')
% %         disp('or a char name of an image.');
% %         error('exiting.');
% % end
% %
% % end
% %
% %
