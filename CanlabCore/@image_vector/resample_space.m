function obj = resample_space(obj, sampleto, varargin)
% Resample the images in an fmri_data object (obj) to the space of another
% image (sampleto; e.g., a mask image). Works for all image_vector objects.
% The object includes only voxels in the in-mask region in the target
% (sampleto) image.
%
% :Usage:
% ::
%
%    obj = resample_space(obj, sampleto, [sampling method])
%
% Sampleto may be one of these:
%   1. a volInfo structure (the image does not have to exist on the path)
%   2. an image_vector, fmri_data, fmri_mask_image object
%   3. a string with the name of an image
%
% Can enter resampling method as optional input. Takes any input to
% interp3:
%       'nearest' - nearest neighbor interpolation
%       'linear'  - linear interpolation (default)
%       'spline'  - spline interpolation
%       'cubic'   - cubic interpolation as long as the data is uniformly
%                   spaced, otherwise the same as 'spline'
%
% :Examples:
% ::
%
%    label_mask = fmri_data(which('atlas_labels_combined.img'));
%    label_mask = resample_space(label_mask, ivec, 'nearest') % resamples and masks label image
%
% ..
%    Programmers' notes:
%
%    1/27/2012 Tor edited to handle .mask field in fmri_data and .sig field in
%    statistic_image.  Was causing errors otherwise...
%           Also changed automatic behavior to reparse contig voxels with
%           'nonempty' in output obj
%
%   2/5/2018    Tor added support for atlas objects - special handling
%
%   2/23/2018   Stephan changed replace_empty(obj) into replace_empty(obj,'voxels') to prevent adding removed images back in
%   5/18/2021   Tor removed line: obj_out = replace_empty(obj_out); for
%   statistic_image objects, as it was causing a voxel mismatch...incorrect
%   removed_voxels due to partially built object
%   10/28/2024   Zizhuang changed the default method to resample .sig field
%   to nearest neighbor, and add related warnings
% ..

n_imgs = size(obj.dat, 2);

if ischar(sampleto)
    sampleto = fmri_data(sampleto);
end

Vto = sampleto.volInfo;
SPACEto = define_sampling_space(Vto, 1);

Vfrom = obj.volInfo;
SPACEfrom = define_sampling_space(Vfrom, 1);

obj_out = obj;
obj_out.dat = [];
obj_out.volInfo = Vto;

obj = replace_empty(obj,'voxels');  % to make sure vox line up
% changed to 'voxels' only. SG 2/23/18

if ~isa(obj, 'atlas')
    
    % Standard image_vector objects
    % -----------------------------------------------------------------------
    
    for i = 1:n_imgs
        
        voldata = iimg_reconstruct_vols(obj.dat(:, i), obj.volInfo);
        
        resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, varargin{:});
        
        resampled_dat = resampled_dat(:);
        
        obj_out.dat(:, i) = resampled_dat(Vto.wh_inmask);
        
    end
    
    % in case of NaN values
    obj_out.dat(isnan(obj_out.dat)) = 0;
    
    
    % Special object subtypes
    % -----------------------------------------------------------------------
    
else % if  isa(obj, 'atlas')
    
    n_prob_imgs = size(obj.probability_maps, 2);
    
    obj_out.probability_maps = [];
    
    % Use probability images if available
    
    for i = 1:n_prob_imgs
        
        voldata = iimg_reconstruct_vols(obj.probability_maps(:, i), obj.volInfo);
        
        resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, varargin{:});
        
        resampled_dat = resampled_dat(:);
        
        obj_out.probability_maps(:, i) = resampled_dat(Vto.wh_inmask);
        
    end
    
    % rebuild .dat from probability images - done below
    %     if n_prob_imgs
    %         obj_out = probability_maps_to_region_index(obj_out);
    %     end
    
    % if no prob images, need to be careful about how to resample integer vector data
    
    if ~n_prob_imgs
        
        % integer_vec = zeros(Vto.n_inmask, 1);
        
        n_index_vals = length(unique(obj.dat(obj.dat ~= 0)));
        
        % create a set of pseudo-"probabilities" for each region, resampled. Then
        % we can take the max prob, so that each voxel gets assigned to the best-matching parcel.
        
        pseudo_prob = zeros(Vto.n_inmask, n_index_vals);
        
        for i = 1:n_index_vals
            
            %myintegervec = i * double(obj.dat(:, 1) == i);
            
            myintegervec = double(obj.dat(:, 1) == i); % 1/0 "pseudo-probability"
            
            voldata = iimg_reconstruct_vols(myintegervec, obj.volInfo);
            
            resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, varargin{:});
            
            resampled_dat = resampled_dat(:);
            resampled_dat = resampled_dat(Vto.wh_inmask);    % take relevant voxels only
            % resampled_dat(~(round(resampled_dat) == i)) = 0; % take only values that round to integer
            
            pseudo_prob(:, i) = resampled_dat;
            %integer_vec = integer_vec + round(resampled_dat);
            
        end
        
        obj_out.probability_maps = pseudo_prob;
        
        % obj_out.dat = integer_vec; % will be rounded later, but should be rounded already here...
        
    end % rebuild integers
    
    % rebuild .dat from probability images
    obj_out = probability_maps_to_region_index(obj_out);
    
end % atlas object


if isa(obj_out, 'statistic_image')
    % Rebuild fields specific to statistic_images
    
%     obj_out = replace_empty(obj_out); % TOR REMOVED 5/18/21, AS IT
%     RESULTS IN VOXEL LIST MISMATCH WITH PARTIALLY BUILT OBJECT FIELDS

    k = size(obj_out.dat, 2);
    
    [obj_out.p, obj_out.ste, obj_out.sig, obj_out.N] = deal([]); % these will be resampled
    
    p = ones(obj.volInfo.nvox, k);
    ste = Inf .* ones(obj.volInfo.nvox, k);
    sig = zeros(obj.volInfo.nvox, k);
    N = zeros(obj.volInfo.nvox, 1);
    
    for i = 1:k
        % this may break if nvox (total in image) is different for 2
        % images...
        
        % Wani: in some cases, obj could have empty p, ste, and sig
        % Tor, 4/2021. Need to resample these appropriately, handle
        % logicals, and add N field
        
        if ~isempty(obj.p)
            
            p(obj.volInfo.wh_inmask, i) = obj.p(:, i);
            
            voldata = iimg_reconstruct_vols(p(:, i), obj.volInfo);
            resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, varargin{:});
            resampled_dat = resampled_dat(:);
            obj_out.p(:, i) = resampled_dat(Vto.wh_inmask);
        end
        
        if ~isempty(obj.ste)
            
            ste(obj.volInfo.wh_inmask, i) = obj.ste(:, i);
            
            voldata = iimg_reconstruct_vols(ste(:, i), obj.volInfo);
            resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, varargin{:});
            resampled_dat = resampled_dat(:);
            obj_out.ste(:, i) = resampled_dat(Vto.wh_inmask);
            
        end
        
        % For .sig, we must convert from logical and threshold back to logical
        % Can't interpolate logical vectors
        if ~isempty(obj.sig)
            
            % ---------------------------------------------
            % Added by Zizhuang Miao 10/28/2024
            % the default linear interpolation during resampling
            % could render .sig field in statistic images invalid,
            % because it could make many voxels with an original .sig = 0 into
            % having a .sig = 1 as long as it is close to a significant voxel
            % generally we won't recommend resampling .sig field
            % give a warning on that
            warning(['Resampling voxel significance can cause false positives or negatives. ' ...
                'Consider resampling data before running statistical analysis.']);
            % ---------------------------------------------

            sig(obj.volInfo.wh_inmask, i) = double(obj.sig(:, i));
            voldata = iimg_reconstruct_vols(sig(:, i), obj.volInfo);

            % -----------------------------------------
            % Edited by Zizhuang Miao 10/28/2024
            % nearest neighbor method can limit false positives,
            % so use that as default unless otherwise specified
            if varargin{:} == 'linear'
                % if users specify using linear interpolation
                warning(['Linear interpolation will lead to many false positives. ' ...
                    'Consider using the default nearest neighbor method.']);
                resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, varargin{:});
            else
                % default to nearest neighbor
                warning(['Using nearest neighbor method to resample .sig field. ' ...
                    'This can limit false positives, but can also cause a mismatch ' ...
                    'between .sig and other fields like .p.'])
                resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, 'nearest');
            end
            % ------------------------------------------

            resampled_dat = resampled_dat(:);
            resampled_dat(isnan(resampled_dat)) = 0;
            obj_out.sig(:, i) = logical(resampled_dat(Vto.wh_inmask));
        end
        
    end
    
        if ~isempty(obj.N)
            
            N(obj.volInfo.wh_inmask, 1) = obj.N(:, 1);
            
            voldata = iimg_reconstruct_vols(N, obj.volInfo);
            resampled_dat = interp3(SPACEfrom.Xmm, SPACEfrom.Ymm, SPACEfrom.Zmm, voldata, SPACEto.Xmm, SPACEto.Ymm, SPACEto.Zmm, varargin{:});
            resampled_dat = resampled_dat(:);
            obj_out.N = resampled_dat(Vto.wh_inmask);
            
        end
        
%     if ~isempty(obj.p), obj_out.p = p(Vto.wh_inmask, :); end
%     if ~isempty(obj.ste), obj_out.ste = ste(Vto.wh_inmask, :); end
%     if ~isempty(obj.sig), obj_out.sig = sig(Vto.wh_inmask, :); end
    
end

% End special object subtypes
% -----------------------------------------------------------------------


% Handle removed voxels
% -----------------------------------------------------------------------

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


if isempty(obj_out)
    % return - nothing more to do
    return   
end

% This stuff below added 1/27/13 by tor

% re-parse clusters
obj = reparse_contiguous(obj, 'nonempty');

obj.history{end+1} = sprintf('Resampled data to space of %s', sampleto.volInfo.fname);

% Special object subtypes
% -----------------------------------------------------------------------

if isa(obj, 'fmri_data')
    % fmri_data has this field, but other image_vector objects do not.
    obj.mask = resample_space(obj.mask, sampleto);
end

% if isa(obj, 'statistic_image')
%     % statistic_image has this field, but other image_vector objects do not.
%     obj.sig = ones(size(obj.dat));
%     disp('.sig field reset. Re-threshold if necessary.');
% end

if isa(obj, 'atlas')
    % Rebuild index so we have integers only. Rebuild if we have prob maps,
    % or round .dat if not.
    n_regions = max([size(obj.probability_maps, 2) length(obj.labels)]); % edit from num_regions method to exclude using .dat, as we are trying to adjust dat for interpolation
    
    has_pmaps = ~isempty(obj.probability_maps) && size(obj.probability_maps, 2) == n_regions;
    
    if has_pmaps
        obj = probability_maps_to_region_index(obj);
    else
        obj.dat = int32(round(obj.dat));
    end
    
    [obj, ~, ~, missing_regions] = check_properties(obj, 'compress_index');  % check. adjust indices and print warning if we have lost regions
    
    if any(missing_regions)
        disp('Some atlas regions lost in resampling:');
        disp(missing_regions');
    end
    
end

% End special object subtypes
% -----------------------------------------------------------------------



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
