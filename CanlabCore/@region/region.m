% 'region' is a class of objects that contains information on groups of
% voxels ("regions") defined in various ways (contiguous voxels above a
% threshold in an analysis, based on atlases, etc.)
%
% the region class replaces the "clusters" structure in the scnlab toolbox,
% and all or almost all of the "clusters" functions should work for regions
% structures as well.  The advantage of 'region' is that the class
% definition can help constrain the use and provide additional error
% checking, etc.
%
% Defining a region and initializing:
% ============================================
% Usage:
% r = region(obj1, [obj2], [keywords] )
%   obj1: [fmri_data/statistic_image object to define regions]
%         Can also be char array of image filename
%   obj2: Optional: fmri_data/statistic_image object to extract data from
%   keywords: Optional: 'unique_mask_values' or 'contiguous_regions'
%
% cl = region;  % generates an empty structure
%               % You can add fields yourself if you want to, but best to
%               define based on existing file name or image_vector object
%
% Define regions based on continuous voxel values (in this example image,
% there is only one set of contiguous voxels, so 1 region...)
%
% mask_image = which('brainmask.nii');
% cl = region(mask_image);
%
% INPUTS
% --------------------------------------------
% There are two ways to define which voxels are grouped into a 'region',
% which becomes an element of the region variable cl.
% enter 'contiguous_regions' -> group by contiguous blobs
% or    'unique_mask_values' -> group by unique values in mask .dat field
%
% Note: 'contiguous_regions' uses contiguity/clustering information stored
% in mask.volInfo.cluster, which may not have veridical contiguity info if
% you have manipulated it or incorporated anatomical information in working with
% the mask object, or if you have borrowed the volInfo structure from
% another source in creating it.
%
% Examples
% --------------------------------------------
% Define regions based on continuous values in an anatomical mask:
%
% mask_image = which('atlas_labels_combined.img');
% cl = region(mask_image, 'unique_mask_values');
%
% Resample mask_image to space of data_comb fmri_data object, and extract
% averages.  Space is defined by data_comb.
% cl = extract_roi_averages(data_comb, mask_image, 'unique_mask_values');
%
% Define regions based on unique voxels values in mask_image, and extract
% data stored in data_comb object, resampled to space of mask_image.
% Space is defined by mask_image:
% cl = region(mask_image, data_comb, 'unique_mask_values');
%
% Reslice mask to space of functional images and define regions based on
% mask values in the functional space (good for extracting data, etc.)
%
% mask_image = which('anat_lbpa_thal.img');
% mask = fmri_mask_image(mask_image);
% mask = resample_to_image_space(mask, image_names(1, :));
% cl = region(mask);
%
%
%
% Methods
% ============================================
% Try typing methods(cl)
%
% methods for 'regions' include:
% Visualization methods:
%     montage, orthviews, surf, etc.
%
%

classdef region
    
    properties
        title = 'Untitled';
        shorttitle = 'region';
        descrip1 = 'Region: a group of voxels.'
        descrip2 ='Some methods: orthviews, montage, surface, extract_data, table';
        
        XYZ
        XYZmm
        val
        val_descrip = 'Description of values for each voxel in val field.';
        
        Z
        Z_descrip = 'Legacy values for each voxel; Z-scores or other max stat';
        
        threshold % legacy, for compatibility
        voxSize
        M
        dim
        numVox
        numpeaks
        
        center
        mm_center
        
        timeseries
        contrastdata
        dat
        all_data
        
        source_images
        
        custom_info1
        custom_info1_descrip
        
        custom_info2
        custom_info2_descrip
        
    end % properties
    
    methods
        
        % class constructor method
        % takes maskinput in one of several forms:
        % 1) mask object
        % 2) mask image file
        % 3) mask vector data (???)
        %
        % takes one of two string inputs for how to define regions (see below)
        % 'contiguous_regions' (default)
        % 'unique_mask_values'
        
        function obj = region(maskinput, varargin)
            
            % initialize empty
            % return if no additional args
            
            obj.title = 'Untitled';
            obj.shorttitle = 'region';
            obj.descrip1 = 'Region: a group of voxels.';
            obj.descrip2 ='Some methods: orthviews, montage, surface, extract_data, table';
            obj.XYZ = [];
            obj.XYZmm = [];
            obj.val = [];
            obj.Z = [];
            obj.val_descrip = 'Values for each voxel, usually max stat.';
            obj.Z_descrip = 'Legacy values for each voxel; Z-scores or other max stat';
            obj.threshold = []; % legacy, for compatibility
            obj.voxSize = [];
            obj.M = [];
            obj.dim = [];
            obj.numVox = 0;
            obj.numpeaks = NaN;
            obj.center = [];
            obj.mm_center = [];
            obj.timeseries = [];
            obj.contrastdata = [];
            obj.dat = [];
            obj.all_data = [];
            obj.source_images = char([]);
            obj.custom_info1 = [];
            obj.custom_info1_descrip = 'Mask image name for region definition';
            obj.custom_info2 = [];
            obj.custom_info2_descrip = char([]);
            
            if nargin == 0
                return
            end
            
            % ---------------------------------
            % define mask object
            % ---------------------------------
            
            if isa(maskinput, 'char') % string file name
                mask = fmri_mask_image(maskinput);
                
            elseif isa(maskinput, 'image_vector')
                % case {'fmri_data', 'fmri_mask_image', 'statistic_image', 'image_vector'}
                
                % special for thresholded stats images: use threshold
                
                % if sig field exists, use sig voxels only
                if isa(maskinput, 'statistic_image') && ~isempty(maskinput.sig)
                    
                    for img = 1:size(maskinput.dat, 2)
                        maskinput.dat(~maskinput.sig(:, img), img) = 0;
                    end
                    
                end
                
                mask = maskinput;
                
            else
                error('region class constructor: unknown mask input type.')
            end
            clear maskinput
            
            % ---------------------------------
            % define what to average over
            % ---------------------------------
            
            average_over = 'contiguous_regions'; %'contiguous_regions'  or 'unique_mask_values';
            
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    switch varargin{varg}
                        
                        % reserved keywords
                        case 'contiguous_regions', average_over = 'contiguous_regions';
                        case 'unique_mask_values', average_over = 'unique_mask_values';
                            
                        otherwise
                            disp('region class constructor: Illegal string value for average_over.');
                            fprintf('You entered ''%s''\n Valid values are %s or %s\n', varargin{varg}, '''contiguous_regions''', '''unique_mask_values''');
                            error('Exiting');
                    end
                else
                    if isa(varargin{varg}, 'image_vector')
                        % data object to extract from at end
                        
                        dataobj = varargin{varg};
                        varargin{varg} = [];
                        
                        % Note: If you have manipulated an image_vector (e.g., fmri_data,
                        % statistic_image) object and eliminated some voxels, in order to create
                        % regions, contiguous voxels are automatically reparsed into regions using
                        % fmri_data.reparse_contiguous
                        
                        dataobj = reparse_contiguous(dataobj, 'nonempty');
                    end
                end
            end
            
            % Mask data can already be reduced to those indexed by
            % wh_inmask or not.
            % In addition, voxels/images with empty values may be removed.
            % So insert those first.
            % Note: .dat can sometimes have 2+ cols, so use only first one
            mask = replace_empty(mask);
            
            if size(mask.dat, 1) == mask.volInfo.n_inmask
                maskData = mask.dat(:, 1);
                
            elseif size(mask.dat, 1) == mask.volInfo.nvox
                % We have a full-length vector
                maskData = mask.dat(mask.volInfo.wh_inmask(:, 1));
            else
                error('Illegal size for mask.dat, because it does not match its volInfo structure.')
            end
            
            % If extracting data, we'll recreate the regions in
            % extract_roi_averages, so do that here and then exit.
            % Otherwise, continue to region definition.
            
            if exist('dataobj', 'var')
                % extract data
                disp('> Found image data, extracting region averages.');
                
                dataobj = replace_empty(dataobj); % may need to do this to get voxels to line up
                
                cs = compare_space(dataobj, mask);
                if cs == 3
                    error('Spaces for data object and mask object line up, but voxel numbers do not. Check.');
                elseif cs
                    disp('> Resampling to mask space first.');
                    dataobj = resample_space(dataobj, mask); % resample data to mask space
                end
                
                obj = extract_roi_averages(dataobj, mask, average_over);
                return
            end
            
            % ---------------------------------
            % get unique values for voxel grouping code
            % ---------------------------------
            maskValues = maskData; % values to save
            
            switch average_over
                
                % Define integer codes for sets of voxels to average over.
                
                case 'unique_mask_values'
                    
                    maskData = round(maskData);
                    u = unique(maskData)'; u(u == 0) = [];
                    nregions = length(u);
                    fprintf('Grouping voxels with unique mask values, assuming integer-valued mask: %3.0f regions\n', nregions);
                    
                    
                case 'contiguous_regions'
                    
                    isinmask = maskData ~= 0 & ~isnan(maskData);
                    
                    % re-make cluster ID for in-mask voxels
                    if sum(isinmask) < 50000
                        mask.volInfo.cluster(isinmask) = spm_clusters(mask.volInfo.xyzlist(isinmask, :)')';
                    else
                        % don't, and print warning
                        mask.volInfo.cluster(isinmask) = ones(sum(isinmask), 1);
                        disp('Warning: spm_cluster will not parse clusters for masks with > 50000 voxels.');
                    end
                    
                    u = unique(mask.volInfo.cluster(isinmask)); u(u == 0) = [];
                    maskData(isinmask) = mask.volInfo.cluster(isinmask);
                    nregions = length(u);
                    fprintf('Grouping contiguous voxels: %3.0f regions\n', nregions);
                    
                otherwise
                    error('This should never happen.');
            end
            
            % ---------------------------------
            % Now define the regions
            % ---------------------------------
            
            obj(1:nregions) = obj;  % fill in all empty fields
            
            for i = 1:nregions
                imgvec = maskData == u(i);
                
                obj(i).title = sprintf('Region %3.0f', i);
                obj(i).shorttitle = sprintf('Region%03d', i);
                obj(i).XYZ = mask.volInfo.xyzlist(imgvec,:)';
                obj(i).XYZmm = voxel2mm(obj(i).XYZ,mask.volInfo.mat);
                
                obj(i).val = maskValues(imgvec); %ones(1,size(obj(i).XYZ,2));
                obj(i).Z = maskValues(imgvec)'; %ones(1,size(obj(i).XYZ,2));
                
                obj(i).voxSize = abs(diag(mask.volInfo.mat(1:3, 1:3)));
                obj(i).M = mask.volInfo.mat;
                obj(i).dim = mask.volInfo.dim;
                obj(i).numVox = size(obj(i).XYZ, 2);
                
                obj(i).center = center_of_mass(obj(i).XYZ, obj(i).Z);
                obj(i).mm_center = center_of_mass(obj(i).XYZmm, obj(i).Z);
                
                [dd, ff, ee] = fileparts(mask.volInfo.fname);
                obj(i).custom_info1 = [ff ee];
                obj(i).custom_info1_descrip = 'Mask image name for region definition';
                
            end
            
            if nregions == 0
                % obligatory things even for empty regions
                obj(1).voxSize = abs(diag(mask.volInfo.mat(1:3, 1:3)));
                obj(1).M = mask.volInfo.mat;
                obj(1).dim = mask.volInfo.dim;
            end
            
        end % class constructor function
        
    end % methods
    
    
end