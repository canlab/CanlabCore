% fmri_data: Data class for storing data matrices and information
%
% 'fmri_data' is a data class containing information about generic fmri
% datasets stored in a structure-like object.  Using this has the
% advantages that the fields and methods are standardized and controlled.
% It also keeps track of the history of what was done to the dataset.
%
% Creating class instances
% -----------------------------------------------------------------------
% You can create an empty object by using:
% fmri_dat = fmri_data
% (fmri_dat is the object)
%
% You can create an object and extract data from a mask (defining many of
% the fields in the object) like this:
% 
% dat = fmri_data(imgs, maskimagename);
% 
% e.g.,
% dat = fmri_data(imgs, which('brainmask.nii'));
%
% Defining the space of the extracted data
% -----------------------------------------------------------------------
% Note: There are two options for defining the space (i.e., coordinates/voxels)
% that the data is mapped to.
% By default, the mask is resliced to the same space as the first image in the
% input image name set (not coregistered; just resliced to the same voxel sizes.
% The images are assumed to be in register.)
% YOU CAN ALSO map the image data to the space of the mask, by entering
% 'sample2mask' as in input argument.
%
% Creating class instances
% -----------------------------------------------------------------------
% The fmri_data object will store image data (.X) also outcome data (.Y)
% Try typing the name of an object (class instance) you create to see its
% properties, and a link to its methods (things you can run specifically
% with this object type).
%
% Extracting ROI data easily
% -----------------------------------------------------------------------
% You can extract image data, and save averages within regions of
% interest, by doing something like this:
% [fmri_dat, cl] = read_image_files(image_names(3, :));
%
% cl is the ROI data in a region object 
% region is a class.  It's data structure is like the older "clusters"
% structure format, with average data values stored in cl.dat
% regions can be defined by EITHER contiguous voxels or based on unique integer
% values in images.
%
% More about working with masks
% -----------------------------------------------------------------------
% You need a mask image to define which voxels are extracted and possibly the
% space of the image data.
% If you do not yet have a mask image, but have data extracted separately,
% you can add mask information (from a mask in the same space) like this:
%
% dat = create(dat, 'mask', fmri_mask_image(maskimg));
%
% More methods
% -----------------------------------------------------------------------
%
% Methods include create, extract_roi_averages
%
% Create lets you add fields/data to a structure (see help
% fmri_data.create)
%
% extract_roi_averages lets you specify a new mask, and extract and average
% data from ROIs defined by the new mask (provided they were in the
% original mask from which you extracted data!)
%
% Examples:
% obj = fmri_data(image_names, maskinput)
% obj = fmri_data(image_names, [], 'noverbose')


classdef fmri_data < image_vector
    
    properties
        % also inherits the properties of image_vector.
        
        source_notes = 'Source notes...';

        X % legacy; temporary, so we can load old objects
        
        mask = fmri_mask_image;
        mask_descrip = 'mask is an fmri_mask_image object that defines the mask.';

        images_per_session
        
        Y = [];
        Y_names;
        Y_descrip = 'Behavioral or outcome data matrix.';
        
        covariates;
        covariate_names = {''};
        covariates_descrip = 'Nuisance covariates associated with data';
        
        history_descrip = 'Cell array: names of methods applied to this data, in order';
        
        additional_info = struct('');
        
    end % properties
    
    methods
        
        % Class constructor
        function obj = fmri_data(image_names, maskinput, varargin)
            %
            % [obj, cl_with_averages] = fmri_data(image_names, mask_image, varargin)
            %
            % Reads a set of image files and a mask image, and returns
            % an fmri_data object with data for all in-mask voxels.
            
            % ---------------------------------
            % Create empty fmri_data object, and return if no additional
            % arguments
            % ---------------------------------
            

            obj.source_notes = 'Info about image source here';
            obj.mask = fmri_mask_image;
            obj.mask_descrip = 'Volume and in-area mask info from iimg_read_img';

            obj.X = []; % legacy; temporary, so we can load old objects
            obj.Y = [];
            obj.Y_names;
            obj.Y_descrip = 'Behavioral or outcome data matrix.';
            obj.covariates;
            obj.covariate_names = {''};
            obj.covariates_descrip = 'Nuisance covariates associated with data';
            
            obj.images_per_session = [];
            
            obj.history = {''};
            obj.history_descrip = 'Cell array of names of methods applied to this data, in order';
            obj.additional_info = struct('');
            
            if nargin == 0
                return
            end
            
            verbose = 1;
            verbosestr = 'verbose';
            sample2mask = 0;
            
            % SET UP OPTIONAL INPUTS
            % -----------------------------------
            
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        
                        case 'verbose', varargin{i} = []; % nothing else needed
                        case 'noverbose', verbose = 0; verbosestr = 'noverbose'; varargin{i} = [];
                        case 'sample2mask', sample2mask = 1; varargin{i} = [];
                            
                        otherwise, warning(['Unknown input string option:' varargin{i}]);
                    end
                end
            end
            
            % ---------------------------------
            % Special: if existing image_vector
            % map into space of fmri_data and return
            % ---------------------------------
            if iscell(image_names), image_names = char(image_names{:}); end
            
            if isa(image_names, 'image_vector')
                % Map fields of input object into fmri_data structure
                
                obj2 = struct(image_names);
                
                N = fieldnames(obj);
                for i = 1:length(N)
                    if isfield(obj2, (N{i}))
                        obj.(N{i}) = obj2.(N{i});
                    end
                end
                
                obj.mask.volInfo = obj2.volInfo;
                return
            end
            
            % ---------------------------------
            % define mask object
            % ---------------------------------
            
            if nargin < 2 || isempty(maskinput)
                maskinput = which('brainmask.nii');
                if verbose, fprintf('Using default mask: %s\n', maskinput); end
                if isempty(maskinput), error('Cannot find mask image!'); end
            end
            
            switch class(maskinput)
                case 'char' % string file name
                    maskobj = fmri_mask_image(maskinput);
                    
                case 'fmri_mask_image'
                    maskobj = maskinput;
                    
                otherwise
                    error('region class constructor: unknown mask input type.')
            end
            clear maskinput
            
            % Either extract the image data in the space of the mask, or
            % vice versa
            %sample2mask = strmatch('sample2mask', varargin);
            
            if sample2mask
                % Read data in mask space; map images to mask
                % ------------------------------------------------ 
                if verbose, fprintf('Expanding image filenames if necessary\n'); end
                
                for i = 1:size(image_names, 1)
                    
                    iinames{i} = expand_4d_filenames(image_names(i, :));
                end
                iinames = char(iinames{:});
                imgdat = zeros(length(maskobj.volInfo.wh_inmask), size(iinames, 1), 'single');
                
                if isempty(iinames)
                    disp('Images do not exist!'); disp(image_names); error('Exiting');
                end
                
                % Now extract the actual data from the mask
                if verbose
                    fprintf('Sampling %3.0f images to mask space: %04d', size(iinames, 1), 0)
                end
                
                for i = 1:size(iinames, 1)
                    
                    if verbose, fprintf('\b\b\b\b%04d', i); end
                    idat = scn_map_image(iinames(i, :), maskobj);
                    imgdat(:, i) = idat(maskobj.volInfo.wh_inmask);
                    
                end
                
                if verbose, fprintf('\n'); end
                
                [dd, ff, ee] = fileparts(maskobj.volInfo.fname);
                maskobj.space_defining_image_name = [ff ee];
                
                obj.volInfo = maskobj.volInfo;
                
            else
                % Read data in image space; map mask to images
                % ------------------------------------------------
                
                % resample to image space if necessary
                space_defining_image = deblank(image_names(1, :));
                maskobj = resample_to_image_space(maskobj, space_defining_image);
                
                
                % Now extract the actual data from the mask   
                switch spm('Ver')
                    
                    case {'SPM8', 'SPM5'}
                        imgdat = iimg_get_data(maskobj.volInfo, image_names, 'single', verbosestr, 'noexpand');
                        
                    case {'SPM2', 'SPM99'}
                        % legacy, for old SPM
                        imgdat = iimg_get_data(maskobj.volInfo, image_names, 'single', verbosestr);
                        
                    otherwise
                        error('Unknown version of SPM! Update code, check path, etc.');
                end
                
                imgdat = imgdat';
                
            end % read data, depending on mask sampling
            
            % add mask object to fmri_data object
            obj = create(obj, 'mask', maskobj, verbosestr);
            
            % append data
            obj = create(obj, 'dat', imgdat, 'image_names', image_names, verbosestr);
            clear imgdat
            
            % append description info
            [dd, ff, ee] = fileparts(maskobj.volInfo.fname);
            mask_image_name = [ff ee];
            %obj = create(obj, 'dat_descrip', sprintf('Data from %s: %s', mask_image_name, obj.dat_descrip));
            obj = create(obj, 'mask_descrip', mask_image_name, verbosestr);

            obj.volInfo = maskobj.volInfo;
            
            obj.history(end+1) = {sprintf('Sampled to space of %s', maskobj.space_defining_image_name)};
            obj.history(end+1) = {['Masked with ' mask_image_name]};
            
            obj = check_image_filenames(obj, verbosestr);
            
        end % constructor function
        
    end % methods
    
    
end