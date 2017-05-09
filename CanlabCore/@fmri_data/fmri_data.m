% fmri_data: Data class for storing data matrices and information
%
% 'fmri_data' is a data class containing information about generic fmri
% datasets stored in a structure-like object.  Using this has the
% advantages that the fields and methods are standardized and controlled.
% It also keeps track of the history of what was done to the dataset.
%
% Basic usage:
% obj = fmri_data(image_names, [maskinput], [other optional inputs]) 
%
% maskinput     :   Name of mask image to use.  Default: 'brainmask.nii', a
%                   brain mask that is distributed with SPM software
%                   Alternative in CANlab tools: which('gray_matter_mask.img')
% 'noverbose'   :   Suppress verbose output
% 'sample2mask' :   Sample images to mask space. Default: Sample mask to
%                   image space, use native image space
%
%
% Creating class instances
% -----------------------------------------------------------------------
% You can create an empty object by using:
% fmri_dat = fmri_data
% - fmri_dat is the object.
% - It will be created with a standard brain mask, brainmask.nii
% - This image should be placed on your Matlab path
% - The space information is stored in fmri_dat.volInfo
% - Data is stored in fmri_dat.dat, in a [voxels x images] matrix
% - You can replace or append data to the fmri_dat.dat field.
%
% You can create an fmri_data object with extacted image data.
% - Let "imgs" be a string array or cell array of image names
% - This command creates an object with your (4-D) image data:
% - fmri_dat = fmri_data(imgs);
% - Only values in the standard brain mask, brainmask.nii, will be included.
% - This saves memory by reducing the number of voxels saved dramatically.
%
% You can specify any mask you'd like to extract data from.
% - Let "maskimagename" be a string array with a mask image name.
% - this command creates the object with data saved in the mask:
% - fmri_dat = fmri_data(imgs, maskimagename);
% - The mask information is saved in fmri_dat.mask
%
% e.g., this extracts data from images within the standard brain mask:
% dat = fmri_data(imgs, which('brainmask.nii'));
%
% Properties and methods
% -----------------------------------------------------------------------
% Properties are data fields associated with an object.
% Try typing the name of an object (class instance) you create to see its
% properties, and a link to its methods (things you can run specifically
% with this object type). For example: After creating an fmri_data object 
% called fmri_dat, as above, type fmri_dat to see its properties.
%
% There are many other methods that you can apply to fmri_data objects to
% do different things.
% - Try typing methods(fmri_data) for a list.
% - You always pass in an fmri_data object as the first argument.
% - Methods include utilities for many functions - e.g.,:
% - resample_space(fmri_dat) resamples the voxels
% - write(fmri_dat) writes an image file to disk (careful not to overwrite by accident!)
% - regress runs multiple regression
% - predict runs cross-validated machine learning/prediction algorithms
%
% Defining the space of the extracted data
% -----------------------------------------------------------------------
% Note: There are two options for defining the space (i.e., coordinates/voxels)
% that the data is mapped to.
% By default, the mask is resliced to the same space as the first image in the
% input image name set (not coregistered; just resliced to the same voxel sizes.
% The images are assumed to be in register.)
%
% Reampling to mask space: YOU CAN ALSO map the image data to the space of the mask, by entering
% 'sample2mask' as in input argument.
% For loading images in different spaces together in one object, use the 'sample2mask' option.
%
% Attaching additional data
% -----------------------------------------------------------------------
% The fmri_data object has a number of fields for appending specific types of data.
%
% - You can replace or append data to the fmri_dat.dat field.
% - The fmri_data object will also store predictor data (.X) also outcome data (.Y)
% - There are many fields for descriptions, notes, etc., like "dat_descrip" and "source_notes"
% - Attach custom descriptions in these fields to document your object.
% - The "history" field stores a cell array of strings with the processing
% history of the object. Some methods add to this history automatically.
%
% Extracting ROI data easily
% -----------------------------------------------------------------------
% You can extract image data, and save averages within regions of
% interest, by doing something like this:
% ****[fmri_dat, cl] = read_image_files(image_names(3, :));
%
% cl is the ROI data in a region object.
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

% Programmers' notes:
% Tor Wager, 1/14/17 : Previously, if you stack names of images in different spaces and load them, 
% fmri_data did not return an error, but it returned distorted/incorrectly
% loaded images.  Now, it returns an error. For loading images in different
% spaces together, use the 'sample2mask' option.

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
            
            % DEFAULT INPUTS
            % -----------------------------------
            
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
            
            if nargin == 0
                % SPECIAL DEFAULT METHOD: Standard empty mask
                % -----------------------------------
                
                % Empty: Define with standard default mask
                [image_names, maskinput] = deal(which('brainmask.nii'));
                
                if isempty(maskinput)
                    disp('Warning: Cannot find brainmask.nii, creating without mask info.');
                    return
                end
                
                verbose = 0; 
                verbosestr = 'noverbose';
                
            end
            
            
            % ---------------------------------
            % Special: if existing image_vector
            % map into space of fmri_data and return
            % ---------------------------------
            if iscell(image_names), image_names = char(image_names{:}); end
            
            
            if isa(image_names, 'image_vector')
                % Map fields of input object into fmri_data structure
                
                warning off
                obj2 = struct(image_names); %image_names;  % tor: struct not needed i think. nope, is needed.
                warning on
                
                N = fieldnames(obj);
                for i = 1:length(N)
                    if isfield(obj2, (N{i}))
                        obj.(N{i}) = obj2.(N{i});
                    end
                end
                
                obj.mask.volInfo = obj2.volInfo;
                return
                
            else
                
                % Handle .gz by unzipping if needed
                image_names = gunzip_image_names_if_gz(image_names);
                
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
                    
                    case {'SPM12','SPM8', 'SPM5'}
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

