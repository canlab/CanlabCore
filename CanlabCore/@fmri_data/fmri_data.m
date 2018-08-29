% fmri_data: Data class for storing and analyzing neuroimaging datasets and associated information
%
% 'fmri_data' is a data class containing information about generic fmri
% datasets stored in a structure-like object.  Using this has the
% advantages that the properties and methods are standardized and controlled.
% It also keeps track of the history of what was done to the dataset.
%
% -------------------------------------------------------------------------
% Features and philosophy:
% -------------------------------------------------------------------------
%
% - Store image data in a flat (2-d), space-efficient voxels x images matrix 
% - Save space: Tools to remove out-of-image voxels and empty (zero/NaN) voxels, store data in single-precision format 
% - Analysis-friendly: 2-d matrices can be read and analyzed in multiple packages/algorithms
% - Meta-data included to convert back to 3-d image volume space, 
%   with easy tools (methods) to reconstruct and visualize images
%   Built-in resampling makes it easy to compare/combine datasets with different voxel sizes and image bounding boxes
%   Reduces overhead for statisticians/data scientists unfamiliar with neuroimaging to apply their algorithms.
% - Multiple images can be stored in a single object
% - Methods have short, intuitive names, and perform high-level functions specialized for neuroimaging:
% - Visualization (plot, orthviews, surface, montage, histogram, isosurface methods)
% - Image manipulation (apply_mask, get_wh_image, resample_space, compare_space, flip, threshold methods)
% - Data extraction (apply_atlas, apply_parcellation, extract_gray_white_csf, extract_roi_averages)
% - Analysis (ica, mahal, image_math, and many more in the fmri_data subclass)
% - Provenance: Ability to track and update history of changes to objects
%
% fmri_data is a subclass of the image_vector object and inherits its
% properties and methods.
%
% -------------------------------------------------------------------------
% To construct/create a new instance of an object:
% -------------------------------------------------------------------------
%
% Basic usage:
% obj = fmri_data(image_names, [maskinput], [other optional inputs]) 
%
% maskinput     :   [optional] name of mask image to use.  Default: 'brainmask.nii', a
%                   brain mask that is distributed with SPM software and in
%                   the CANlab Core Tools
%                   Alternative in CANlab tools: which('gray_matter_mask.img')
% 'noverbose'   :   Suppress verbose output during image loading
% 'sample2mask' :   Sample images to mask space. Default: Sample mask to
%                   image space, use native image space                            
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
% You can create an object by assembling an image_vector object from parts
% (entering fields) and converting using fmri_obj = fmri_data(image_vec_obj)
%
% You can create an fmri_data object with extacted image data.
% - Let "imgs" be a string array or cell array of image names
% - This command creates an object with your (4-D) image data:
% - fmri_dat = fmri_data(imgs);
% - Images can be zipped (.gz) as well. fmri_data() will unpack them.
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
% -----------------------------------------------------------------------
% Properties and methods
% -----------------------------------------------------------------------
% Properties are data fields associated with an object.
% Type the name of an object (class instance) you create to see its
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
% - regress(fmri_dat) runs multiple regression
% - predict(fmri_dat) runs cross-validated machine learning/prediction algorithms
%
% Key properties and methods (a partial list; type doc fmri_data for more):
% -------------------------------------------------------------------------
% fmri_data Properties (a partial list; type doc fmri_data for more):
%   dat                     - Image data, a [voxels x images] matrix, single-format
%   fullpath                - List of image names loaded into object with full paths 
%   history                 - History of object processing, for provenance 
%   image_names             - List of image names loaded into object, no paths  
%   removed_images          - Vector of images that have been removed (saves space; see remove_empty.m, replace_empty.m) 
%   removed_voxels          - Vector of empty in-mask voxels that have been removed (saves space; see remove_empty.m, replace_empty.m) 
%   volInfo                 - Structure with info on brain mask (saves space) and mapping voxels to brain space
%
% fmri_data Methods (a partial list; type doc fmri_data for more):
%   General:
% . 	descriptives        -  Get descriptives for an fmri_data or other image_vector object 
%       enforce_variable_types	-  Re-casts variables in objects into standard data types, which can save 
%       flip                - Flips images stored in an object left-to-right
%     	history             - Display history for image_vector object 
%       write               - Write an image_vector object to hard drive as an Analyze image (uses .fullpath field for image names)     
%   Data extraction:
%       apply_atlas         - Computes the mean value or pattern expression for each reference region specified in an atlas object 
%       apply_mask          - Apply a mask image (image filename or fmri_mask_image object) to an image_vector object 
%       apply_parcellation  - Computes the mean value or pattern expression for each parcel specified in a data object 
%       extract_gray_white_csf	- Extracts mean values (values) and top 5 component scores (components) 
%       extract_roi_averages	- This image_vector method a extracts and averages data stored in an fmri_data object 
%         
%   Handling brain space and image selection:
%       compare_space           - Compare spaces of two image_vector objects 
%       get_wh_image            - For an image_vector with multiple images (cases, contrasts, etc.), select a subset. 
%       reconstruct_image       - Reconstruct a 3-D or 4-D image from image_vector object obj 
%       remove_empty            - remove vox: logical vector of custom voxels to remove, VOX x 1 
%       reparse_contiguous      - Re-construct list of contiguous voxels in an image based on in-image 
%       replace_empty           - Replace empty/missing values in an image data object 
%       resample_space          - Resample the images in an fmri_data object (obj) to the space of another 
%
%   Display and visualization:
%   	display_slices      - Creates 3 separate montage views - ax, cor, sagg in a special figure window 
%       histogram           - Create a histogram of image values or a series of histograms for each 
%       image_similarity_plot - Associations between images in object and set of 'spatial basis function' images (e.g., 'signatures' or pre-defined maps)
%       isosurface          - Create and visualize an isosurface created from the boundaries in an image object. 
%       montage             - Create a montage of an image_vector (or statistic_image or fmri_data) 
%       orthviews               - display SPM orthviews for CANlab image_vector (or fmri_data, statistic_image) object 
%       pattern_surf_plot_mip	- axial maximum intensity projection pattern surface plot 
%   	sagg_slice_movie	- Movie of successive differences (sagittal slice)
%       slices              - Create a montage of single-slice results for every image in an image_vector object
%       surface             - Render image data on brain surfaces; options for cutaways and canonical surfaces
%       wedge_plot_by_atlas	- Plot a data object or 'signature' pattern divided into local regions 
%
%   Data processing and analysis:
%   	ica                 - Spatial ICA of an fmri_data object 
%       image_math          - Perform simple mathematical and boolean operations on image objects (see also plus, minus, power)
%       mahal               - Mahalanobis distance for each image in a set compared to others in the set
%       mean                - Mean across a set of images. Returns a new image_vector object. 
%       preprocess          - Preprocesses data in an image_vector (e.g., fmri_data) object; many options for filtering and outlier id 
%       qc_metrics_second_level	- Quality metrics for a 2nd-level analysis (set of images from different subjects) 
%       searchlight         - Run searchlight multivariate prediction/classification on an image_vector 
%       threshold           - Threshold image_vector (or fmri_data or fmri_obj_image) object based on raw threshold values
%       union               - ...and intersection masks for two image_vector objects 
%   
% -------------------------------------------------------------------------
% Examples and help:
% -------------------------------------------------------------------------
%
% To list properties and methods for this object, type:
% doc fmri_data, methods(fmri_data)
%
% Example 1: Load images (and run a simple analysis)
%
% Load a sample dataset into an fmri_data object (subclass of image_vector)
% This loads one of a set of named image collections used in demos/help:
% data_obj = load_image_set('emotionreg');
%
% You can load the same images manually, by locating the files, listing
% their names in a character array (or 1 x n cell array of strings), and
% then passing those into fmri_data:
%
% data_obj = fmri_data(which('Wager_2008_emo_reg_vs_look_neg_contrast_images.nii.gz'));
%
% filedir = what(fullfile('CanlabCore', 'Sample_datasets', 'Wager_et_al_2008_Neuron_EmotionReg'));
% image_names = filenames(fullfile(filedir.path, '*img'));
% data_obj = fmri_data(image_names);
%
% Now you can interact with the object.  Try, e.g.,:
% methods(data_obj)                               % List methods for object type
% descriptives(data_obj);                         % Print summary of descriptive statistics for the dataset
% plot(data_obj)                                  % Custom fmri_data specific plots
% t = ttest(data_obj);                            % Perform a voxel-wise one-sample t-test across images
% t = threshold(t, .005, 'unc', 'k', 10);         % Re-threshold with extent threshold of 10 contiguous voxels
% r = region(t);                                  % Turn t-map into a region object with one element per contig region
%
% Example 2: Extract data averaged over regions of interest:
%
% First run Example 1.  Now you have a thresholded t-statistic map.
% Extract averages (across voxels) for each subject in each contiguous
% region by typing:
%
% r = extract_roi_averages(data_obj, t);
%
% This returns r, which is another object--a "region"-class object.
% Region objects contain vectors, one element per pre-defined region in the
% image (in this case, significant blobs from our analysis). 
% Each element contains info that describes the region, including the voxels included,
% and average and voxel-by-voxel data if extracted from images and attached.   
% Type "doc region" for more info.  
% r has properties that hold multiple types of data:
% - .dat holds generic extracted data
% - .all_data holds voxel-by-voxel extracted data
%
% For more examples and walkthroughs, see the CANlab_help_examples
% repository at https://github.com/canlab/CANlab_help_examples
%
% Some example tutorials:
% canlab_help_1_installing_tools
% canlab_help_2_load_a_sample_dataset
% canlab_help_3_voxelwise_t_test_walkthrough
% canlab_help_4_write_data_to_image_file_format
% canlab_help_5_regression_walkthrough


% Programmers' notes:
% Tor Wager, 1/14/17 : Previously, if you stack names of images in different spaces and load them, 
% fmri_data did not return an error, but it returned distorted/incorrectly
% loaded images.  Now, it returns an error. For loading images in different
% spaces together, use the 'sample2mask' option.
% 
% Stephan Geuter, 5/16/17: Processing of masks in line 270 had been
% changed, masks were ignored. modified check in line 270 to include
% masking again

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
                            
                        case 'native_image_space' % do nothing, for convenience in calling scripts
                            
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
                
            elseif isempty(image_names)
               % We must have image names to run
               disp('fmridata: image_names is empty. No images to load.');
               return
               
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
                
                obj = run_checks_and_fixes(obj, verbosestr);
                
                return
                
            else
                % It's a file. Unzip if needed and check that image_names exist
                
                % Handle .gz by unzipping if needed
                [image_names, was_gzipped] = gunzip_image_names_if_gz(image_names);
                
            end
            
            % ---------------------------------
            % define mask object
            % ---------------------------------
            
            % if all keywords, like 'noverbose', then varargin will be empty
            % SG changed check, maskinput can be defined and varargin empty. 
            % maskinput is not part of varargin. 5/16/17                                                           
  
            % 5/16/17 Tor: if all varargin cells are empty, varargin may not
            % be... also, deal with special case where 2nd input is not
            % mask, but keyword.
            
            % Special case: 2nd argument is keyword, not mask.  This is
            % improper usage, but allow it for legacy reasons and
            % usability.
            if exist('maskinput', 'var') && ~isempty(maskinput) && ischar(maskinput)
                
                switch maskinput
                    
                    case 'verbose', maskinput = []; % nothing else needed
                    case 'noverbose', verbose = 0; verbosestr = 'noverbose'; maskinput = [];
                    case 'sample2mask', sample2mask = 1; maskinput = [];
                end
            end
            
            % Empty mask: use default
            if (nargin < 2 || isempty(maskinput)) % && isempty(varargin)  
                
                maskinput = which('brainmask.nii');
                if verbose, fprintf('Using default mask: %s\n', maskinput); end
                if isempty(maskinput), error('Cannot find mask image!'); end
                
            end
            
            switch class(maskinput)
                case 'char' % string file name
                    maskobj = fmri_mask_image(maskinput);
                    
                case {'fmri_mask_image'} % , 'fmri_data', 'statistic_image'} % others will not work with resample_to_space_defining... below
                        
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
                
            else % Not sample2mask, so resample mask to data image space
                
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
            
            % re-zip images if they were originally zipped
            % add .gz back to file names.
            if any(was_gzipped)
                
                image_names = re_zip_images(image_names, was_gzipped);
                
            end
            
            % add mask object to fmri_data object
            obj = create(obj, 'mask', maskobj, verbosestr);
            
            % append data
            obj = create(obj, 'dat', imgdat, 'image_names', image_names, verbosestr);
            clear imgdat
            
            % append description info
            [~, ff, ee] = fileparts(maskobj.volInfo.fname);
            mask_image_name = [ff ee];
            
            %obj = create(obj, 'dat_descrip', sprintf('Data from %s: %s', mask_image_name, obj.dat_descrip));
            obj = create(obj, 'mask_descrip', mask_image_name, verbosestr);
            
            obj.volInfo = maskobj.volInfo;
            
            obj.history(end+1) = {sprintf('Sampled to space of %s', maskobj.space_defining_image_name)};
            obj.history(end+1) = {['Masked with ' mask_image_name]};
            
            % Checks and fixes
            % -------------------------------------------------------------
            
            obj = run_checks_and_fixes(obj, verbosestr);
            
            
        end % constructor function
        
    end % methods
    
    
end

% -------------------------------------------------------------
% -------------------------------------------------------------

% Sub-functions

% -------------------------------------------------------------
% -------------------------------------------------------------

function obj = run_checks_and_fixes(obj, verbosestr)

obj = check_image_filenames(obj, verbosestr);

if isempty(obj.mask) && ~isempty(obj)  % isempty is an object method here
    % fix/create mask
    
    obj.mask.dat = any(obj.dat, 2);
    obj.mask.removed_voxels = obj.removed_voxels;
    obj.mask.volInfo_descrip = 'Generic mask built from .dat. Any voxel with a value in .dat is in-mask';
    
end

end % function

% -------------------------------------------------------------

function image_names_out = re_zip_images(image_names, was_gzipped)
% in : char, out: char

image_names_out = cellstr(image_names);

for i = 1:size(image_names, 1)
    
    if was_gzipped(i)
        
       % Use system to remove unzipped version after zipping.
       % Will wait for input, and not overwrite, if images exist
%         [status, result] = system(['gzip ' image_names(i, :)]);
        gzip(deblank(image_names(i, :)));
        image_names_out{i} = [deblank(image_names_out{i}) '.gz'];
        
    end
    
end

image_names_out = char(image_names_out{:}); 

end


