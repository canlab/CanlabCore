classdef image_vector
% image_vector: An object that allows for storage and manipulation of neuroimaging data
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
% image_vector is a superclass of several other object types: 
% fmri_data
% statistic_image
% atlas
%
% These object types have additional analysis methods and properties, and
% are most often used. You will rarely create an image_vector object
% directly. Typically you would create instances of one of the classes above.
%
% -------------------------------------------------------------------------
% To construct/create a new instance of an object:
% -------------------------------------------------------------------------
%
% As with any neuroimaging-data object class in this toolbox 
% (image_vector, fmri_data, statistic_image, atlas), you can create an object by
% specifying the names of a set of images you want to load in and passing them into the
% function with the same name as the object class (e.g., fmri_data). This
% will create an object with image data and image space information attached.
%
% Use fmri_data instead of image_vector to load image data:
% data_obj = fmri_data(image_names);
%
% ...where image_names is a character array or 1 x n cell array of strings
%
% You can also manually construct an object by entering its constituent
% fields. For the object to be valid, these much match one another in size/space/etc, and you
% may get errors later if they do not match. 
%
% For example, this code first loads a sample atlas object and then creates an image_vector object
% from its fields. Then, we re-cast it as an fmri_data object, which is
% generally more useful and has more associated methods.
%
% obj = load_atlas('cerebellum');  
% obj_with_region_indices = image_vector('dat', single(obj.dat), 'volInfo', obj.volInfo, 'removed_voxels', obj.removed_voxels, 'removed_images', obj.removed_images, 'image_names', obj.image_names, 'noverbose');
% obj_with_region_indices = fmri_data(obj_with_region_indices);
% orthviews(obj_with_region_indices)
%
% -------------------------------------------------------------------------
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
% Key properties and methods (a partial list; type doc image_vector for more):
% -------------------------------------------------------------------------
% image_vector Properties:
%   dat                     - Image data, a [voxels x images] matrix, single-format
%   fullpath                - List of image names loaded into object with full paths 
%   history                 - History of object processing, for provenance 
%   image_names             - List of image names loaded into object, no paths  
%   removed_images          - Vector of images that have been removed (saves space; see remove_empty.m, replace_empty.m) 
%   removed_voxels          - Vector of empty in-mask voxels that have been removed (saves space; see remove_empty.m, replace_empty.m) 
%   volInfo                 - Structure with info on brain mask (saves space) and mapping voxels to brain space
%
% image_vector Methods:
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
% doc image_vector, methods(image_vector)
%
% Example 1: 
% Load a sample dataset into an fmri_data object (subclass of image_vector)
% This loads one of a set of named image collections used in demos/help:
% data_obj = load_image_set('emotionreg');
%
% You can load the same images manually, by locating the files, listing
% their names in a character array (or 1 x n cell array of strings), and
% then passing those into fmri_data:
%
% filedir = what(fullfile('CanlabCore', 'Sample_datasets', 'Wager_et_al_2008_Neuron_EmotionReg'));
% image_names = filenames(fullfile(filedir.path, '*img'));
% data_obj = fmri_data(image_names);
%
% Now you can interact with the object.  Try, e.g.,:
% methods(data_obj)                               % List methods for object type
% plot(data_obj)                                  % Custom fmri_data specific plots
% t = ttest(data_obj);                            % Perform a voxel-wise one-sample t-test across images
% t = threshold(t, .005, 'unc', 'k', 10);         % Re-threshold with extent threshold of 10 contiguous voxels
% r = region(t);                                  % Turn t-map into a region object with one element per contig region
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
%
% License and authors: 
% -------------------------------------------------------------------------
% By Tor Wager and CANlab members and collaborators
% GNU General Public License; see http://www.gnu.org/licenses/ 
% See Programmers' Notes in code for more info
% 
%
% -------------------------------------------------------------------------
% See also: 
% fmri_data
% statistic_image
% atlas
% region

% ..
%     Author and copyright information:
%
%     Copyright (C) 2010 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..


    properties
        
        dat                             % Image data, a [voxels x images] matrix (single datatype; saves space)
        dat_descrip                     % String description of the dataset
        
        % Structure with info on brain mask (saves space) and mapping voxels to brain space
        % 'mat' field contains origin, voxel sizes and affine transform matrix
        % xyzlist contains voxel coordinates for in-mask voxels (see voxel2mm.m)
        % cluster contains indices for contiguous clusters
        volInfo                     
        
         % Vector of empty in-mask voxels that have been removed (saves space; see remove_empty.m, replace_empty.m)
        removed_voxels = logical(0);   

        % Vector of images that have been removed (saves space; see remove_empty.m, replace_empty.m)
        removed_images = logical(0);
                
        image_names                     % List of image names loaded into object, no paths  
        
        % List of image names loaded into object with full paths
        % write(image_vector) uses this field for filenames to write to
        fullpath 
        
        files_exist                     % Check of whether files exist on disk, logical vector
        history                         % History of object processing, for provenance
         
    end
    
    methods
        
        % Class constructor
        function obj = image_vector(varargin)
            % Enter fieldname', value pairs in any order to create class
            % instance
            
            % ---------------------------------
            % Create empty image_vector object, and return if no additional
            % arguments
            % ---------------------------------
            
            obj.dat = [];
            obj.volInfo = [];
            obj.image_names = [];
            obj.fullpath = char([]);
            obj.files_exist = false;
            obj.history = {};
            
            % The code below can be generic to any class definition
            % It parses 'fieldname', value pairs of inputs
            % and returns a warning if unexpected strings are found.
            
            if nargin == 0
                return
            end
            
            % all valid fieldnames
            valid_names = fieldnames(obj);
            
            control_args = {};
            
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    
                    % known control strings (ignore but pass in)
                    if strcmp(varargin{i}, 'noverbose')
                        control_args{end + 1} = varargin{i};
                        
                        continue % skip the rest in this loop
                    end
                    
                    % Look for a field (attribute) with the input name
                    wh = strmatch(varargin{i}, valid_names, 'exact');
                    
                    % behaviors for valid fields
                    if ~isempty(wh)
                        
                        obj.(varargin{i}) = varargin{i + 1};
                        
                        % eliminate strings to prevent warnings on char
                        % inputs
                        if ischar(varargin{i + 1})
                            varargin{i + 1} = [];
                        end
                        
                        % special methods for specific fields
                        switch varargin{i}
                            
                        end
                        
                    else
                        warning('inputargs:BadInput', sprintf('Unknown field: %s', varargin{i}));
                        
                    end
                end % string input
            end % process inputs
            
            % load data, if we don't have data yet
            % But we do have valid image names.
            % ----------------------------------------
            
            obj = check_image_filenames(obj, control_args{:});
            
            if isempty(obj.dat) && any(obj.files_exist)
                obj = read_from_file(obj);
                
            elseif isempty(obj.dat)
                disp('Warning: .dat is empty and files cannot be found.  No image data in object.');
            end
            
                
        end % class constructor
        
    end  % methods
    
end  % classdef

