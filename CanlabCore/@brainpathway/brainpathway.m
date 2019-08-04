% brainpathway: Data class for storing and analyzing brain connections and pathways
%
% -------------------------------------------------------------------------
% Features and philosophy:
% -------------------------------------------------------------------------
%
% 'brainpathway' is a data class containing information about fmri or other
% neuroimaging datasets stored in a structure-like object.  Using this has the
% advantages that the properties and methods are standardized and controlled.
% It also keeps track of the history of what was done to the dataset.
% Brainpathway objects have specialized methods for calculating
% connectivity in a priori defined regions or networks and across large
% sets of regions. Some relatively unique features include the capacity to:
% (1) create partitions among sets of regions (e.g, cortex, basal ganglia)
% to restrict connectivity estimation and plotting to pathways that connect
% sets (e.g., cortico-striatal loops).
% (2) create multiple, potentially overlapping networks
% (3) separate the definition of nodes and regions, so that a region may
% have multiple nodes, each defined as a (potentially) complex function of
% voxel data in a region. e.g., a region may have two nodes defined by
% different multivariate patterns across the region's voxels.
%
% A brainpathway object has these parts:
% - region_atlas:                An atlas-class object defining k regions
% - WEIGHTS object ->  n fmri_data objects, one per network,
% - connections:            A series of matrices specifying [k x k] bivariate connections
% - connections.apriori:    [1/0] logical matrices specifying existing connections, k x k x n for n networks
% k x n latent variable weights - voxel weights. {1....k} cell k has {v x n}, v voxels x n networks
% OR: weights could be in n fmri_data objects, one per network, with
% - connections.est:        Estimated connectivity strengths
% - connections.se:         Standard error of estimated connectivity strengths
% - connections.metric      Metric type [r, cos_sim, tau, partial_r]
% - graph_properties:
% - partitions:             An integer vector of partition labels for each
%                           node, which define blocks of nodes. This can be used to, say, identify nodes in Block A maximally connected to those in Block B.
%
% - DATA STORAGE:
%   brainpathway objects can potentially store data at several spatial scales, and 
%   also represent the multi-level (hierarchical) nature of data typical
%   for neuroimaging experiments. Spatial scales include, in approximate order from fine to coarse:
%   - voxel_level:      Stored in .voxel_dat, voxels x images/observations
%   - node_level:       Stored in .node_dat, images/observations x nodes 
%                       There are 1 or more nodes per brain region
%                       Nodes are often patterns/linear combinations across voxels within regions 
%                       Region-to-node mappings are stored in a cell array,
%                       1 cell per region, in a [voxels x nodes] matrix
%                       within each cell. A private attribute (not
%                       accessible by users) maps regions to nodes in a
%                       [nodes x 1] integer vector, with every node coded
%                       as a unique integer. 
%   - region_level:     Stored in .region_dat, images/observations x regions 
%                       Stores averages across in-region voxels by default
%                       Mappings from voxels to regions are stored in region_atlas.dat [integer vector]
%   - network_level:    Stored in .network_dat, images/observations x networks 
%                       Networks are generalized 'pathways', each
%                       consisting of connections among k nodes.
%                       .network_dat stores the averages across in-pathway
%                       node-level data
%                       Mappings from voxels to regions are stored in region_atlas.dat [integer vector]
% 	- partition_level:  Stored in .partition_dat, images/observations x partitions. 
%                       Stores averages across in-partition voxels by default 
%
%   Hierarchical data structure:
%   Data within each field is stacked across participants/replicate blocks. 
%   For example, if we have t time points for each of 10 participants, the 
%   number of rows in .dat and columns in .node_dat would be t * 10.                    

% brainpathway Methods (a partial list; type doc brainpathway for more):
%   General:
% . 	descriptives        -  Get descriptives for an fmri_data or other image_vector object
%       enforce_variable_types	-  Re-casts variables in objects into standard data types, which can save
%       flip  


% - timeseries_data:        Level-1 (time series, t) cell array with one cell per subject, [t x k] data matrix in each cell, NaN is missing
% - person_data:            Level-2 (person-level, s) data for each region, [s x k] padded matrix, NaN is missing
% - data_properties:        Provenance for what has been done to data
%
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
%
% Some example tutorials:
% canlab_help_1_installing_tools
% canlab_help_2_load_a_sample_dataset
% canlab_help_3_voxelwise_t_test_walkthrough
% canlab_help_4_write_data_to_image_file_format
% canlab_help_5_regression_walkthrough


% Programmers' notes:
% Tor Wager, 7/25/2019 : initial creation

classdef brainpathway
    
    properties
        % also inherits the properties of image_vector.
        
        region_atlas (1, 1) atlas; % An atlas-class object defining k regions
        
        voxel_dat (:, :) single;   
        node_dat  (:, :) single; 
        region_dat (:, :) single; 
        network_dat (:, :) single; 
        partition_dat (:, :) single; 

        
        node_weights (1, :) cell;           % . fmri_data = []; % A series of n fmri_data objects, one per network, whose data field defines pattern weights
        
        connections (1, 1) struct = struct(''); % A series of matrices specifying [k x k] bivariate connections
        
        connections_apriori (:, :, :) logical = false(1, 1, 1);    % [1/0] logical matrices specifying existing connections, k x k x n for n networks
        
        % k x n latent variable weights - voxel weights. {1....k} cell k has {v x n}, v voxels x n networks
        % OR: weights could be in n fmri_data objects, one per network, with
        % - connections.est:        Estimated connectivity strengths
        % - connections.se:         Standard error of estimated connectivity strengths
        % - connections.metric      Metric type [r, cos_sim, tau, partial_r]
        % - graph_properties:
        % - timeseries_data:        Level-1 (time series, t) cell array with one cell per subject, [t x k] data matrix in each cell, NaN is missing
        % - person_data:            Level-2 (person-level, s) data for each region, [s x k] padded matrix, NaN is missing
        % - data_properties:        Provenance for what has been done to data
        %
        % - partitions:             An integer vector of partition labels for each
        %
        
        additional_info (1, 1) struct = struct(''); % A flexible structure defining user-specified additional information.
        
    end % properties
    
    methods
        
        % Class constructor
        function obj = brainpathway(varargin)
            
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS
            % -------------------------------------------------------------------------

            % Assign argument following ANY valid property name
            % allowable_inputs = properties(brainpathway)';       % should
            % be row cell vector of property names. This is recursive if we
            % do this though...
            allowable_inputs = {'region_atlas'    'voxel_dat'    'node_dat'    'region_dat'    'network_dat'    'partition_dat' 'node_weights'    'connections'    'connections_apriori'    'additional_info'};
            
            % optional inputs with default values - each keyword entered will create a variable of the same name
            
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        
                        case allowable_inputs
                            
                            eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                            
                        case 'noverbose'
                            % Other special inputs/control strings
                            
                        otherwise, warning(['Unknown input string option:' varargin{i}]);
                    end
                end
            end


%             input_atlas = false;
            
            for i = 1:length(varargin)
                
                switch class(varargin{i})
                    
                    case 'atlas', obj.region_atlas = varargin{i}; % input_atlas = true;
                        
                    case 'char'
                        switch varargin{i}
                            
                            case 'verbose', varargin{i} = []; % nothing else needed
                            case 'noverbose', verbose = 0; verbosestr = 'noverbose'; varargin{i} = [];
                            case 'sample2mask', sample2mask = 1; varargin{i} = [];
                                
                            case 'native_image_space' % do nothing, for convenience in calling scripts
                                
                            otherwise, warning(['Unknown input string option:' varargin{i}]);
                        end
                        
                end % switch class
                
            end % process varargin
            

            
            % -------------------------------------------------------------------------
            % SPECIAL COMMANDS/PROCESSES
            % -------------------------------------------------------------------------

            isatlas = cellfun(@(x) isa(x, 'atlas'), varargin);
            if ~any(isatlas) % load a default atlas if no atlas was passed in
                
                obj.region_atlas = load_atlas('canlab2018_2mm');
                
            end
            
            
            
             % initialize_nodes: Initialize nodes for each region, with weights of 1 if no other information is available
            if isempty(obj.node_weights)
                
                obj = intialize_nodes(obj);
                
            end
            
            % update_region_dat : Take voxel-level data and get region averages
            
            % update_node_dat : Take voxel-level data and get node activity 
            
            % need listeners: 
            % when vox data is assigned (if not empty), re-calculate region averages.  
            % when node weights are assigned, re-calculate node response
            % and connectivity
            
            % -------------------------------------------------------------------------
            % VALIDATE ATTRIBUTES OF INPUTS
            % -------------------------------------------------------------------------
            k = num_regions(obj.region_atlas);
            
            validateattributes(obj.region_atlas,{'atlas'},{},'brainpathway','.region_atlas');
            
            validateattributes(obj.node_weights,{'cell'},{'size', [1 k]},'brainpathway','.node_weights');
            
            % node weights must be voxels x nodes, so rows == region vox, for each region 
            
        end % class constructor
        
        
    end % methods
    
end % classdef

function obj = update_region_dat(obj)



end

function obj = intialize_nodes(obj)

k = num_regions(obj.region_atlas);

obj.node_weights = cell(1, k);

% Get number of voxels for each region

region_idx = obj.region_atlas.dat;
u = unique(region_idx);
nvox = zeros(1, k);

for i = 1:k
    nvox(i) = sum(region_idx == i);
end

% initialize
for i = 1:k
    
    obj.node_weights{i} = ones(nvox(i), 1);
    
end

end % function


% % Extracts local patterns from an image_vector
% % ---------------------------------------------------------
% pain_regions_pdm1 = extract_data(pain_regions, pdm1);
% % pain_regions_pdm1(1).all_data -> weights are stored in in all_data
% % save in .val field, which extract_data will use to extract
% for i = 1:length(pain_regions_pdm1)
%     pain_regions_pdm1(i).val = pain_regions_pdm1(i).all_data';
%     pain_regions_pdm1(i).Z = pain_regions_pdm1(i).all_data;
% end
% k = length(pain_regions_pdm1);


