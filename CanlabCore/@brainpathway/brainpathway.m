% brainpathway: Data class for storing and analyzing brain connectivity and pathways
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
% Brainpathway is a handle class, which means...
% It includes listeners, which automatically recalculate derivative data
% values if you attach new data. For example, if you update voxel_dat,
% region averages, node data, connectivity, etc. will all be recalculated.
%
% A brainpathway object has these parts:
% - region_atlas:                An atlas-class object defining k regions
% - WEIGHTS object ->  n fmri_data objects, one per network,
% - connectivity:            A series of connectivity matrices specifying
% [k x k] bivariate connections for various types of data stored in the object
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
% brainpathway Properties (a partial list; type doc brainpathway for more):
%   xxx                     - xxxxxxxx
%
% fmri_data Methods (a partial list; type doc fmri_data for more):
%   	xxx                 - xxxx
%
%   Data extraction:
%   	xxx                 - xxxx
%
%
%   Handling brain space and image selection:
%   	xxx                 - xxxx
%
%
%   Display and visualization:
%   	xxx                 - xxxx
%
%
%   Data processing and analysis:
%   	xxx                 - xxxx
%
%
% -------------------------------------------------------------------------
% Examples and help:
% -------------------------------------------------------------------------
%
% To list properties and methods for this object, type:
% doc fmri_data, methods(fmri_data)
%
% b = brainpathway(pain_pathways);  % Construct a brainpathway object from an atlas object, here "pain_pathways"
% b.region_atlas = pain_pathways;   % Alternate way of assigning a region atlas, or changing the atlas
% b.voxel_dat = randn(352328, 20);
% b.voxel_dat = [];                 % Remove original data, leave calculated derivatives intact
%
% b = brainpathway(pain_pathways); % Construct a brainpathway object from an atlas object
% b.region_dat = ST_cleaned.big_regions;
% plot_connectivity(b, 'notext')
%
% Add clusters:
% b = cluster_regions(b);
% plot_connectivity(b, 'notext', 'partitions', b.node_clusters, 'partition_labels', {'Thal', 'Bstem/Amy' 'S2/insula', 'S1');
%
% Could do (but not needed - will be done automatically when voxel data or node weights are assigned): 
% b.update_node_data(b);
% b.update_node_connectivity(b);
%
% b.node_dat = ST_cleaned.cPDM;
% plot_connectivity(b, 'notext')
% 
% *** to-do : add external variables (e.g., temp, pain)
% *** partitions

% Programmers' notes:
% Tor Wager, 7/25/2019 : initial creation

classdef brainpathway < handle
    
    properties (SetObservable = true)
        
        region_atlas (1, 1) atlas;          % An atlas-class object defining k regions
        
        voxel_dat (:, :) single;            % A [voxels x images/observations] matrix of data
        node_dat  (:, :) single;
        region_dat (:, :) single;
        network_dat (:, :) single;
        partition_dat (:, :) single;
        
        node_weights (1, :) cell;           %  A series of n cells, one per node. Each cell contains a vector of pattern weights across voxels
        node_labels (1, :) cell;           %  A series of n cells, one per node. Each cell contains a char array name for the node.
        node_clusters (1, :) int32;         % n integers indicating cluster membership (see cluster_regions)
        
        region_indx_for_nodes (1, :) int32 = []; 
        
        connectivity (1, 1) struct = struct('regions', [], 'nodes', []); % A series of matrices specifying [k x k] bivariate connections
        
        %       Specify a function handle and optional arguments to the
        %       function (in addition to data). This allows connectivity_properties to be defined in a very flexible way, using multiple functions and inputs.
        %       For example, the default is:
        %       obj.connectivity_properties = struct('c_fun_han', @corr, 'c_fun_arguments', {})
        %           ...which uses Pearson's correlations
        %       An alternative using Spearman's correlations (rank) would be:
        %       obj.connectivity_properties = struct('c_fun_han', @corr, 'c_fun_arguments', {'type', 'Spearman'})
        %       Or, for partial correlations:
        %       obj.connectivity_properties = struct('c_fun_han', @partialcorr, 'c_fun_arguments', {'type', 'Spearman'})
        
        connectivity_properties (1, 1) struct = struct('c_fun_han', @corr, 'c_fun_arguments', {{}});  % 'metric', 'r', 'rank', false, 'robust', false, 'partialcorr', false);
        
        connections_apriori (:, :, :) logical = logical([]);    % [1/0] logical matrices specifying existing connections, k x k x n for n networks
        
        % k x n latent variable weights - voxel weights. {1....k} cell k has {v x n}, v voxels x n networks
        % OR: weights could be in n fmri_data objects, one per network, with
        % - connectivity.est:        Estimated connectivity strengths
        % - connectivity.se:         Standard error of estimated connectivity strengths
        % - connectivity.metric      Metric type [r, cos_sim, tau, partial_r]
        % - graph_properties:
        % - timeseries_data:        Level-1 (time series, t) cell array with one cell per subject, [t x k] data matrix in each cell, NaN is missing
        % - person_data:            Level-2 (person-level, s) data for each region, [s x k] padded matrix, NaN is missing
        % - data_properties:        Provenance for what has been done to data
        %
        % - partitions:             An integer vector of partition labels for each
        %
        
        additional_info (1, 1) struct = struct(''); % A flexible structure defining user-specified additional information.
        
        listeners = [];         % listeners
        
        verbose = true;
        
    end % properties
    
    events
        % Set events for listeners here
        
    end % events
    
    
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
            allowable_inputs = {'region_atlas'    'voxel_dat'    'node_dat'    'region_dat'    'network_dat'    'partition_dat' 'node_weights'    'connectivity'    'connections_apriori'    'additional_info'};
            
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
            
            % -------------------------------------------------------------------------
            % LISTENERS: Check properties and recalculate/update data fields
            % -------------------------------------------------------------------------
            

            
            % When voxel_dat is set/updated, ...
            % ------------------------------------------------------------
            % update region_dat
            obj.listeners = addlistener(obj,'voxel_dat', 'PostSet',  @(src, evt) brainpathway.update_region_data(obj, src, evt));

            % update node_dat
            obj.listeners(end+1) = addlistener(obj,'voxel_dat', 'PostSet', @(src, evt) brainpathway.update_node_data(obj, src, evt));
            
            
            % When region_dat is set/updated...
            % ------------------------------------------------------------
            % update region connectivity
            obj.listeners(end+1) = addlistener(obj,'region_dat', 'PostSet',  @(src, evt) brainpathway.update_region_connectivity(obj, src, evt));
  
            % When node_weights are set/updated...
            % ------------------------------------------------------------
            % update node_dat
            obj.listeners(end+1) = addlistener(obj,'node_weights', 'PostSet', @(src, evt) brainpathway.update_node_data(obj, src, evt));
            
            
            % When node_dat is set/updated...
            % ------------------------------------------------------------
            % update node connectivity
            obj.listeners(end+1) = addlistener(obj,'node_dat', 'PostSet', @(src, evt) brainpathway.update_node_connectivity(obj, src, evt));
            
            % When connectivity_properties are set/updated...
            % ------------------------------------------------------------
            % update node connectivity
            obj.listeners(end+1) = addlistener(obj,'connectivity_properties', 'PostSet', @(src, evt) brainpathway.update_node_connectivity(obj, src, evt));
            
            % update region connectivity
            obj.listeners(end+1) = addlistener(obj,'connectivity_properties', 'PostSet', @(src, evt) brainpathway.update_region_connectivity(obj, src, evt));
            
            % ------------------------------------------------------------
            % When region_atlas is set/updated ....
            % update region_dat
            obj.listeners(end+1) = addlistener(obj,'region_atlas', 'PostSet', @(src, evt) brainpathway.update_region_data(obj, src, evt));
            
            % initialize nodes
            %obj.listeners(end+1) = addlistener(obj,'region_atlas', 'PostSet', @(src, evt) brainpathway.intialize_nodes(obj, src, evt)); % this should update nodes too...not yet...
                        
            
        end % class constructor
        
        
        
        function obj = intialize_nodes(obj)
            % Initialize nodes with 1 node per region
            % (It is possible to assign multiple nodes per region)
            
            k = num_regions(obj.region_atlas);
            
            obj.node_weights = cell(1, k);
            
            % Get number of voxels for each region
            nvox = count_vox_per_region(obj);
            
            %obj = brainpathway.update_node_data(obj);
            
            
            % initialize
            for i = 1:k
                
                obj.node_weights{i} = ones(nvox(i), 1) ./ nvox(i);
                
                obj.node_labels{i} = obj.region_atlas.labels{i};
                
            end
            
            obj.region_indx_for_nodes = get_node_info(obj);
            
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
        
        
        function obj = cluster_regions(obj, varargin)
            
            n_clusters = max(8, num_regions(obj.region_atlas));
            
            if length(varargin) > 0
                n_clusters = varargin{1};
            end
            
            obj.node_clusters = (clusterdata(obj.region_dat', 'linkage', 'ward', 'maxclust', n_clusters))';
            
        end
        
        function plot_connectivity(obj, varargin)
        % Takes any optional input arguments to plot_correlation_matrix
        
            S = struct('r', obj.connectivity.regions.r, 'p', obj.connectivity.regions.p, 'sig', obj.connectivity.regions.p < 0.05);
            
            Xlabels = format_strings_for_legend(obj.region_atlas.labels);
                    
            figtitle = 'brainpathway_connectivity_view';
            
            if isempty(obj.connectivity.nodes)
                create_figure(figtitle);
            else
                create_figure(figtitle, 1, 2);
            end
            
            OUT = plot_correlation_matrix(S, 'nofigure', ...
                'var_names', Xlabels, varargin{:});

            num_nodes = length(obj.region_indx_for_nodes);
            if num_nodes == num_regions(obj.region_atlas)
                
                node_labels = Xlabels;
                
            else node_labels = [];
            end
            
            if isempty(obj.connectivity.nodes)
                return
            end
            
            subplot(1, 2, 2);
            
            S = struct('r', obj.connectivity.nodes.r, 'p', obj.connectivity.nodes.p, 'sig', obj.connectivity.nodes.p < 0.05);
            
            OUT = plot_correlation_matrix(S, 'nofigure', ...
                'var_names', node_labels, varargin{:});
            
            
            %             Xpartitions = ones(size(S.r, 2), 1);
%             partitionlabels = {'Regions'};
    

%             OUT = plot_correlation_matrix(S, 'nofigure', ...
%                 'var_names', Xlabels, 'partitions', Xpartitions, 'partitionlabels', partitionlabels);

        end
        
    end % methods
    
    methods (Static)
        
        function obj = update_region_data(obj, src, evt)
            % Update region average data by extracting from obj.voxel_dat
            
            % Notes: this will not do anything fancy with zero voxels or
            % NaNs. Input is validated to be non-NaN.
            % An alternative would be atlas.apply_parcellation
            
            v = size(obj.region_atlas.dat, 1); % num vox
            k = num_regions(obj.region_atlas);
            
            if isempty(obj.voxel_dat)
                
                % Special case: leave existing region_dat alone
                %obj.region_dat = [];
                return
            end
            
            if obj.verbose, fprintf('Updating region averages.\n'); end
            
            validateattributes(obj.voxel_dat,{'numeric'},{'nrows', v, '2d', 'nonnan'}, 'brainpathway.update_region_data','.voxel_dat');
            
            % nvox = count_vox_per_region(obj);
            
            parcel_indic = condf2indic(obj.region_atlas.dat, 'integers', k);
            
            %for computing means, scale each column of parcels to sum to 1
            mydat = bsxfun(@rdivide, parcel_indic, nansum(parcel_indic));
            
            %matrix products will give us the mean now...
            obj.region_dat = obj.voxel_dat' * mydat;
            
        end
        
        
        function obj = update_node_data(obj, src, evt)
            % Update node response data by applying node weights for each region
            
            % Notes: this will not do anything fancy with zero voxels or
            % NaNs. Input is validated to be non-NaN.
            % An alternative would be atlas.apply_parcellation
            
            %             1 cell per region, in a [voxels x nodes] matrix
            % %                       within each cell. A private attribute (not
            % %                       accessible by users) maps regions to nodes in a
            % %                       [nodes x 1] integer vector, with every node coded
            % %                       as a unique integer.
            
            %simfun = @dot;

            v = size(obj.region_atlas.dat, 1); % num vox
            k = num_regions(obj.region_atlas);
            
            if isempty(obj.voxel_dat)
                % Special case: leave existing node_dat alone
                
                return
            end
            
            if obj.verbose, fprintf('Updating node response data.\n'); end
            
            validateattributes(obj.voxel_dat,{'numeric'},{'nrows', v, '2d', 'nonnan'}, 'brainpathway.update_region_data','.voxel_dat');
            validateattributes(obj.node_weights,{'cell'},{}, 'brainpathway.update_node_data','.node_weights');
            
            % nvox = count_vox_per_region(obj);
            
            [obj.region_indx_for_nodes, ~, node_start, node_end] = get_node_info(obj);
            
            % index into voxels in region
            parcel_indic = logical(condf2indic(obj.region_atlas.dat, 'integers', k));
            
            for i = 1:k
                
                mydat = obj.voxel_dat(parcel_indic(:, i), :);
                
                myvals = (obj.node_weights{i}' * mydat)' ;  % obs x nodes for this region
                
                allnodedat(:, [node_start(i):node_end(i)]) = myvals;
                
            end
            
            obj.node_dat = allnodedat;
        end
        
        
        function obj = update_region_connectivity(obj, src, evt)
            
            
            dat = obj.region_dat;
            outputfield = 'regions';
            
            if obj.verbose, fprintf('Updating obj.connectivity.%s.\n', outputfield); end
            
            fhan = obj.connectivity_properties.c_fun_han;
            in_args = obj.connectivity_properties.c_fun_arguments;
            
            if isempty(dat)
                [obj.connectivity.(outputfield).r,  obj.connectivity.(outputfield).p] = deal([]);
              
                return
            end
            
            % note: will work for corr, partialcorr now, not other methods,
            % e.g., multilevel methods like ttest3d.
            [r, p] = fhan(dat, in_args{:});
            
            obj.connectivity.(outputfield).r = r;
            obj.connectivity.(outputfield).p = p;
            
        end
        
        
        function obj = update_node_connectivity(obj, src, evt)
            
            dat = obj.node_dat;
            outputfield = 'nodes';
            
            if obj.verbose, fprintf('Updating obj.connectivity.%s.\n', outputfield); end
            
            fhan = obj.connectivity_properties.c_fun_han;
            in_args = obj.connectivity_properties.c_fun_arguments;
            
            if isempty(dat)
                [obj.connectivity.(outputfield).r,  obj.connectivity.(outputfield).p] = deal([]);
                
                return
            end
            
            % note: will work for corr, partialcorr now, not other methods,
            % e.g., multilevel methods like ttest3d.
            [r, p] = fhan(dat, in_args{:});
            
            obj.connectivity.(outputfield).r = r;
            obj.connectivity.(outputfield).p = p;
            
        end
        
        
    end % static methods
    
end % classdef


% Get number of voxels for each region
function nvox = count_vox_per_region(obj)

k = num_regions(obj.region_atlas);

region_idx = obj.region_atlas.dat;
u = unique(region_idx);
nvox = zeros(1, k);

for i = 1:k
    nvox(i) = sum(region_idx == i);
end

end


function [region_indx_for_nodes, sz, node_start, node_end] = get_node_info(obj)

% Count nodes per region
sz = cellfun(@size, obj.node_weights, 'UniformOutput', false);
sz = cat(1, sz{:});
sz = sz(:, 2);

% Index into starting and ending nodes for each region
node_start = cumsum(sz);
node_end = [node_start(2:end) - 1 ; sum(sz)];

k = num_regions(obj.region_atlas);
region_indx_for_nodes = cell(1, k);

for i = 1:k
    region_indx_for_nodes{i} = i * ones(1, sz(i));
end

region_indx_for_nodes = cat(2, region_indx_for_nodes{:});

end

