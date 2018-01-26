% atlas: Subclass of image_vector that stores information about brain atlases and parcellations
%
% 'atlas' is a data class containing information about brain atlases and parcellations
% stored in a structure-like object.  It inherits the properties and
% methods of fmri_data and image_vector objects.
%
% Creating an atlas object requires either images with probability maps (a 4-d image) 
% Or an integer-valued image with one integer per atlas region. 
% For full functionality, the atlas also requires both probability maps and 
% text labels, one per region, in a cell array. But some functionality will
% work without these.
%
% Basic usage for creating a new atlas image:
% obj = atlas(image_names, [maskinput], [other optional inputs]) 
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
% obj = atlas
% - obj is the object.
% - It will be created with a standard brain mask, brainmask.nii
% - This image should be placed on your Matlab path
% - The space information is stored in obj.volInfo
% - Data is stored in obj.dat, in a [voxels x images] matrix
% - You can replace or append data to the obj.dat field.
%
% You can create an atlas object with extacted image data.
% - Let "imgs" be a string array or cell array of image names
% - This command creates an object with your (4-D) image data:
% - fmri_dat = atlas(imgs);
% - Only values in the standard brain mask, brainmask.nii, will be included.
% - This saves memory by reducing the number of voxels saved dramatically.
%
% You can specify any mask you'd like to extract data from.
% - Let "maskimagename" be a string array with a mask image name.
% - this command creates the object with data saved in the mask:
% - fmri_dat = atlas(imgs, maskimagename);
% - The mask information is saved in fmri_dat.mask
%
% e.g., this extracts data from images within the standard brain mask:
% dat = atlas(imgs, which('brainmask.nii'));
%
% Properties and methods
% -----------------------------------------------------------------------
% Properties are data fields associated with an object.
% Try typing the name of an object (class instance) you create to see its
% properties, and a link to its methods (things you can run specifically
% with this object type). For example: After creating an atlas object 
% called fmri_dat, as above, type fmri_dat to see its properties.
%
% There are many other methods that you can apply to atlas objects to
% do different things.
% - Try typing methods(atlas) for a list.
% - You always pass in an atlas object as the first argument.
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
% The atlas object has a number of fields for appending specific types of data.
%
% - You can replace or append data to the fmri_dat.dat field.
% - Attach custom descriptions in several fields to document your object.
% - The "history" field stores a cell array of strings with the processing
% history of the object. Some methods add to this history automatically.
%
%
% Examples:
% parcellation_file = 'CIT168toMNI152_prob_atlas_bilat_1mm.nii';  %'prob_atlas_bilateral.nii';
% labels = {'Put' 'Cau' 'NAC' 'BST_SLEA' 'GPe' 'GPi' 'SNc' 'RN' 'SNr' 'PBP' 'VTA' 'VeP' 'Haben' 'Hythal' 'Mamm_Nuc' 'STN'};
% dat = atlas(which(parcellation_file), 'labels', labels, ...
%               'space_description', 'MNI152 space', ...
%               'references', 'Pauli 2018 Bioarxiv: CIT168 from Human Connectome Project data', 'noverbose');
% % Display:
% orthviews(dat);
% figure; montage(dat);
%
% % Convert to region object and display:
% r = atlas2region(dat);
% orthviews(r)
% montage(r);

% Programmers' notes:
% Tor Wager, 1/14/18 : Created from hybrid of image_vector (parent) and some fmri_data code

classdef atlas < image_vector
    
    properties
        % also inherits the properties of image_vector.
        
        atlas_name              % a short description or name of the atlas
        probability_maps        % voxels x regions matrix with probability values for each region
        labels
        label_descriptions      % a regions x 1 cell array of long-form descriptions for labels
        labels_2
        labels_3
        labels_4
        labels_5
        references
        space_description = '';
        property_descriptions = { ...
            'atlas_name: a short description or name of the atlas' ...
            'probability_maps: voxels x regions matrix with probability values for each region' ...
            'dat: Integer vector of voxels x 1, with region index numbers' ...
            'labels: Cell array of text strings for atlas labels, one per region/index number' ...
            'label_descriptions: a regions x 1 cell array of long-form descriptions for labels' ...
            'labels_2, etc: Cell array of text strings for atlas labels, one per region/index number' ...
            'references: String matrix of associated publications' ...
            'space_description: Description of atlas space/template, e.g., MNI152 space, original origin' ...
            };
        
        additional_info = struct('');
        
    end % properties
    
    methods
        
        % Class constructor
        function obj = atlas(image_names, varargin)
            %
            % [obj, cl_with_averages] = atlas(image_names, varargin)
            %
            % Reads a set of image files, and returns an atlas object with stored data.
            
            % ---------------------------------
            % Create empty object, and return if no additional
            % arguments
            % ---------------------------------
            
            obj.probability_maps = single([]);
            obj.dat = int32(0);
            [obj.labels obj.labels_2 obj.labels_3 obj.labels_4 obj.labels_5] = deal({''});
            obj.history = {''};
            obj.additional_info = struct('');
            
            % DEFAULT INPUTS
            % -----------------------------------
            
            maskinput = [];
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
                            
                        case 'mask', maskinput = varargin{i+1}; varargin{i+1} = [];
                                
                        case 'atlas_name', obj.atlas_name = varargin{i+1}; varargin{i+1} = [];
                        case 'labels', obj.labels = varargin{i+1}; varargin{i+1} = [];
                            
                        case 'label_descriptions', obj.label_descriptions = varargin{i+1}; varargin{i+1} = [];

                        case 'labels_2', obj.labels_2 = varargin{i+1}; varargin{i+1} = [];
                        case 'labels_3', obj.labels_3 = varargin{i+1}; varargin{i+1} = [];
                        case 'labels_4', obj.labels_4 = varargin{i+1}; varargin{i+1} = [];
                        case 'labels_5', obj.labels_5 = varargin{i+1}; varargin{i+1} = [];
                            
                        case 'references', obj.references = varargin{i+1}; varargin{i+1} = [];
                        case 'space_description', obj.space_description = varargin{i+1}; varargin{i+1} = [];
                                
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
               error('fmridata: image_names is empty. Invalid input image names.');
               
            end
            
            if iscell(image_names), image_names = char(image_names{:}); end
            
            % ---------------------------------
            % Special: if existing image_vector
            % map into space of atlas and return
            % ---------------------------------
            
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
        
            % Empty mask: use default
            if isempty(maskinput)   
                
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
            %obj = create(obj, 'mask', maskobj, verbosestr);
            
            obj = create(obj, 'image_names', image_names, verbosestr);
            
            % append data, depending on what it looks like
            % --------------------------------------------------------------------------------
            is_prob_image = size(imgdat, 2) >= 1 && all(imgdat(:) >=0) && all(imgdat(:) <=1);
            is_parcel_index = size(imgdat, 2) == 1 && all(abs(imgdat(:) - round(imgdat(:))) < 100 * eps);
            
            if is_prob_image
                
                obj.probability_maps = imgdat;
                
                obj = probability_maps_to_region_index(obj);
                
            elseif is_parcel_index
                
                obj.dat = int32(round(imgdat));
                
            else
                error('I don''t recognize this type of image. Enter probability maps, component weights, or parcel index');
            
            end
            
            clear imgdat
            
            % append description info
            % --------------------------------------------------------------------------------
            [dd, ff, ee] = fileparts(maskobj.volInfo.fname);
            mask_image_name = [ff ee];
            %obj = create(obj, 'dat_descrip', sprintf('Data from %s: %s', mask_image_name, obj.dat_descrip));
            %obj = create(obj, 'mask_descrip', mask_image_name, verbosestr);
            
            obj.volInfo = maskobj.volInfo;
            
            obj.history(end+1) = {sprintf('Sampled to space of %s', maskobj.space_defining_image_name)};
            obj.history(end+1) = {['Masked with ' mask_image_name]};
            
            obj = check_image_filenames(obj, verbosestr);
            
            % Check properties and enforce variable types
            % --------------------------------------------------------------------------------  
            obj = check_properties(obj);
            
        end % constructor function
        
    end % methods
    
    
end

