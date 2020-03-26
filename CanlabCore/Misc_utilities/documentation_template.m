% Standard text for function documentation
%
% First line: One-line summary description of function
%
% :Usage:
% ::
%
%     [list outputs here] = function_name(list inputs here, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) <year>  <name of author>
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
%
% :Inputs:
%
%   **param1:**
%        description of param1
%
%   **param2:**
%        description of param2
%
% :Optional Inputs:
%   **param1:**
%        description of param1
%
%   **param2:**
%        description of param2
%
% :Outputs:
%
%   **out1:**
%        description of out1
%
%   **out2:**
%        description of out2
%
% :Examples:
% ::
%
%    % give examples of code here
%    param1 = abc();
%    param2 = xyz();
%    [out1,out2] = func_call(param1, param2)
%
% :References:
%   CITATION(s) HERE
%
% :See also:
%   - list other functions related to this one, and alternatives*
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% BELOW IS A STANDARD TEMPLATE FOR DEFINING VARIABLE (OPTIONAL) INPUT
% ARGUMENTS. MANY FUNCTIONS NEED TO PARSE OPTIONAL ARGS, SO THIS MAY BE
% USEFUL. I (TOR) RECOMMEND COPYING AND PASTING THE DOCUMENTATION ABOVE AND
% ALL THE STUFF THROUGH ARGUMENT VALIDATION, BELOW, AND MODIFYING IT FOR
% YOUR FUNCTION.

% SEE BELOW FOR TEMPLATE FOR OBJECT CONSTRUCTORS. THIS HAS STRUCTURED
% DOCUMENTATION WITH SOME TEXT FORMATTING RULES THAT WILL MAKE IT
% COMPATIBLE WITH MATLAB'S DOC COMMAND, WHICH WILL THEN CREATE AN HTML PAGE
% LISTING METHODS AND ATTRIBUTES. I (TOR) RECOMMEND COPYING AND PASTING THE 
% DOCUMENTATION BELOW AND MODIFYING IT FOR YOUR FUNCTION.

% ..
%    DEFAULTS AND INPUTS
% ..

% Some useful things
n_cols = 80;
sep_str = repmat('_', 1, n_cols);  % see textwrap


% ----------------------------------------------------------------------
% Parse inputs
% ----------------------------------------------------------------------
% This 2019 version uses the inputParser object. Older schemes are below.
% Note: With this, you can pass in EITHER keyword, value pairs OR a
% structure with fields that define keywords.
% Note: You can also make this a subfunction at the end of your function,
% e.g., function ARGS = parse_inputs(varargin)...end, and call it using:

ARGS = parse_inputs(varargin{:});

% If you want to distribute arguments back out to variables, use this:

IN = p.Results;
fn = fieldnames(IN);

for i = 1:length(fn)
    str = sprintf('%s = IN.(''%s'');', fn{i}, fn{i});
    eval(str)
end

% Logical flags
% ----------------------------------------------------------------------
if any(strcmp(varargin, 'noverbose')), doverbose = false; end
if any(strcmp(varargin, 'noplots')), doplots = false; end

% ======== This goes at the end of the file / after the main function ========
function ARGS = parse_inputs(varargin)

p = inputParser;

% Validation functions - customized for each type of input
% ----------------------------------------------------------------------

valfcn_scalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar'});

valfcn_number = @(x) validateattributes(x, {'numeric'}, {'nonempty'}); % scalar or vector

% Validation: Region object, structure, or [x1 x2 x3] triplet 
valfcn_custom = @(x) isstruct(x) || isa(x, 'region') || (~isempty(x) && all(size(x) - [1 3] == 0) && all(isnumeric(x)));

% Validation: [x1 x2 x3] triplet 
valfcn_xyz = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'size', [1 3]});

valfcn_logical = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar', '>=', 0, '<=', 1}); % could enter numeric 0,1 or logical


% Required inputs 
% ----------------------------------------------------------------------
p.addRequired('x', valfcn_custom);
p.addRequired('y', valfcn_custom);

% Optional inputs 
% ----------------------------------------------------------------------
% Pattern: keyword, value, validation function handle

p.addParameter('color', [.9 .2 0], valfcn_xyz);
p.addParameter('bendpercent', .15, valfcn_number); % can be scalar or vector
p.addParameter('thickness', .1, valfcn_scalar);
p.addParameter('nstreamlines', 30, valfcn_scalar);

% Parse inputs and distribute out to variable names in workspace
% ----------------------------------------------------------------------
% e.g., p.parse([30 1 0], [-40 0 10], 'bendpercent', .1);
p.parse(varargin{:});

ARGS = p.Results;

end % parse_inputs



% ==================== end pattern =====================================


% Next is the previously used scheme for CANlab tools for several years


% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------

rowsz = [];
doplot = 0;
basistype = 'spm+disp';
% initalize optional variables to default values here.

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

% This is a compact way to assign multiple variables. The input argument
% names and variable names must match, however:

allowable_inputs = {'var_names' 'p_thr' 'dospearman' 'dopartial' 'dofdr' 'dofigure' 'doimage' 'docircles' 'dotext' 'colorlimit' 'text_x_offset' 'text_y_offset' 'text_fsize' 'text_nonsig_color' 'text_sig_color' 'partitions' 'partitioncolors' 'partitionlabels'};

% optional inputs with default values - each keyword entered will create a variable of the same name

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case allowable_inputs
                
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% This pattern will flexibly assign arguments based on keywords. 
% The names of the input keyword and variable created do not need to match.
% Multiple input keywords can be mapped to the same variable.

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            case {'rows', 'rowsz'}, rowsz = varargin{i+1}; varargin{i+1} = [];
            case 'plot', doplot = 1; 
            case 'basistype', basistype = varargin{i+1}; varargin{i+1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% VALIDATE ATTRIBUTES OF INPUTS
% -------------------------------------------------------------------------

validateattributes(X,{'numeric'},{'2d'},'plot_correlation_matrix','X', 1);

logical_args = {'dofdr' 'false' 'dospearman' 'dofigure' 'doimage' 'docircles' 'dotext'};
for i = 1:length(logical_args)
    
    my_arg = eval([logical_args{i} ';']);
    validateattributes(my_arg,{'logical'},{'scalar'},'plot_correlation_matrix',logical_args{i});

end

cell_args = {'partitioncolors' 'partitionlabels' 'var_names'};
for i = 1:length(cell_args)
    
    my_arg = eval([cell_args{i} ';']);
    validateattributes(my_arg,{'cell'},{},'plot_correlation_matrix',cell_args{i});

end

vector_args = {'text_sig_color' 'text_nonsig_color'};
for i = 1:length(vector_args)
    
    my_arg = eval([vector_args{i} ';']);
    validateattributes(my_arg,{'numeric' 'vector' 'nonnegative' 'nonnan'},{},'plot_correlation_matrix',vector_args{i});

end

supportedImageClasses    = {'int8','uint8','int16','uint16','int32','uint32','single','double','logical'};
supportedImageAttributes = {'real','nonsparse','nonempty'};
validateattributes(V,supportedImageClasses,supportedImageAttributes,mfilename,'V');


% SEE:
%  validateattributes Check validity of array.
%     validateattributes(A,CLASSES,ATTRIBUTES) validates that array A belongs
%     to at least one of the specified CLASSES and has all of the specified
%     ATTRIBUTES. 
%
% See also inputParser.m


% ==================== end pattern =====================================



% THIS IS ANOTHER WAY THAT IS MORE STREAMLINED, BUT ACTUALLY LESS INTUITIVE

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

% Defaults
% -----------------------------------
% Set color maps for + / - values
poscm = colormap_tor([1 0 .5], [1 1 0], [.9 .6 .1]);  %reddish-purple to orange to yellow
negcm = colormap_tor([0 0 1], [0 1 1], [.5 0 1]);  % cyan to purple to dark blue

% optional inputs with default values
% -----------------------------------
% - allowable_args is a cell array of argument names
% - avoid spaces, special characters, and names of existing functions
% - variables will be assigned based on these names
%   i.e., if you use an arg named 'cl', a variable called cl will be
%   created in the workspace

allowable_args = {'cl', 'ycut_mm', 'pos_colormap', 'neg_colormap', ...
    'surface_handles', 'existingfig'};

default_values = {[], [], poscm, negcm, ...
    [], 0};

% define actions for each input
% -----------------------------------
% - cell array with one cell for each allowable argument
% - these have special meanings in the code below
% - allowable actions for inputs in the code below are: 'assign_next_input' or 'flag_on'

actions = {'assign_next_input', 'assign_next_input', 'assign_next_input', 'assign_next_input', ...
    'assign_next_input', 'flag_on'};

% logical vector and indices of which inputs are text
textargs = cellfun(@ischar, varargin);
whtextargs = find(textargs);

for i = 1:length(allowable_args)
    
    % assign default
    % -------------------------------------------------------------------------
    
    eval([allowable_args{i} ' =  default_values{i};']);
    
    wh = strcmp(allowable_args{i}, varargin(textargs));
    
    if any(wh)
        % Optional argument has been entered
        % -------------------------------------------------------------------------
        
        wh = whtextargs(wh);
        if length(wh) > 1, warning(['input ' allowable_args{i} ' is duplicated.']); end
        
        switch actions{i}
            case 'assign_next_input'
                eval([allowable_args{i} ' = varargin{wh(1) + 1};']);
                
            case 'flag_on'
                eval([allowable_args{i} ' = 1;']);
                
            otherwise
                error(['Coding bug: Illegal action for argument ' allowable_args{i}])
        end
        
    end % argument is input
end

% END DEFAULTS AND INPUTS
% -------------------------------------------------------------------------


% Below are some default headers

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% Sub-functions
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% OBJECT FUNCTION DEFINITION - followed by one-line description
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

