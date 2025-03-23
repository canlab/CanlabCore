% Define validation functions for each object property

% X: should be numeric (empty allowed)
validate_X = @(x) validateattributes(x, {'numeric'}, {});

% mask: must be a 1×1 fmri_mask_image object
validate_mask = @(x) isa(x, 'fmri_mask_image');

% mask_descrip: must be a char string
validate_mask_descrip = @(x) ischar(x);

% images_per_session: should be numeric (empty allowed)
validate_images_per_session = @(x) validateattributes(x, {'numeric'}, {});

% Y: should be numeric (empty allowed)
validate_Y = @(x) validateattributes(x, {'numeric'}, {});

% Y_names: must be a char array
validate_Y_names = @(x) ischar(x);

% Y_descrip: must be a char string
validate_Y_descrip = @(x) ischar(x);

% covariates: should be numeric (empty allowed)
validate_covariates = @(x) validateattributes(x, {'numeric'}, {});

% covariate_names: must be a cell array of char strings
validate_covariate_names = @(x) iscell(x) && all(cellfun(@ischar, x));

% covariates_descrip: must be a char string
validate_covariates_descrip = @(x) ischar(x);

% history_descrip: must be a char string
validate_history_descrip = @(x) ischar(x);

% additional_info: must be a structure (empty struct allowed)
validate_additional_info = @(x) isstruct(x);

% metadata_table: must be a table
validate_metadata_table = @(x) isa(x, 'table');

% source_notes: must be a char string
validate_source_notes = @(x) ischar(x);

% dat: should be numeric (e.g., single, double; empty allowed)
validate_dat = @(x) validateattributes(x, {'numeric'}, {});

% dat_descrip: should be a char string or empty
validate_dat_descrip = @(x) ischar(x) || isempty(x);

% volInfo: must be a struct or empty
validate_volInfo = @(x) (isstruct(x) || isempty(x));

% removed_voxels: should be numeric (scalar or vector)
validate_removed_voxels = @(x) validateattributes(x, {'numeric'}, {});

% removed_images: should be numeric (scalar or vector)
validate_removed_images = @(x) validateattributes(x, {'numeric'}, {});

% image_names: must be a char array (not a cell array)
validate_image_names = @(x) ischar(x);

% fullpath: must be a char array (not a cell array)
validate_fullpath = @(x) ischar(x);

% files_exist: must be a logical array
validate_files_exist = @(x) islogical(x);

% history: must be a cell array of char strings
validate_history = @(x) iscell(x) && all(cellfun(@ischar, x));

prop_names = {'X', 'mask', 'mask_descrip', 'images_per_session', 'Y', 'Y_names', 'Y_descrip', 'covariates', 'covariate_names', 'covariates_descrip', 'history_descrip', 'additional_info', 'metadata_table', 'source_notes', 'dat', 'dat_descrip', 'volInfo', 'removed_voxels', 'removed_images', 'image_names', 'fullpath', 'files_exist', 'history'};


p = inputParser;
p.addParameter('X', [], validate_X);
p.addParameter('mask', [], validate_mask);
p.addParameter('mask_descrip', '', validate_mask_descrip);
p.addParameter('images_per_session', [], validate_images_per_session);
p.addParameter('Y', [], validate_Y);
p.addParameter('Y_names', '', validate_Y_names);
p.addParameter('Y_descrip', '', validate_Y_descrip);
p.addParameter('covariates', [], validate_covariates);
p.addParameter('covariate_names', {''}, validate_covariate_names);
p.addParameter('covariates_descrip', '', validate_covariates_descrip);
p.addParameter('history_descrip', '', validate_history_descrip);
p.addParameter('additional_info', struct(), validate_additional_info);
p.addParameter('metadata_table', table(), validate_metadata_table);
p.addParameter('source_notes', '', validate_source_notes);
p.addParameter('dat', [], validate_dat);
p.addParameter('dat_descrip', '', validate_dat_descrip);
p.addParameter('volInfo', struct(), validate_volInfo);
p.addParameter('removed_voxels', 0, validate_removed_voxels);
p.addParameter('removed_images', 0, validate_removed_images);
p.addParameter('image_names', '', validate_image_names);
p.addParameter('fullpath', '', validate_fullpath);
p.addParameter('files_exist', false, validate_files_exist);
p.addParameter('history', {''}, validate_history);

p.parse(varargin{:});
ARGS = p.Results;

