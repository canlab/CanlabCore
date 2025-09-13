function obj = validate_object(obj)
% validate_object Validate the properties of an fmri_data object.
%
% :Usage:
% ::
%     obj.validate_object()
%
% :Description:
%     Checks that key properties of the fmri_data object are correctly sized and
%     of the correct type. In particular, the method validates:
%
%       - The number of rows in obj.dat equals (obj.volInfo.n_inmask - sum(obj.removed_voxels)).
%       - If obj.images_per_session is non-empty, its values sum to (number of columns in obj.dat - sum(obj.removed_images)).
%       - obj.removed_voxels is empty, a scalar zero, or a vector of ones/zeros (numeric or logical)
%         with length equal to obj.volInfo.n_inmask.
%       - obj.removed_images is empty, a scalar zero, or a vector of ones/zeros (numeric or logical).
%
%     Additionally, the volInfo structure is validated:
%
%       - volInfo.fname is a non-empty character array.
%       - volInfo.dim is a 1×3 integer vector.
%       - volInfo.dt is a 1×2 numeric integer vector.
%       - volInfo.pinfo is a 3×1 numeric vector.
%       - volInfo.mat is a 4×4 numeric matrix.
%       - volInfo.n is a 1×2 numeric integer vector.
%       - volInfo.descrip is a non-empty character array.
%       - volInfo.private is either empty or a 1×1 nifti object.
%       - volInfo.nvox is a scalar integer.
%       - volInfo.image_indx is a logical column vector with length equal to nvox.
%       - volInfo.wh_inmask is a numeric column vector with length equal to n_inmask.
%       - volInfo.n_inmask is a scalar integer.
%       - volInfo.xyzlist is a numeric matrix of size [n_inmask × 3].
%       - volInfo.cluster is a numeric column vector with length equal to nvox.
%
% :Example:
% ::
%     fmri_obj.validate_object();
%
% Author: Your Name
% Date: YYYY-MM-DD
% License: GNU General Public License v3 or later

% -------------------------------------------------------------------------
% Check specific sizes of variables that must add up
% -------------------------------------------------------------------------


%% Validate dimensions of obj.dat relative to volInfo and removed_voxels
% Compute expected number of voxels after removal.
if isempty(obj.removed_voxels)
    removedVoxels = 0;

elseif isscalar(obj.removed_voxels) && obj.removed_voxels == 0
    removedVoxels = 0;

else

    % Validate removed_voxels if not empty or scalar zero.
    validateattributes(obj.removed_voxels, {'numeric','logical'}, ...
        {'vector', 'numel', obj.volInfo.n_inmask}, mfilename, 'removed_voxels');
    removedVoxels = sum(obj.removed_voxels(:));

end

expectedVoxels = obj.volInfo.n_inmask - removedVoxels;
actualVoxels = size(obj.dat, 1);

if actualVoxels ~= expectedVoxels
    error('fmri_data:InvalidDimensions', ...
        'Number of rows in obj.dat (%d) must equal volInfo.n_inmask (%d) minus sum(removed_voxels) (%d).', ...
        actualVoxels, obj.volInfo.n_inmask, removedVoxels);
end

%% Validate images_per_session relative to obj.dat and removed_images
numImages = size(obj.dat, 2);

if isempty(obj.removed_images)
    removedImages = 0;
elseif isscalar(obj.removed_images) && obj.removed_images == 0
    removedImages = 0;
else
    validateattributes(obj.removed_images, {'numeric','logical'}, {'vector'}, mfilename, 'removed_images');
    removedImages = sum(obj.removed_images(:));
end

if ~isempty(obj.images_per_session)

    validateattributes(obj.images_per_session, {'numeric'}, {'vector'}, mfilename, 'images_per_session');

    if sum(obj.images_per_session) ~= (numImages - removedImages)
        error('fmri_data:InvalidDimensions', ...
            'Sum of images_per_session (%d) must equal (number of columns in obj.dat (%d) minus sum(removed_images) (%d)).', ...
            sum(obj.images_per_session), numImages, removedImages);
    end
end

% -------------------------------------------------------------------------
% Validate volInfo structure fields
% -------------------------------------------------------------------------

volInfo = obj.volInfo;

% fname: non-empty char array.
validateattributes(volInfo.fname, {'char'}, {'nonempty'}, mfilename, 'volInfo.fname');

% dim: 1x3 integer vector.
validateattributes(volInfo.dim, {'numeric'}, {'size',[1,3], 'integer'}, mfilename, 'volInfo.dim');

% dt: 1x2 numeric integer vector.
validateattributes(volInfo.dt, {'numeric'}, {'size',[1,2], 'integer'}, mfilename, 'volInfo.dt');

% pinfo: 3x1 numeric vector.
validateattributes(volInfo.pinfo, {'numeric'}, {'size',[3,1]}, mfilename, 'volInfo.pinfo');

% mat: 4x4 numeric matrix.
validateattributes(volInfo.mat, {'numeric'}, {'size',[4,4]}, mfilename, 'volInfo.mat');

% n: 1x2 numeric integer vector.
validateattributes(volInfo.n, {'numeric'}, {'size',[1,2], 'integer'}, mfilename, 'volInfo.n');

% descrip: non-empty char array.
validateattributes(volInfo.descrip, {'char'}, {'nonempty'}, mfilename, 'volInfo.descrip');

% private: empty or 1x1 nifti class object.
if ~isempty(volInfo.private)
    validateattributes(volInfo.private, {'nifti'}, {'scalar'}, mfilename, 'volInfo.private');
end

% nvox: scalar integer.
validateattributes(volInfo.nvox, {'numeric'}, {'scalar','integer'}, mfilename, 'volInfo.nvox');

% image_indx: logical column vector, length must equal nvox.
validateattributes(volInfo.image_indx, {'logical'}, {'column'}, mfilename, 'volInfo.image_indx');
if numel(volInfo.image_indx) ~= volInfo.nvox
    error('fmri_data:InvalidDimensions', ...
        'volInfo.image_indx must have length equal to volInfo.nvox (%d).', volInfo.nvox);
end

% wh_inmask: numeric column vector, length must equal n_inmask.
validateattributes(volInfo.wh_inmask, {'numeric'}, {'column'}, mfilename, 'volInfo.wh_inmask');
if numel(volInfo.wh_inmask) ~= volInfo.n_inmask
    error('fmri_data:InvalidDimensions', ...
        'volInfo.wh_inmask must have length equal to volInfo.n_inmask (%d).', volInfo.n_inmask);
end

% n_inmask: scalar integer.
validateattributes(volInfo.n_inmask, {'numeric'}, {'scalar','integer'}, mfilename, 'volInfo.n_inmask');

% xyzlist: numeric matrix with size [n_inmask x 3].
validateattributes(volInfo.xyzlist, {'numeric'}, {'size',[volInfo.n_inmask,3]}, mfilename, 'volInfo.xyzlist');

% cluster: numeric column vector, length must equal nvox.
validateattributes(volInfo.cluster, {'numeric'}, {'column'}, mfilename, 'volInfo.cluster');

if numel(volInfo.cluster) ~= volInfo.n_inmask
    error('fmri_data:InvalidDimensions', ...
        'volInfo.cluster must have length equal to volInfo.nvox (%d).', volInfo.nvox);
end



% -------------------------------------------------------------------------
% Validate the fields in image_metadata
% -------------------------------------------------------------------------

% For logical flag fields: they must be a logical scalar if not NaN.
logicalFields = {'is_timeseries', 'is_single_trial_series', 'is_first_level_maps', ...
    'is_MNI_space', 'is_HP_filtered', 'covariates_removed'};

for i = 1:length(logicalFields)
    val = obj.image_metadata.(logicalFields{i});
    % If not NaN, check that the value is logical.
    if ~isnan(val)
        validateattributes(val, {'logical'}, {'scalar'}, mfilename, logicalFields{i});
    end
end

% For numeric fields: they must be a numeric scalar (NaN is acceptable).
numericFields = {'TR_in_sec', 'HP_filter_cutoff_sec'};
for i = 1:length(numericFields)
    val = obj.image_metadata.(numericFields{i});
    validateattributes(val, {'numeric'}, {'scalar'}, mfilename, numericFields{i});
end



%% If all validations pass, print a confirmation (if verbose mode is enabled)
if isprop(obj, 'verbose') && obj.verbose
    fprintf('fmri_data object validated successfully.\n');
end



end % main function


