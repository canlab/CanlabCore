function obj = enforce_variable_types(obj)
% enforce_variable_types Re-cast image_vector fields into standard data types to save space.
%
% Re-casts variables in objects into standard data types, which can save
% space in memory. Also removes empty voxels as a space-saving device.
%
% Needs testing to ensure that doing this does not break any other
% display/computation functions.
%
% :Usage:
% ::
%
%     obj = enforce_variable_types(obj)
%
% :Inputs:
%
%   **obj:**
%        An image_vector or subclass (fmri_data, statistic_image, atlas)
%        object.
%
% :Outputs:
%
%   **obj:**
%        The input object with empty voxels removed and the following
%        type conversions applied:
%
%        - .dat -> single (or double for statistic_image with type 'T')
%        - .volInfo.wh_inmask -> uint32
%        - .volInfo.xyzlist -> uint16
%        - .volInfo.cluster -> uint32
%        - For fmri_data, .mask is similarly down-cast.
%        - For statistic_image, .sig -> logical and .p -> single
%          (except for 'T' type, which retains higher precision .p).
%
%        If the object lacks an 'image_metadata' property, a default
%        image_metadata struct is initialized.
%
% :Examples:
% ::
%
%     dat = enforce_variable_types(dat);
%
% :See also:
%   - remove_empty
%   - replace_empty
%
% ..
%    Programmers' notes:
%    Edited 02/26/2026 (Zizhuang Miao): force the .dat field of a
%    statistic_image of type 'T' to be double, to enable more accurate
%    estimates of p values.
% ..

obj = remove_empty(obj);
% EDITED: force the .dat field of a statistic_image to be double
% to enbale more accurate estimates of p values
% (Zizhuang Miao, 02/26/2026)
if isa(obj, 'statistic_image') && strcmp(obj.type, 'T')
    obj.dat = double(obj.dat);
else
    obj.dat = single(obj.dat);
end

obj.volInfo.wh_inmask = uint32(obj.volInfo.wh_inmask);  % 4 billion positive valued integers
obj.volInfo.xyzlist = uint16(obj.volInfo.xyzlist);
obj.volInfo.cluster = uint32(obj.volInfo.cluster);

if isa(obj, 'fmri_data')
    % future: remove .mask altogether...?
    obj.mask.dat = single(obj.mask.dat);
    
    obj.mask = remove_empty(obj.mask);
    
    obj.mask.volInfo.wh_inmask = uint32(obj.mask.volInfo.wh_inmask);  % 4 billion positive valued integers
    obj.mask.volInfo.xyzlist = uint16(obj.mask.volInfo.xyzlist);
    obj.mask.volInfo.cluster = uint32(obj.mask.volInfo.cluster);
    
end

if isa(obj, 'statistic_image')
    
    obj.sig = logical(obj.sig);
    
    if ~strcmp(obj.type, 'T')
        obj.p = single(obj.p);
    end
    
end

if ~isprop(obj, 'image_metadata')

    % Define a structure "image_metadata" with default values for the fmri_data object.
    image_metadata = struct( ...
        'is_timeseries',           NaN, ...  % Logical flag: true, false, or NaN (unknown)
        'is_single_trial_series',  NaN, ...  % Logical flag: true, false, or NaN (unknown)
        'is_first_level_maps',     NaN, ...  % Logical flag: true, false, or NaN (unknown)
        'is_MNI_space',            NaN, ...  % Logical flag: true, false, or NaN (unknown)
        'is_HP_filtered',          NaN, ...  % Logical flag: true, false, or NaN (unknown)
        'covariates_removed',      NaN, ...  % Logical flag: true, false, or NaN (unknown)
        'TR_in_sec',               NaN, ...  % Numeric value (real number in seconds) or NaN if not set
        'HP_filter_cutoff_sec',    NaN);     % Numeric value (real number in seconds) or NaN if not set

end


end % main function

