function obj = enforce_variable_types(obj)
% Re-casts variables in objects into standard data types, which can save space
% in memory.  Also removes empty voxels as a space-saving device.
%
% obj = enforce_variable_types(obj)
%
% Needs testing to ensure that doing this does not break any other
% display/computation functions.

obj = remove_empty(obj);
obj.dat = single(obj.dat);

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
    
    obj.p = single(obj.p);
    
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

