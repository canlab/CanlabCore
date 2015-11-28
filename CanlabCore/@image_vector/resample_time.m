function obj = resample_time(obj, source_TR, target_TR, varargin)
% Resample the time-series images (source_time_interval) in an fmri_data object (obj) 
% to the different time series (target_time_interval). Works for all image_vector objects.
%
%   - obj = resample_time(obj, source_time_interval, target_time_interval, varargin)
%
% :Optional Inputs:
%
%   **meth (Interpolation methods):**
%        You can enter resampling method as optional input. Takes any input to
%          - 'nearest' - nearest neighbor interpolation
%          - 'linear'  - linear interpolation (default)
%          - 'spline'  - spline interpolation
%          - 'cubic'   - cubic interpolation as long as the data is uniformly
%                       spaced, otherwise the same as 'spline'
%   **slice:**
%      A fraction of the slice timing correction.
%      The default is 0.5, meaning if your TR is 2s, the time point of your TR image
%      will be considered as the middle point of the TR bins. You can use this option
%      to use different time points. If you are upsampling your data (i.e.,
%      your target TR is shorter than your source TR), you need to discard the
%      first column of your data. This function will return the first time point data as NaN. 
%
%
% :Examples:
% ::
%
%    dat = fmri_data('/Volumes/RAID1/labdata/current/BMRK3/Imaging/spatiotemp_biomarker/STmarker1.img');
%    dat = resample_time(dat, 2, 1.3) 
%
%    % with options:
%    dat = resample_time(dat, 2, 1.3, 'meth', 'linear', 'slice', .3)
%

slice_frac = .5;
intp_meth = [];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'meth', 'methods'}
                intp_meth = varargin{i+1};
            case {'slice'}
                slice_frac = varargin{i+1};
            otherwise
        end
    end
end

n_ts = size(obj.dat, 2);
n_vox = size(obj.dat,1);
obj_out = obj;
obj_out.dat = [];

sc_start_time = source_TR.*slice_frac;
tg_start_time = target_TR.*slice_frac;

end_time = source_TR.*n_ts-sc_start_time;

[X, Y] = meshgrid(sc_start_time:source_TR:end_time, 1:n_vox);
[Xq, Yq] = meshgrid(tg_start_time:target_TR:end_time, 1:n_vox);

if isempty(intp_meth)
    obj_out.dat = interp2(X, Y, obj.dat, Xq, Yq);
else
    obj_out.dat = interp2(X, Y, obj.dat, Xq, Yq, intp_meth);
end

% obj_out.dat(:,(sum(isnan(obj_out.dat)) ~= 0)) = [];

obj = obj_out;

end
