function [fwd, mean_fwd, est_outliers, mvmt_mtx_mm] = framewise_displacement(mvmt_mtx, varargin)
% FRAMEWISE_DISPLACEMENT Calculate framewise displacement from 6 movement parameters.
%
% This function computes the framewise displacement (FD) based on Power et al. (2012, 2019)
% to estimate movement outliers in fMRI data. Framewise displacement quantifies the
% amount of head movement between successive frames using both translation and rotation parameters.
%
% CAUTION!! Rotations must be 1st 3 columns of mvmt_mtx, then translations
%
% Usage:
%   mvmt_file = 'path_to_BIDS_confounds_timeseries.tsv';
%   mvmt_data = importBIDSfile(mvmt_file);  % Import data from fMRIprep-style confounds_timeseries.tsv
%   mvmt_mtx = [mvmt_data.rot_x, mvmt_data.rot_y, mvmt_data.rot_z, mvmt_data.trans_x, mvmt_data.trans_y, mvmt_data.trans_z]; 
%   [mvmt_mtx, corr_out, uncorr_out] = framewise_displacement(mvmt_mtx);
%
% Inputs:
%   mvmt_mtx - A T x 6 matrix containing rotation (rot_x, rot_y, rot_z) in radians
%              and translation (trans_x, trans_y, trans_z) in mm.
%
% Outputs:
%   fwd - T x 1 vector of overall framewise displacement estimates
%   mean_fwd - mean framewise displacement 
%   est_outliers - Logical vector of estimated outliers after correction
%   mvmt_mtx_mm - The input matrix with rotation parameters converted to mm.
%
% References:
%   - Power et al. (2012), "Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion."
%   - Power et al. (2019), "Respiration pseudomotion in functional connectivity MRI."
%
% :Optional Inputs:
%
%   **'thresh'**
%   Threshold in mm for determining outliers; Default = 0.5 mm
%
% Author: Michael Sun, Ph.D. 10/1/2024

% -------------------------------------------------------------------------
% DEFAULT ARGUMENT VALUES
% -------------------------------------------------------------------------
thresh = 0.5; % mm 

% -------------------------------------------------------------------------
% OPTIONAL INPUTS
% -------------------------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'thresh'
                eval([varargin{i} ' = varargin{i+1}; varargin{i+1} = [];']);

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Number of time points (T is the number of rows in mvmt_mtx)
T = size(mvmt_mtx, 1);  

% Set radius for rotation adjustment (in mm)
radius = 50;  

% Extract translation and rotation parameters (assuming mvmt_mtx is in the order: rot_x, rot_y, rot_z, trans_x, trans_y, trans_z)
mvmt_mtx_mm = mvmt_mtx(:, 1:6);

% Adjust rotational displacements (columns 1-3) by converting radians to displacement in mm (arc length on a sphere of radius 50 mm)
% Formula: arc length = (2 * pi / 360) * radius * rotation (in degrees)
mvmt_mtx_mm(:, 1:3) = (2 * radius * pi / 360) * mvmt_mtx_mm(:, 1:3);

% Calculate differences between successive time points (frame-wise differences)
dts = diff(mvmt_mtx_mm);  % This computes T-1 differences

% Compute framewise displacement (FD) as the sum of absolute differences between consecutive frames
% Add a zero at the beginning since there is no displacement for the first time point
fwd = [0; sum(abs(dts), 2)];

% Initialize the output logical vectors to mark outliers based on a threshold (> 0.25 mm)
est_outliers = false(T, 1);

% Mark outliers based on threshold
est_outliers(fwd > thresh) = true;
mean_fwd = nanmean(fwd);

end
