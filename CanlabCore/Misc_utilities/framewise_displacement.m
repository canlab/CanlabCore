function [mvmt_mtx, est_outliers_corr, est_outliers_uncorr] = framewise_displacement(mvmt_mtx)
% FRAMEWISE_DISPLACEMENT Calculate framewise displacement from 6 movement parameters.
%
% This function computes the framewise displacement (FD) based on Power et al. (2012, 2019)
% to estimate movement outliers in fMRI data. Framewise displacement quantifies the
% amount of head movement between successive frames using both translation and rotation parameters.
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
%   mvmt_mtx - The input matrix with rotation parameters converted to mm.
%   est_outliers_corr - Logical vector of estimated outliers after correction (> 0.25 mm).
%   est_outliers_uncorr - Logical vector of estimated outliers before correction (> 0.25 mm).
%
% References:
%   - Power et al. (2012), "Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion."
%   - Power et al. (2019), "Respiration pseudomotion in functional connectivity MRI."
%
% Author: Michael Sun, Ph.D. 10/1/2024


% Number of time points (T is the number of rows in mvmt_mtx)
T = size(mvmt_mtx, 1);  

% Set radius for rotation adjustment (in mm)
radius = 50;  

% Extract translation and rotation parameters (assuming mvmt_mtx is in the order: rot_x, rot_y, rot_z, trans_x, trans_y, trans_z)
ts = mvmt_mtx(:, 1:6);

% Adjust rotational displacements (columns 1-3) by converting radians to displacement in mm (arc length on a sphere of radius 50 mm)
% Formula: arc length = (2 * pi / 360) * radius * rotation (in degrees)
ts(:, 1:3) = (2 * radius * pi / 360) * ts(:, 1:3);

% Calculate differences between successive time points (frame-wise differences)
dts = diff(ts);  % This computes T-1 differences

% Compute framewise displacement (FD) as the sum of absolute differences between consecutive frames
% Add a zero at the beginning since there is no displacement for the first time point
fwd = [0; sum(abs(dts), 2)];

% Initialize the output logical vectors to mark outliers based on a threshold (> 0.25 mm)
est_outliers_corr = false(T, 1);
est_outliers_uncorr = false(T, 1);

% Mark outliers based on threshold
est_outliers_corr(fwd > 0.25) = true;
est_outliers_uncorr(fwd > 0.25) = true;

end