function [mvmt_matrix mvmt_regs_24, mvmt_table] = movement_regressors(mvmt_file)
% Read a realignment parameter file and create a set of movement regressors with quadratic/derivative transformations
%
% :Usage:
% ::
%     [mvmt_matrix, mvmt_regs_24, mvmt_table] = movement_regressors(mvmt_file)
%
% :Inputs:
%
%   **mvmt_file:**  
%        String. Full path to the movement parameter file (e.g., rp_sub-sid001567_task-pinel_acq-s1p2_run-03_bold.txt)
%        generated after realignment (motion correction). This file contains the head motion 
%        estimates (translations and rotations) required to align each image to a reference image.
%
% :Outputs:
%
%   **mvmt_matrix:**  
%        Numeric matrix. The raw movement parameters imported from mvmt_file.
%
%   **mvmt_regs_24:**  
%        Numeric matrix. A matrix containing 24 movement-related regressors for each image.
%        These regressors include the standardized movement parameters, their squared values,
%        first differences (derivatives), and squared first differences.
%
%   **mvmt_table:**  
%        Table object. A table version of mvmt_regs_24 with one variable per regressor.
%        The variable names are as follows:
%            - 'mvmt_x', 'mvmt_y', 'mvmt_z', 'mvmt_r', 'mvmt_p', 'mvmt_y'
%            - 'mvmt_x^2', 'mvmt_y^2', 'mvmt_z^2', 'mvmt_r^2', 'mvmt_p^2', 'mvmt_y^2'
%            - 'mvmt_xdiff', 'mvmt_ydiff', 'mvmt_zdiff', 'mvmt_rdiff', 'mvmt_pdiff', 'mvmt_ydiff'
%            - 'mvmt_xdiff^2', 'mvmt_ydiff^2', 'mvmt_zdiff^2', 'mvmt_rdiff^2', 'mvmt_pdiff^2', 'mvmt_ydiff^2'
%
% :Examples:
% ::
%     % Example: Compute movement regressors from a realignment parameters file.
%     mvmt_file = fullfile(pwd, 'rp_sub-sid001567_task-pinel_acq-s1p2_run-03_bold.txt');
%     [mvmt_matrix, mvmt_regs_24, mvmt_table] = movement_regressors(mvmt_file);
%     % Display the table of movement regressors.
%     disp(mvmt_table);
%
% :References:
%   Power, J.D. et al. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion.
%   Power, J.D. (2019). [Details regarding the correction of respiratory pseudomotion artifacts].
%
% :See also:
%   framewise_displacement, outliers, spm_realign, realignment
%

% More notes:
% Power (2012) framewise displacement: Uses the average absolute deviation across 6 movement parameters.
% Power 2012: "Spurious but systematic correlations in functional connectivity MRInetworks arisefrom subject motion"
% "FDi = |Δdix| + |Δdiy| + |Δdiz| + |Δαi| + |Δβi| + |Δγi|, where Δdix = d(i − 1)x − dix, 
% and similarly for the other rigid body parameters [dix diy diz αi βi γi]. 
% Rotational displacements were converted from degrees to millimeters by calculating 
% displace- ment on the surface of a sphere of radius 50 mm, which is approximately 
% the mean distance from the cerebral cortex to the center of the head."  Also: FD 
% computed with respect to the volume ~2 seconds previous, with respiration pseudo-motion 
% filtered out of head motion traces (Power 2019). 
% Power 2019: 
% "Breathing can cause the head to move, but changes in lung volume can also cause a particular 
% kind of artifact called pseudomotion, which manifests as a shift of the brain when the lung 
% expands (Brosch et al., 2002; Durand et al., 2001; Raj et al., 2001)."
%
% Yoni's protocol for OLP:
% spike criteria: FD > .25mm OR abs(zscore(DVARS)) > 3, and the following 4 volumes (following ~2 sec)
% bandpass filter, passband = [.01 .1] Hz

mvmt_matrix = importdata(mvmt_file);

% fill out 24 movement-related parameters per run
% ----------------------------------------------------------------
mvmt = zscore(mvmt_matrix);
mvmt2 = mvmt .^ 2;
mvmt3 = [zeros(1, 6); diff(mvmt)];
mvmt4 = zscore(mvmt3) .^ 2;

% mvmt2 = zscore(mvmt) .^ 2;                  % Square the Z-scored data
% % mvmt3 = [zeros(1, 6); diff(zscore(mvmt))];  % Append zeroes as the first row
% mvmt3 = [zeros(1, size(mvmt, 2)); diff(zscore(mvmt))];  
% mvmt4 = zscore(mvmt3) .^ 2;                 % 

mvmt_regs_24 = [zscore(mvmt) mvmt2 mvmt3 mvmt4];

names = {'mvmt_x' 'mvmt_y' 'mvmt_z' 'mvmt_roll' 'mvmt_pitch' 'mvmt_yaw' };
names = [names 'mvmt_x^2' 'mvmt_y^2' 'mvmt_z^2' 'mvmt_roll^2' 'mvmt_pitch^2' 'mvmt_yaw^2'];
names = [names 'mvmt_xdiff' 'mvmt_ydiff' 'mvmt_zdiff' 'mvmt_rolldiff' 'mvmt_pitchdiff' 'mvmt_yawdiff'];
names = [names 'mvmt_xdiff^2' 'mvmt_ydiff^2' 'mvmt_zdiff^2' 'mvmt_rolldiff^2' 'mvmt_pitchdiff^2' 'mvmt_yawdiff^2'];

% Convert into a table using tasknames as variable names
mvmt_table = array2table(mvmt_regs_24, 'VariableNames', names);

% % simple geometric mean across translation and rotation; not absolute movement, which varies by voxel
% % ----------------------------------------------------------------
% geo_displacement = [0; sum(diff(mvmt) .^ 2, 2) .^ .5];
% 
% % adjust outliers: add geo displacement
% 
% high_movement_timepoints(geo_displacement > 0.25) = true;

%% New Code

% T = size(mvmt, 1);  % Number of time points
% % len = length(names);  % Number of files
% 
% % FDM = zeros(T, len);  % Preallocate forward displacement matrix
% radius = 50;  % Radius for rotation adjustment (in mm)
% 
% % for i = 1:len
%     % Read each file (assumed to contain motion parameters)
%     % ts = readmatrix(names{i});
% 
% 
%     % This requires a file with X motion correction parameters from *_LR*.txt
%     % Select the first 6 columns (translation and rotation) xyz-rpyaw
%     ts = ts(:, 1:6);
% 
%     % Adjust rotation columns (4-6) using the formula for arc length
%     ts(:, 4:6) = (2 * radius * pi / 360) * ts(:, 4:6);
% 
%     % Calculate the difference between consecutive time points
%     dts = diff(ts);
% 
%     % Compute forward displacement as the sum of absolute differences
%     fwd = sum(abs(dts), 2);
% 
%     % Store the forward displacement in FDM, ensuring it fits in the matrix
%     if length(fwd) == (T - 1)
%         FDM(2:T, i) = fwd;  % Start from the second row, as fwd is one element shorter
%     end
% % end
% 
% % Display the FDM matrix
% disp(FDM);




end
