function [obj, true_obj, noise_obj] = sim_data(obj, varargin)
% Create simulated data and attach it to an fmri_data object (obj)
%
% :Usage:
% ::
%
%     obj = sim_data(obj, [optional inputs])
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2016 Tor Wager
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
%   **obj:**
%       an fmri_data object
%
% :Optional Inputs:
%
% **'n':**
%       Number of images/observations (e.g., subjects)
%
% **'d'**
%       Cohen's d for true regions
%
% **'true_region_mask'**
%       Followed by mask of which regions to add true signal
%       Default is the 'frontoparietal' network from Yeo et al. 2011
%
% **'noise_sigma'**
%       Followed by noise sigma (std) value
%
% **'smoothness'**
%       Followed by smoothing value (default = 5)
%
%
% :Outputs:
%
%   **obj:**
%        Object with simulated data, signal + noise
%        obj.additional_info{1} = structure with simulation PARAMS
%        obj.additional_info{2} = structure with SIGNAL calculations and indicators;
%
%
%   **true_obj:**
%        Signal data object
%
%   **noise_obj:**
%        Noise data object
%
% :Examples:
% ::
%
% Create a simulated dataset with default values:
% [obj, true_obj, noise_obj] = sim_data(fmri_data);
%
% Create one with default values but 30 subjects:
%[obj, true_obj, noise_obj] = sim_data(fmri_data, 'n', 30);
%
% Display the mean true signal:
% [obj, true_obj, noise_obj] = sim_data(fmri_data, 'n', 30);
% o2 = canlab_results_fmridisplay([], 'compact2');
% o2 = addblobs(o2, region(mean(true_obj)), 'color', [0 0 1]);
%
% Show some summary information about the mean effects in true-signal and
% false-signal regions:
% obj.additional_info{2}
%
% Generate some data with true mean shift and true correlation with a
% simulated predictor, and show summary stats:
%
% N = 20;
% [obj, true_obj, noise_obj] = sim_data(fmri_data, 'n', N, 'r', .5);
% obj.additional_info{2}
%
% Generate some data no mean signal but a true correlation with a
% simulated predictor, and show summary stats:
%
% [obj, true_obj, noise_obj] = sim_data(fmri_data, 'n', N, 'd', 0, 'r', .5);
% obj.additional_info{2}
%
% Generate null data with no true signal, and show summary stats:
% Note: there is still a default mask with 'true' and 'false' regions, but
% effect sizes are zero.
% 
% [obj, true_obj, noise_obj] = sim_data(fmri_data, 'n', N, 'd', 0);
% OR:
% [obj, true_obj, noise_obj] = sim_data(fmri_data, 'n', N, 'null');
%
% Re-calculate the voxel-wise correlation with the outcome obj.Y:
% Uses a stored Anonymous Function handle to get vector-with-matrix correlation
% corr_matrix = obj.additional_info{1}.corr_matrix; 
% r = corr_matrix(obj.Y, obj.dat');
%
% :References:
%   None.
%
% :See also:
%   -
%
%
% ..
%    Programmers' notes:
%   Documentation and update, July 2016, Tor Wager
%   Update for simulations, Feb 2019
% ..

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

rng('shuffle');         % seeds the random number generator 

[v, n] = size(obj.dat); % voxels x imagess

true_region_mask = [];  % to use default mask
d = 0.5;                % Cohen's d for true regions (true effect size)
r = 0;                  % Correlation with outcome vector; individual diffs

noise_sigma = 1;        % Noise standard deviation (Gaussian)
smoothness = 5;         % Data smoothing

doplot = 0;

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'n', n = varargin{i+1}; varargin{i+1} = [];
            case 'd', d = varargin{i+1}; varargin{i+1} = [];
  
            case {'correlation', 'r'}, r = varargin{i+1}; varargin{i+1} = [];
                    
            case 'true_region_mask', true_region_mask = varargin{i+1}; varargin{i+1} = [];
            case 'noise_sigma', noise_sigma = varargin{i+1}; varargin{i+1} = [];
            case 'smoothness', smoothness = varargin{i+1}; varargin{i+1} = [];
                
            case 'null', d = 0; r = 0;
                
            case 'plot', doplot = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% Build a structure with hyperparameters to save in sim objects
% -------------------------------------------------------------------------

PARAMS = struct('v', v, 'n', n, 'd', d, 'r', r, 'noise_sigma', noise_sigma, 'smoothness', smoothness);

% function for calculating voxel-wise correlation with predictor
corr_matrix = @(a, b) ((a-mean(a))' * (b-mean(b)) ./ (length(a) - 1))' ./ (std(a)*std(b)');
PARAMS.corr_matrix = corr_matrix;


% Create noise object
% -------------------------------------------------------------------------

% noise
obj.dat = noise_sigma .* randn(v, n);

% smooth noise
% and since we have smoothed away noise var, renormalize noise.
if smoothness
    obj = preprocess(obj, 'smooth', smoothness);
    obj = rescale(obj, 'zscoreimages');
    obj.dat = noise_sigma .* obj.dat;
end

% Predictor (regressor[s])
Y = randn(n, 1);
obj.X = Y;

% Outcome (for multivariate prediction)
obj.Y = Y;

% Add PARAMS
obj.additional_info{1} = PARAMS;

noise_obj = obj;

% Signal(H1) = mean activity + indiv diffs + error
% Signal X   = d + r * Y + e
% d = Cohen's d for mean shift
% r = Pearson's r
% Y = "true" signal of std = 1, N(0, 1)

% Signal: True data
% ----------------------------------------------------------
if ~isa(true_region_mask, 'image_vector') % do not use isempty(true_region_mask) because fmri_data.isempty checks for 0 data values
    
    true_obj = load_image_set('bucknerlab');
    true_obj = get_wh_image(true_obj, 6);  % 'Frontoparietal'
    
else
    true_obj = fmri_data(true_region_mask);
end

true_obj = resample_space(true_obj, obj);

wh = logical(true_obj.dat(:, 1));

true_obj.dat = zeros(v, n);
true_obj.dat(wh, :) = d .* noise_sigma; % Add in mean signal.  
% Multiply by noise_sigma so that effect size d is constant for any value of noise_sigma entered. 

% Add in variance with predictor
% This will add error variance if not modeled in analysis
% But it is not strictly error, as it could be accounted for by the known
% predictor.  Thus, entering r ~= 0 will reduce the mean effect size (d) if
% the predictor is not regressed out.
true_obj.dat(wh, :) = true_obj.dat(wh, :) + repmat(r .* obj.X', sum(wh), 1); 
true_obj.Y = obj.Y;

% Final data: noise + signal
% ----------------------------------------------------------
obj.dat = obj.dat + true_obj.dat;


% Add PARAMS and helpful data vectors with true signal

d_obs = (mean(obj.dat') ./ std(obj.dat'))';  
r_obs = corr_matrix(obj.Y, obj.dat');

mean_d_true = mean(d_obs(wh));
mean_d_false = mean(d_obs(~wh));
mean_r_true = mean(r_obs(wh));
mean_r_false = mean(r_obs(~wh));

descrip = 'wh_true = indicator for true effect voxels. d_obs=voxelwise d. r_obs = voxel_wise r with obj.Y';

SIGNAL = struct('wh_true', wh, 'd_obs', d_obs, 'r_obs', r_obs, ...
    'mean_d_true', mean_d_true, 'mean_d_false', mean_d_false, 'mean_r_true', mean_r_true, 'mean_r_false', mean_r_false);

SIGNAL.descrip = descrip;

obj.additional_info{1} = PARAMS;
obj.additional_info{2} = SIGNAL;


% Optional plot
% ----------------------------------------------------------
if doplot
    
    o2 = canlab_results_fmridisplay([], 'multirow', 3);
    
    o2 = addblobs(o2, region(mean(true_obj)), 'color', [0 0 1], 'wh_montages', 5:6);
    o2 = addblobs(o2, region(get_wh_image(noise_obj, 1)), 'mincolor', [.2 .2 .2], 'maxcolor', [.7 .7 .7], 'wh_montages', 3:4);
    
    o2 = addblobs(o2, region(get_wh_image(obj, 1)), 'mincolor', [.2 .2 .2], 'maxcolor', [.7 .7 .7], 'wh_montages', 1:2, 'trans');
    o2 = addblobs(o2, region(mean(true_obj)), 'color', [0 0 1], 'wh_montages', 1:2, 'outline', 'linewidth', 3);
    
end


end % function


