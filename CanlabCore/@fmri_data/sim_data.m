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
%   **dat:**
%        Object with simulated data, signal + noise
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
% ..

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------

rng('shuffle');         % seeds the random number generator 

[v, n] = size(obj.dat); % voxels x imagess

true_region_mask = [];  % to use default mask
d = 0.5;                % Cohen's d for true regions (true effect size)

noise_sigma = 1;        % Noise standard deviation (Gaussian)
smoothness = 5;         % Data smoothing

doplot = 0;

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'n', n = varargin{i+1}; varargin{i+1} = [];
            case 'd', d = varargin{i+1}; varargin{i+1} = [];
            case 'true_region_mask', true_region_mask = varargin{i+1}; varargin{i+1} = [];
            case 'noise_sigma', noise_sigma = varargin{i+1}; varargin{i+1} = [];
            case 'smoothness', smoothness = varargin{i+1}; varargin{i+1} = [];
                
            case 'plot', doplot = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


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
obj.X = randn(n, 1);

% Outcome (for multivariate prediction)
% Now all noise; could easily add signal
obj.Y = randn(n, 1);

noise_obj = obj;

% Signal: True data
% ----------------------------------------------------------
if isempty(true_region_mask)
    
    true_obj = load_image_set('bucknerlab');
    true_obj = get_wh_image(true_obj, 6);  % 'Frontoparietal'
    
else
    true_obj = fmri_data(true_region_mask);
end

true_obj = resample_space(true_obj, obj);

wh = logical(true_obj.dat(:, 1));

true_obj.dat = zeros(v, n);
true_obj.dat(wh, :) = d .* noise_sigma;

% Final data: noise + signal
% ----------------------------------------------------------
obj.dat = obj.dat + true_obj.dat;

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


