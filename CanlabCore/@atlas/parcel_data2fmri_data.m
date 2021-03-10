function obj = parcel_data2fmri_data(atlas_obj, data_values)
% Utility that transforms parcel-wise data values into an fmri_data object in voxelwise image space
% Expands parcel values so that they are replicated for each voxel in the parcel.
%
% :Usage:
% ::
%
% obj = parcel_data2fmri_data(atlas_obj, data_values)
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2021 Tor Wager
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
%   **atlas_obj:**
%        An atlas-class object with k distinct regions or parcels
%
%   **data_values:**
%        n x k matrix of k data values for each of n parcels. Each of k
%        will be saved as a different image in the output obj
%
% :Optional Inputs:
%
% :Outputs:
%
%   **obj:**
%        An fmri_data object containing voxelwise data. .dat contains
%        the values. Parcel values will be replicated for each voxel in the
%        parcel.
%
%
% :Examples:
% ::
%
% atl = load_atlas('canlab2018_2mm');
% n = length(atl.labels); % number of parcels
% data_values = randn(n, 10); % create 10 images with random data, 1 unique value per parcel
% test_obj = parcel_data2fmri_data(atl, data_values);
%
% :References:
%   None 
%
% :See also:
%   parcel_stats2statistic_image
%

% Create placeholder statistic image

atlas_obj = replace_empty(atlas_obj);   
k = size(data_values, 2);
placeholder_vec = ones(atlas_obj.volInfo.n_inmask, k);
all_parcel_idx = double(atlas_obj.dat);
u = unique(all_parcel_idx); u(u == 0) = [];

% initialize output object with correct size

obj = fmri_data( ...
    struct('dat', single(0 .* placeholder_vec), ...
    'volInfo', atlas_obj.volInfo, ...
    'removed_voxels', atlas_obj.removed_voxels) ...
    );

% number of parcels 

n = size(data_values, 1); % num parcels

if n ~= length(u)
        error('Parcel indices in atlas_obj and number of rows in parcel-wise data_values do not match. Check inputs.')
end



obj.dat = zeros(size(placeholder_vec, 1), k); % voxels x conditions


    for j = 1:length(u)
        % For each parcel, fill in statistic_image object
        % ----------------------------------------------------
        parcelidx = u(j);
        
        wh_vox = all_parcel_idx == parcelidx;
        
        % map parcels to voxels
        obj.dat(wh_vox, :) = repmat(data_values(j, :), sum(wh_vox), 1);
                
    end
       
obj = enforce_variable_types(obj);  % space-saving: 5/24/17

% input_struct.t_statistic_obj = obj;

end % function
