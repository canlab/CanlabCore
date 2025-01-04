function t_obj = parcel_stats2statistic_image(atlas_obj, tscores, pvalues, dfe, sig)
% Utility that transforms parcel-wise t-statistics and p-values into a statistic_image object in voxelwise image space
% Expands parcel values so that they are replicated for each voxel in the parcel.
%
% :Usage:
% ::
%
% t_obj = parcel_stats2statistic_image(atlas_obj, tscores, pvalues, dfe, sig)
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
%   **tscores:**
%        n x k matrix of k t-scores for each of n parcels, where k is,
%        e.g., the number of predictors tested in a multiple regression
%
%   **pvalues:**
%        n x k matrix of k p-values for each of n parcels
%
%   **dfe:**
%        n x 1 matrix of error degrees of freedom for each of n parcels
%
%   **sig:**
%        n x k matrix of k significant tests for each of n parcels
%        This may be, e.g., corrected for multiple comparisons.
%
% :Optional Inputs:
%
% :Outputs:
%
%   **t_obj:**
%        A statistic_image object containing t-maps. .dat contains
%        t-scores, .p contains p-values. Use statistic_image methods like
%        montage(t_obj), orthviews(t_obj), threshold(t_obj) to interact
%        with this object.
%
%
% :Examples:
% ::
%
% :References:
%   None 
%
% :See also:
%   parcel_data2fmri_data
%

% Create placeholder statistic image

atlas_obj = replace_empty(atlas_obj);   
k = size(tscores, 2);
placeholder_vec = ones(atlas_obj.volInfo.n_inmask, k);
all_parcel_idx = double(atlas_obj.dat);
u = unique(all_parcel_idx); u(u == 0) = [];

if size(tscores, 2) ~= size(sig, 2)
    error('Number of images (columns) in tscores and sig do not match')
end

if size(tscores, 2) ~= size(pvalues, 2)
    error('Number of images (columns) in tscores and pvalues do not match')
end

% initialize variables with correct size

t_obj = statistic_image('dat', single(0 .* placeholder_vec), ...
    'p', placeholder_vec, ...
    'sig', logical(placeholder_vec), ...
    'type', 'T', ...
    'dfe', placeholder_vec, ...
    'volInfo', atlas_obj.volInfo);

% number of parcels 

n = size(tscores, 1); % num parcels

if n ~= length(u)
        error('Parcel indices in atlas_obj and extracted parcel-wise t-values do not match. Check code.')
end



t_obj.dat = zeros(size(placeholder_vec, 1), k); % voxels x conditions


    for j = 1:length(u)
        % For each parcel, fill in statistic_image object
        % ----------------------------------------------------
        parcelidx = u(j);
        
        wh_vox = all_parcel_idx == parcelidx;
        
        % map parcels to voxels
        t_obj.dat(wh_vox, :) = repmat(tscores(j, :), sum(wh_vox), 1);
        t_obj.p(wh_vox, :) = repmat(pvalues(j, :), sum(wh_vox), 1);
        t_obj.dfe(wh_vox, :) = repmat(dfe(j, 1), sum(wh_vox), k);       % replicate dfe
        t_obj.sig(wh_vox, :) = repmat(sig(j, :), sum(wh_vox), 1);
                
    end
       
t_obj = enforce_variable_types(t_obj);  % space-saving: 5/24/17

% input_struct.t_statistic_obj = t_obj;

end % function
